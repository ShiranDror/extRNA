"""BAM reading, filtering and per-strand coverage extraction.

This is the ONLY module that imports pysam. Keeping it isolated means the
discovery / classification logic (and the test suite) run on any platform,
including Windows where pysam has no wheels. Install pysam via bioconda on
Linux / macOS / WSL to run the full pipeline.

STAR-specific notes:
  * STAR sets MAPQ 255 for uniquely-mapped reads and 3/1/0 for multimappers,
    so the default --min-mapq 20 keeps unique reads and drops multimappers.
  * Spliced reads carry N (ref-skip) CIGAR ops; get_blocks() splits on these,
    so introns are never counted as covered by a junction-spanning read while
    genuine continuous/unspliced RNA still accumulates full coverage.
"""

from __future__ import annotations

from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

try:  # pragma: no cover - import guard exercised only where pysam is present
    import pysam
except ImportError:  # pragma: no cover
    pysam = None

from .config import Config
from .utils import flip_strand, get_logger


def _require_pysam() -> None:
    if pysam is None:
        raise ImportError(
            "pysam is required for BAM processing but is not installed. "
            "Install it with `conda install -c bioconda pysam` (Linux/macOS/WSL). "
            "The discovery/classification modules and tests do not need pysam."
        )


def transcription_strand(read, strandedness: str) -> str:
    """Return the transcription strand ('+'/'-') a read supports.

    ``strandedness`` is one of forward | reverse | unstranded. For unstranded
    libraries we fall back to the read's *alignment* strand and the caller is
    warned that gDNA strand-symmetry detection is unreliable in that case.
    """
    aligned = "-" if read.is_reverse else "+"
    if strandedness == "unstranded":
        return aligned

    is_first = read.is_read1 or (not read.is_paired)
    if strandedness == "forward":
        # fr-secondstrand: read1 same strand as transcript.
        return aligned if is_first else flip_strand(aligned)
    # reverse (fr-firststrand / dUTP): read1 opposite strand of transcript.
    return flip_strand(aligned) if is_first else aligned


def read_category(read, cfg: Config):
    """Classify a read as 'unique', 'multi', or None (dropped).

    * unique : primary alignment with MAPQ >= min_mapq (STAR unique = 255).
    * multi  : mapped but a multimapper — primary alignment with MAPQ < min_mapq
               (STAR 3/1/0), or a secondary alignment (if count_secondary).
    * None   : unmapped / qcfail / duplicate / supplementary / secondary when
               secondary counting is disabled.
    """
    if read.is_unmapped or read.is_qcfail or read.is_supplementary:
        return None
    if read.is_duplicate and not cfg.keep_duplicates:
        return None
    if read.is_secondary:
        return "multi" if cfg.count_secondary else None
    # primary alignment
    if read.mapping_quality >= cfg.min_mapq:
        return "unique"
    return "multi"


def _passes_filters(read, cfg: Config) -> bool:
    """True only for uniquely-mapped reads (used by strandedness inference)."""
    return read_category(read, cfg) == "unique"


def open_bam(path: str):
    """Open (and validate) a coordinate-sorted, indexed BAM."""
    _require_pysam()
    bam = pysam.AlignmentFile(path, "rb")
    try:
        if not bam.has_index():
            raise ValueError(
                f"BAM {path!r} has no index. Run `samtools index {path}` first."
            )
    except (ValueError, AttributeError):
        raise ValueError(
            f"BAM {path!r} is missing an index (.bai/.csi). "
            f"Run `samtools index {path}`."
        )
    hdr = bam.header.to_dict()
    so = hdr.get("HD", {}).get("SO")
    if so and so != "coordinate":
        raise ValueError(
            f"BAM {path!r} sort order is {so!r}; a coordinate-sorted BAM is "
            f"required. Run `samtools sort`."
        )
    return bam


def get_chrom_sizes_from_bam(path: str) -> Dict[str, int]:
    """Return {chrom: length} from the BAM header."""
    bam = open_bam(path)
    sizes = dict(zip(bam.references, bam.lengths))
    bam.close()
    return sizes


def strand_coverage_for_chrom(
    path: str,
    chrom: str,
    length: int,
    cfg: Config,
    strandedness: str,
    annotated_mask_plus: "np.ndarray | None" = None,
    annotated_mask_minus: "np.ndarray | None" = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Return (plus_unique, minus_unique, multi, read_stats).

    ``plus_unique`` / ``minus_unique`` are strand-specific coverage from
    uniquely-mapped reads (these DEFINE candidate regions). ``multi`` is
    strand-agnostic coverage from multimapped/secondary reads.

    ``read_stats`` counts reads (not coverage) for QC: total unique reads, unique
    reads whose midpoint falls in an annotated feature ON THAT READ'S
    transcription strand, and multimapped reads. Passing the same positional mask
    as both plus/minus reproduces strand-agnostic counting. Counting by read
    midpoint keeps this O(1) per read within the existing pass.

    Uses difference-array accumulation over aligned blocks (fast path) or a
    per-base loop when base-quality filtering is requested.
    """
    _require_pysam()
    bam = open_bam(path)

    use_baseq = cfg.min_baseq > 0
    size = length if use_baseq else length + 1
    plus = np.zeros(size, dtype=np.int32)
    minus = np.zeros(size, dtype=np.int32)
    multi = np.zeros(size, dtype=np.int32)

    n_unique = 0
    n_unique_annotated = 0
    n_multi = 0

    def _add_blocks(arr, read):
        if use_baseq:
            quals = read.query_qualities
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos is None or rpos >= length:
                    continue
                if quals is not None and quals[qpos] < cfg.min_baseq:
                    continue
                arr[rpos] += 1
        else:
            for bstart, bend in read.get_blocks():
                if bstart >= length:
                    continue
                bend = min(bend, length)
                if bend <= bstart:
                    continue
                arr[bstart] += 1
                arr[bend] -= 1

    for read in bam.fetch(chrom):
        category = read_category(read, cfg)
        if category is None:
            continue
        if category == "multi":
            n_multi += 1
            _add_blocks(multi, read)
            continue
        # unique read
        n_unique += 1
        strand = transcription_strand(read, strandedness)
        strand_mask = annotated_mask_plus if strand == "+" else annotated_mask_minus
        if strand_mask is not None:
            rs = read.reference_start
            re = read.reference_end or (rs + 1)
            mid = (rs + re) // 2
            if 0 <= mid < length and strand_mask[mid]:
                n_unique_annotated += 1
        _add_blocks(plus if strand == "+" else minus, read)

    bam.close()

    read_stats = {
        "n_unique_reads": n_unique,
        "n_unique_reads_annotated": n_unique_annotated,
        "n_multi_reads": n_multi,
    }

    if use_baseq:
        return plus, minus, multi, read_stats
    return (
        np.cumsum(plus[:-1]).astype(np.int32),
        np.cumsum(minus[:-1]).astype(np.int32),
        np.cumsum(multi[:-1]).astype(np.int32),
        read_stats,
    )


def count_unique_reads_in_intervals(
    path: str,
    intervals: List[Tuple[str, int, int]],
    cfg: Config,
    strandedness: str = "unstranded",
    annotation=None,
    stranded_masking: bool = False,
) -> List[int]:
    """Count uniquely-mapped reads assigned to each candidate region.

    ``intervals`` are (chrom, start, end). One targeted fetch per interval
    (candidate regions are few, so this avoids a second full pass). A read is
    assigned by its midpoint, so it counts toward at most one region.

    A read is EXCLUDED from a region if it is already counted as annotated on its
    own transcription strand (same rule used for the annotated read tally), so
    the categories partition the reads without double counting: e.g. an antisense
    region over an exon does not also count the host gene's sense reads (those are
    annotated). Both strands of a genuine bidirectional / gDNA region are counted,
    because those positions are unannotated.
    """
    _require_pysam()
    bam = open_bam(path)

    # Lazily build per-chrom sorted mask starts/ends for a point-in-mask test.
    cache: dict = {}

    def _mask_arrays(chrom):
        if chrom not in cache:
            if stranded_masking and annotation is not None:
                pv = annotation.mask_intervals_plus.get(chrom, [])
                mv = annotation.mask_intervals_minus.get(chrom, [])
            elif annotation is not None:
                pos = annotation.mask_intervals.get(chrom, [])
                pv = mv = pos
            else:
                pv = mv = []
            cache[chrom] = {
                "+": (np.array([s for s, _ in pv]), np.array([e for _, e in pv])),
                "-": (np.array([s for s, _ in mv]), np.array([e for _, e in mv])),
            }
        return cache[chrom]

    def _annotated_on_strand(chrom, pos, strand):
        if annotation is None:
            return False
        starts, ends = _mask_arrays(chrom)[strand if strand in ("+", "-") else "+"]
        if starts.size == 0:
            return False
        i = int(np.searchsorted(starts, pos, side="right")) - 1
        return i >= 0 and pos < ends[i]

    counts: List[int] = []
    for chrom, start, end in intervals:
        n = 0
        try:
            for read in bam.fetch(chrom, start, end):
                if read_category(read, cfg) != "unique":
                    continue
                rs = read.reference_start
                re = read.reference_end or (rs + 1)
                mid = (rs + re) // 2
                if not (start <= mid < end):
                    continue
                strand = transcription_strand(read, strandedness)
                if _annotated_on_strand(chrom, mid, strand):
                    continue  # already counted as annotated on its own strand
                n += 1
        except ValueError:
            pass
        counts.append(n)
    bam.close()
    return counts


def iter_reads_over_intervals(
    path: str, intervals: List[Tuple[str, int, int]], cfg: Config
) -> Iterator:
    """Yield filtered reads overlapping the given (chrom, start, end) intervals.

    Used by strandedness inference to sample reads over annotated exons.
    """
    _require_pysam()
    bam = open_bam(path)
    for chrom, start, end in intervals:
        try:
            for read in bam.fetch(chrom, start, end):
                if _passes_filters(read, cfg):
                    yield read
        except ValueError:
            # Chromosome absent from BAM header; skip.
            continue
    bam.close()
