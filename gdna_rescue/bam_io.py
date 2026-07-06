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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (plus_unique, minus_unique, multi) int32 coverage arrays.

    ``plus_unique`` / ``minus_unique`` are strand-specific coverage from
    uniquely-mapped reads (these DEFINE candidate regions). ``multi`` is
    strand-agnostic coverage from multimapped/secondary reads, used only to
    estimate the uniquely-mapped fraction of a region.

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
            _add_blocks(multi, read)
            continue
        strand = transcription_strand(read, strandedness)
        _add_blocks(plus if strand == "+" else minus, read)

    bam.close()

    if use_baseq:
        return plus, minus, multi
    return (
        np.cumsum(plus[:-1]).astype(np.int32),
        np.cumsum(minus[:-1]).astype(np.int32),
        np.cumsum(multi[:-1]).astype(np.int32),
    )


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
