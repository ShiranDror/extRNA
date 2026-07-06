"""Central configuration objects and defaults for the gDNA-vs-novel pipeline.

Every tunable threshold lives here in one dataclass so the classifier stays
transparent: a reviewer can read the defaults, and the CLI simply overrides
fields on this object. Nothing else in the codebase hard-codes a threshold.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Optional


# Recognised values for library strandedness.
STRANDEDNESS_CHOICES = ("auto", "forward", "reverse", "unstranded")

# What counts as "annotated" when masking coverage.
ANNOTATION_MODES = ("exon", "transcript", "gene", "all")


@dataclass
class Config:
    """All run parameters and thresholds.

    Coordinates are handled internally as 0-based half-open intervals
    (BED convention). Conversion to 1-based inclusive happens only at GTF/TSV
    write time.
    """

    # --- Inputs / outputs -------------------------------------------------
    bam: str = ""
    gtf: str = ""
    out_prefix: str = "sample_analysis"
    fai: Optional[str] = None            # .fai or genome-sizes file (optional)

    # --- Read filtering ---------------------------------------------------
    min_mapq: int = 20                   # STAR uniquely-mapped reads are MAPQ 255,
                                         # multimappers are 3/1/0, so >=20 keeps
                                         # unique reads and drops multimappers.
    min_baseq: int = 0                   # 0 => fast block-based coverage; >0 uses
                                         # per-base filtering (slower).
    keep_duplicates: bool = False        # by default PCR/optical duplicates dropped.

    # --- Multimapper handling --------------------------------------------
    # Reads with MAPQ >= min_mapq (STAR unique = 255) build the coverage that
    # DEFINES regions. Reads below that (STAR multimappers, MAPQ 3/1/0) and
    # secondary alignments are tracked SEPARATELY as a noise signal: a locus
    # swamped by multimapped reads is likely a repeat/alignment artifact, not a
    # real transcript. A genuine region must be established by a majority of
    # uniquely-mapped reads.
    count_secondary: bool = True         # count secondary alignments in the
                                         # multimapper track (repeat pile-ups).
    min_unique_fraction: float = 0.50    # region kept only if uniquely-mapped
                                         # reads are >= this fraction of total
                                         # (unique + multimapped) coverage.

    # --- Region discovery -------------------------------------------------
    min_depth: int = 10                  # min combined (both-strand) per-base depth
                                         # for a base to count as "covered".
    strand_min_depth: int = 3            # min per-strand per-base depth for a base
                                         # to count as "supported" on that strand.
    max_gap: int = 50                    # merge covered runs separated by <= this
                                         # many *unannotated* bases.
    min_region_length: int = 200         # discard candidate intervals shorter than this.
    min_covered_bases: int = 100         # discard regions with fewer covered bases.
    min_covered_fraction: float = 0.7    # discard regions whose covered fraction is
                                         # below this (sparse/punctate signal).

    # --- Annotation -------------------------------------------------------
    annotation_mode: str = "exon"        # exon | transcript | gene | all
    nearest_feature_window: int = 10000  # window used when reporting nearest feature
                                         # / deciding intergenic vs near-gene context.
    stranded_masking: bool = True        # mask annotation per-strand (stranded libs
                                         # only) so antisense-over-feature signal
                                         # stays discoverable; auto-falls back to
                                         # positional masking for unstranded libs.

    # --- Library strandedness --------------------------------------------
    library_strandedness: str = "auto"   # auto | forward | reverse | unstranded
    strand_infer_max_reads: int = 200000 # reads sampled for auto inference.
    strand_infer_min_confidence: float = 0.80  # fraction needed to call fwd/rev.

    # --- gDNA classification thresholds ----------------------------------
    gdna_min_dual_strand_fraction: float = 0.60   # fraction of interval covered on
                                                  # BOTH strands.
    gdna_max_strand_length_ratio_diff: float = 0.25  # |plus_len-minus_len|/max <= this.
    gdna_min_profile_correlation: float = 0.70    # Pearson r of per-base +/- profiles.
    gdna_min_covered_fraction: float = 0.50       # gDNA should be broadly continuous.
    gdna_min_depth_balance: float = 0.50          # min(plus,minus)/max mean depth;
                                                  # gDNA has balanced strand depth.
    gdna_flat_cv_threshold: float = 0.40          # if BOTH strands' coverage is this
                                                  # flat (coeff. of variation), the
                                                  # uniform broad profile is itself a
                                                  # gDNA signature and correlation is
                                                  # not required (randomly-fragmented
                                                  # gDNA gives flat, uncorrelated cov).

    # --- Novel-transcript threshold --------------------------------------
    novel_min_dominant_strand_fraction: float = 0.80  # one strand carries >= this
                                                      # fraction of signal => novel.

    # --- Performance / misc ----------------------------------------------
    threads: int = 4                     # chromosome-level parallelism.
    emit_bed: bool = True                # write candidate_regions.bed
    emit_bedgraph: bool = False          # write per-strand bedGraph of candidates.
    emit_multiqc: bool = True            # write *.gdna_mqc.tsv for MultiQC.
    sample_name: Optional[str] = None    # label used in the MultiQC row (defaults
                                         # to the out-prefix basename).
    verbose: bool = False

    # Populated at run time (not user-set):
    inferred_strandedness: Optional[str] = None
    strandedness_metrics: dict = field(default_factory=dict)

    def to_serialisable(self) -> dict:
        """Return a JSON-friendly dict of the run parameters."""
        return asdict(self)

    def validate(self) -> None:
        """Sanity-check parameter combinations, raising ValueError on problems."""
        if self.library_strandedness not in STRANDEDNESS_CHOICES:
            raise ValueError(
                f"--library-strandedness must be one of {STRANDEDNESS_CHOICES}"
            )
        if self.annotation_mode not in ANNOTATION_MODES:
            raise ValueError(f"--annotation-mode must be one of {ANNOTATION_MODES}")
        if self.min_depth < 1:
            raise ValueError("--min-depth must be >= 1")
        if self.strand_min_depth < 1:
            raise ValueError("--strand-min-depth must be >= 1")
        if not (0.0 <= self.min_covered_fraction <= 1.0):
            raise ValueError("--min-covered-fraction must be in [0, 1]")
        if not (0.0 <= self.gdna_min_dual_strand_fraction <= 1.0):
            raise ValueError("--gdna-min-dual-strand-fraction must be in [0, 1]")
        if not (0.0 <= self.gdna_max_strand_length_ratio_diff <= 1.0):
            raise ValueError("--gdna-max-strand-length-ratio-diff must be in [0, 1]")
        if not (-1.0 <= self.gdna_min_profile_correlation <= 1.0):
            raise ValueError("--gdna-min-profile-correlation must be in [-1, 1]")
        if not (0.0 <= self.gdna_min_depth_balance <= 1.0):
            raise ValueError("--gdna-min-depth-balance must be in [0, 1]")
        if self.gdna_flat_cv_threshold < 0.0:
            raise ValueError("--gdna-flat-cv-threshold must be >= 0")
        if not (0.5 <= self.novel_min_dominant_strand_fraction <= 1.0):
            raise ValueError(
                "--novel-min-dominant-strand-fraction must be in [0.5, 1]"
            )
        if not (0.0 <= self.min_unique_fraction <= 1.0):
            raise ValueError("--min-unique-fraction must be in [0, 1]")
        if self.threads < 1:
            raise ValueError("--threads must be >= 1")
