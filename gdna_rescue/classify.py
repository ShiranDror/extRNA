"""Region-level metrics and the transparent rule-based classifier.

This module is intentionally free of pysam / I/O so that the biological logic
can be unit-tested with plain numpy arrays and swapped for a statistical model
later. The public entry points are:

    compute_region_metrics(plus_cov, minus_cov, cfg) -> RegionMetrics
    classify_region(metrics, cfg)                     -> (label, reason, flags)

``plus_cov`` / ``minus_cov`` are 1-D numpy arrays of per-base transcription-strand
coverage over exactly the candidate interval (same length).

Biological rationale (see README for the long form):
  * Genuine RNA — even continuous / unspliced / intronic / antisense — is
    transcribed from ONE template strand, so in a stranded library the reads
    resolve predominantly to a single transcription strand.
  * Genomic-DNA contamination derives from double-stranded DNA, so fragments
    map to both strands roughly symmetrically and their per-base coverage
    profiles are correlated across the same span.
  * Therefore the discriminating signal is strand *symmetry + profile
    correlation over a broad continuous span*, NOT the mere presence of two
    strands (real loci can be bidirectional) and NOT splice junctions.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, Tuple

import numpy as np

from .config import Config
from .utils import safe_pearson, binary_entropy

# Classification labels.
LIKELY_GDNA = "likely_gDNA"
POSSIBLE_BIDIRECTIONAL = "possible_bidirectional_RNA"
LIKELY_NOVEL = "likely_novel_transcript"
LIKELY_MULTIMAPPER = "likely_multimapper_artifact"

_EPS = 1e-9


@dataclass
class RegionMetrics:
    """All per-region quantities used by the classifier and the TSV report."""

    length: int
    total_coverage: float
    covered_bases: int
    covered_fraction: float
    avg_depth: float
    max_depth: float

    plus_covered_len: int
    minus_covered_len: int
    plus_mean_depth: float
    minus_mean_depth: float
    plus_total: float
    minus_total: float

    unique_total: float              # total coverage from uniquely-mapped reads
    multi_total: float               # total coverage from multimapped reads
    unique_fraction: float           # unique / (unique + multi), 0..1

    pm_depth_ratio: float            # min/max of strand mean depths (0..1)
    covered_len_ratio: float         # min/max of strand covered lengths (0..1)
    strand_length_ratio_diff: float  # |plus_len-minus_len|/max (0..1)
    dominant_strand_fraction: float  # max strand total / total (0.5..1)
    strand_entropy: float            # 1 = balanced, 0 = single strand
    strand_overlap_jaccard: float    # |plus & minus| / |plus | minus|
    dual_strand_fraction: float      # |plus & minus| / length
    profile_correlation: float       # Pearson r of +/- per-base profiles
    plus_cv: float                   # coeff. of variation of + coverage
    minus_cv: float                  # coeff. of variation of - coverage
    dominant_strand: str             # '+', '-' or '.'

    def to_dict(self) -> Dict:
        return asdict(self)


def compute_region_metrics(
    plus_cov: np.ndarray,
    minus_cov: np.ndarray,
    cfg: Config,
    multi_cov: np.ndarray | None = None,
) -> RegionMetrics:
    """Compute all region metrics from per-base coverage arrays.

    ``plus_cov`` / ``minus_cov`` are uniquely-mapped strand coverage.
    ``multi_cov`` is multimapped coverage over the same interval (defaults to
    zeros, e.g. in unit tests of the pure strand logic).
    """
    plus_cov = np.asarray(plus_cov, dtype=np.float64)
    minus_cov = np.asarray(minus_cov, dtype=np.float64)
    assert plus_cov.shape == minus_cov.shape, "strand arrays must be same length"
    if multi_cov is None:
        multi_cov = np.zeros_like(plus_cov)
    else:
        multi_cov = np.asarray(multi_cov, dtype=np.float64)

    length = int(plus_cov.size)
    combined = plus_cov + minus_cov

    covered_mask = combined >= cfg.min_depth
    covered_bases = int(covered_mask.sum())
    covered_fraction = covered_bases / length if length else 0.0

    total_coverage = float(combined.sum())
    avg_depth = total_coverage / length if length else 0.0
    max_depth = float(combined.max()) if length else 0.0

    plus_sup = plus_cov >= cfg.strand_min_depth
    minus_sup = minus_cov >= cfg.strand_min_depth
    plus_covered_len = int(plus_sup.sum())
    minus_covered_len = int(minus_sup.sum())

    plus_total = float(plus_cov.sum())
    minus_total = float(minus_cov.sum())
    # Mean depth taken over the whole interval (includes zeros) so it reflects
    # both intensity and breadth on that strand.
    plus_mean_depth = plus_total / length if length else 0.0
    minus_mean_depth = minus_total / length if length else 0.0

    def _ratio(a: float, b: float) -> float:
        hi = max(a, b)
        return (min(a, b) / hi) if hi > 0 else 0.0

    pm_depth_ratio = _ratio(plus_mean_depth, minus_mean_depth)
    covered_len_ratio = _ratio(plus_covered_len, minus_covered_len)
    hi_len = max(plus_covered_len, minus_covered_len)
    strand_length_ratio_diff = (
        abs(plus_covered_len - minus_covered_len) / hi_len if hi_len > 0 else 0.0
    )

    strand_total = plus_total + minus_total
    dominant_strand_fraction = (
        max(plus_total, minus_total) / strand_total if strand_total > 0 else 1.0
    )
    p_plus = plus_total / strand_total if strand_total > 0 else 0.5
    strand_entropy = binary_entropy(p_plus)

    dual_mask = plus_sup & minus_sup
    union_mask = plus_sup | minus_sup
    dual_bases = int(dual_mask.sum())
    union_bases = int(union_mask.sum())
    strand_overlap_jaccard = dual_bases / union_bases if union_bases > 0 else 0.0
    dual_strand_fraction = dual_bases / length if length else 0.0

    profile_correlation = safe_pearson(plus_cov, minus_cov)

    unique_total = plus_total + minus_total
    multi_total = float(multi_cov.sum())
    denom = unique_total + multi_total
    unique_fraction = (unique_total / denom) if denom > 0 else 1.0

    def _cv(arr: np.ndarray) -> float:
        mean = arr.mean()
        return float(arr.std() / mean) if mean > 0 else 0.0

    plus_cv = _cv(plus_cov)
    minus_cv = _cv(minus_cov)

    if plus_total > minus_total * 1.0 and plus_total > 0 and plus_total != minus_total:
        dominant_strand = "+"
    elif minus_total > plus_total:
        dominant_strand = "-"
    else:
        dominant_strand = "."

    return RegionMetrics(
        length=length,
        total_coverage=total_coverage,
        covered_bases=covered_bases,
        covered_fraction=covered_fraction,
        avg_depth=avg_depth,
        max_depth=max_depth,
        plus_covered_len=plus_covered_len,
        minus_covered_len=minus_covered_len,
        plus_mean_depth=plus_mean_depth,
        minus_mean_depth=minus_mean_depth,
        plus_total=plus_total,
        minus_total=minus_total,
        unique_total=unique_total,
        multi_total=multi_total,
        unique_fraction=unique_fraction,
        pm_depth_ratio=pm_depth_ratio,
        covered_len_ratio=covered_len_ratio,
        strand_length_ratio_diff=strand_length_ratio_diff,
        dominant_strand_fraction=dominant_strand_fraction,
        strand_entropy=strand_entropy,
        strand_overlap_jaccard=strand_overlap_jaccard,
        dual_strand_fraction=dual_strand_fraction,
        profile_correlation=profile_correlation,
        plus_cv=plus_cv,
        minus_cv=minus_cv,
        dominant_strand=dominant_strand,
    )


def classify_region(
    m: RegionMetrics, cfg: Config
) -> Tuple[str, str, Dict[str, bool]]:
    """Apply the transparent rule set. Returns (label, reason, flags).

    Decision order:
      1. Clear single-strand dominance  -> likely_novel_transcript
      2. Both strands + symmetric + correlated + broad -> likely_gDNA
      3. Otherwise (both strands, but asymmetric/uncorrelated) ->
         possible_bidirectional_RNA
    """
    # --- rule 0: multimapper artifact (checked first) ---------------------
    # A genuine region must be established by a majority of uniquely-mapped
    # reads. If multimapped reads dominate, the locus is almost certainly a
    # repeat / alignment artifact and must NOT be rescued, regardless of its
    # strand pattern. (Local multimapper-only stretches are fine as long as the
    # region as a whole is majority-unique.)
    unique_majority = m.unique_fraction >= cfg.min_unique_fraction
    if not unique_majority:
        reason = (
            f"uniquely-mapped fraction {m.unique_fraction:.2f} < "
            f"{cfg.min_unique_fraction:.2f} (multi coverage {m.multi_total:.0f} vs "
            f"unique {m.unique_total:.0f}) -> multimapper/repeat artifact."
        )
        return LIKELY_MULTIMAPPER, reason, {"unique_reads_majority": False}

    # --- individual boolean tests (reported for transparency) --------------
    single_strand = m.dominant_strand_fraction >= cfg.novel_min_dominant_strand_fraction

    gdna_dual = m.dual_strand_fraction >= cfg.gdna_min_dual_strand_fraction
    gdna_symmetric = m.strand_length_ratio_diff <= cfg.gdna_max_strand_length_ratio_diff
    gdna_broad = m.covered_fraction >= cfg.gdna_min_covered_fraction
    gdna_balanced = m.pm_depth_ratio >= cfg.gdna_min_depth_balance

    # gDNA can appear two ways over a broad dual-strand span:
    #   (a) correlated per-base profiles (shared mappability structure), or
    #   (b) FLAT coverage on both strands (randomly fragmented gDNA is uniform,
    #       so profiles have little variance and do not correlate). Requiring
    #       correlation alone would miss the very common flat case, so a flat
    #       symmetric profile counts as consistent too.
    gdna_correlated = m.profile_correlation >= cfg.gdna_min_profile_correlation
    both_flat = (m.plus_cv <= cfg.gdna_flat_cv_threshold
                 and m.minus_cv <= cfg.gdna_flat_cv_threshold)
    gdna_pattern = gdna_correlated or both_flat

    flags = {
        "unique_reads_majority": unique_majority,
        "single_strand_dominant": single_strand,
        "gdna_dual_strand": gdna_dual,
        "gdna_symmetric_length": gdna_symmetric,
        "gdna_broad_continuous": gdna_broad,
        "gdna_depth_balanced": gdna_balanced,
        "gdna_profile_correlated": gdna_correlated,
        "gdna_both_strands_flat": both_flat,
        "gdna_pattern_consistent": gdna_pattern,
    }

    # --- rule 1: single dominant strand -> genuine (novel) transcript ------
    if single_strand:
        reason = (
            f"dominant strand carries {m.dominant_strand_fraction:.2f} of signal "
            f">= {cfg.novel_min_dominant_strand_fraction:.2f}; "
            f"dual-strand fraction {m.dual_strand_fraction:.2f}, "
            f"profile r={m.profile_correlation:.2f} -> single-strand transcription."
        )
        return LIKELY_NOVEL, reason, flags

    # --- rule 2: symmetric, balanced, broad dual strand + consistent pattern
    if gdna_dual and gdna_symmetric and gdna_broad and gdna_balanced and gdna_pattern:
        pattern_txt = (
            f"profile r={m.profile_correlation:.2f} >= "
            f"{cfg.gdna_min_profile_correlation:.2f}"
            if gdna_correlated
            else f"both strands flat (CV {m.plus_cv:.2f}/{m.minus_cv:.2f} <= "
                 f"{cfg.gdna_flat_cv_threshold:.2f})"
        )
        reason = (
            f"dual-strand fraction {m.dual_strand_fraction:.2f} "
            f">= {cfg.gdna_min_dual_strand_fraction:.2f}, "
            f"strand length diff {m.strand_length_ratio_diff:.2f} "
            f"<= {cfg.gdna_max_strand_length_ratio_diff:.2f}, "
            f"depth balance {m.pm_depth_ratio:.2f} >= "
            f"{cfg.gdna_min_depth_balance:.2f}, "
            f"covered fraction {m.covered_fraction:.2f} "
            f">= {cfg.gdna_min_covered_fraction:.2f}, {pattern_txt} "
            f"-> symmetric dual-strand genomic pattern."
        )
        return LIKELY_GDNA, reason, flags

    # --- rule 3: both strands present but not gDNA-like --------------------
    failed = []
    if not gdna_dual:
        failed.append(
            f"dual fraction {m.dual_strand_fraction:.2f} < "
            f"{cfg.gdna_min_dual_strand_fraction:.2f}"
        )
    if not gdna_symmetric:
        failed.append(
            f"strand length diff {m.strand_length_ratio_diff:.2f} > "
            f"{cfg.gdna_max_strand_length_ratio_diff:.2f}"
        )
    if not gdna_balanced:
        failed.append(
            f"depth balance {m.pm_depth_ratio:.2f} < "
            f"{cfg.gdna_min_depth_balance:.2f}"
        )
    if not gdna_pattern:
        failed.append(
            f"neither correlated (r={m.profile_correlation:.2f} < "
            f"{cfg.gdna_min_profile_correlation:.2f}) nor flat on both strands "
            f"(CV {m.plus_cv:.2f}/{m.minus_cv:.2f} > "
            f"{cfg.gdna_flat_cv_threshold:.2f})"
        )
    if not gdna_broad:
        failed.append(
            f"covered fraction {m.covered_fraction:.2f} < "
            f"{cfg.gdna_min_covered_fraction:.2f}"
        )
    reason = (
        "both strands contribute but the pattern is not gDNA-symmetric ("
        + "; ".join(failed)
        + f"); dominant-strand fraction {m.dominant_strand_fraction:.2f} "
        f"< {cfg.novel_min_dominant_strand_fraction:.2f} -> asymmetric bidirectional."
    )
    return POSSIBLE_BIDIRECTIONAL, reason, flags


def is_kept(label: str) -> bool:
    """A region is rescued unless it is gDNA or a multimapper artifact."""
    return label not in (LIKELY_GDNA, LIKELY_MULTIMAPPER)
