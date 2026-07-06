"""Unannotated region discovery and per-region context labelling.

Pure numpy + the Annotation index (no pysam), so this is unit-testable. The BAM
reading layer (bam_io) produces the per-strand coverage arrays that feed here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .config import Config
from .classify import (
    RegionMetrics,
    compute_region_metrics,
    classify_region,
    is_kept,
)
from .gtf_io import Annotation, Gene
from .utils import find_true_runs, merge_runs_with_gap

# Context labels.
CTX_INTERGENIC = "intergenic"
CTX_INTRONIC = "intronic"
CTX_ANTISENSE = "antisense_to_gene"
CTX_SENSE_GENIC = "sense_genic_unannotated"
CTX_NEAR_GENE = "near_gene"


@dataclass
class Candidate:
    """A discovered unannotated interval with metrics and classification."""

    region_id: str
    chrom: str
    start: int          # 0-based
    end: int            # half-open
    metrics: RegionMetrics
    label: str
    reason: str
    flags: dict
    context_label: str
    nearest_feature_id: Optional[str]
    nearest_feature_distance: int
    kept: bool = field(default=False)
    unknown_transcript_name: Optional[str] = None

    @property
    def length(self) -> int:
        return self.end - self.start


def discover_intervals(
    plus_cov: np.ndarray,
    minus_cov: np.ndarray,
    annotated_mask: np.ndarray,
    cfg: Config,
) -> List[tuple]:
    """Find candidate unannotated intervals on one chromosome.

    Steps:
      1. combined depth >= min_depth defines "covered".
      2. remove annotated positions.
      3. merge covered runs across gaps <= max_gap (never bridging annotation).
      4. filter by length, covered bases and covered fraction.
    Returns a list of (start, end) 0-based half-open tuples.
    """
    combined = plus_cov + minus_cov
    covered = combined >= cfg.min_depth
    covered &= ~annotated_mask

    runs = find_true_runs(covered)
    merged = merge_runs_with_gap(runs, max_gap=cfg.max_gap, blocked_mask=annotated_mask)

    out = []
    for start, end in merged:
        length = end - start
        if length < cfg.min_region_length:
            continue
        cov_bases = int(covered[start:end].sum())
        if cov_bases < cfg.min_covered_bases:
            continue
        if length and cov_bases / length < cfg.min_covered_fraction:
            continue
        out.append((start, end))
    return out


def _context_label(
    region_strand: str,
    overlapping: List[Gene],
    exon_mask: np.ndarray,
    start: int,
    end: int,
    nearest_distance: int,
    cfg: Config,
) -> str:
    """Classify genomic context of a region relative to annotation."""
    if overlapping:
        # Inside one or more gene spans (but region itself is unannotated by the
        # chosen mask, e.g. intronic when annotation-mode == exon).
        strands = {g.strand for g in overlapping}
        # Does it overlap an exon at all? (only meaningful if mask wasn't exon)
        region_exonic = bool(exon_mask[start:end].any())
        same_strand = region_strand in strands or region_strand == "."
        opposite = any(
            (region_strand == "+" and g.strand == "-")
            or (region_strand == "-" and g.strand == "+")
            for g in overlapping
        )
        if opposite and not same_strand:
            return CTX_ANTISENSE
        if region_exonic:
            return CTX_SENSE_GENIC
        return CTX_INTRONIC
    if 0 <= nearest_distance <= cfg.nearest_feature_window:
        return CTX_NEAR_GENE
    return CTX_INTERGENIC


def build_candidates_for_chrom(
    chrom: str,
    plus_cov: np.ndarray,
    minus_cov: np.ndarray,
    multi_cov: np.ndarray,
    annotation: Annotation,
    cfg: Config,
    annotated_mask: "np.ndarray | None" = None,
    exon_mask: "np.ndarray | None" = None,
) -> List[Candidate]:
    """Discover intervals on a chromosome and attach metrics + context.

    ``region_id`` values are provisional (per-chrom); the pipeline renumbers
    them globally and assigns unknown_transcript names after gathering all
    chromosomes so ordering is deterministic. Masks may be passed in to avoid
    rebuilding them when the caller already has them.
    """
    length = plus_cov.size
    if annotated_mask is None:
        annotated_mask = annotation.mask_array(chrom, length)
    if exon_mask is None:
        exon_mask = annotation.exon_mask_array(chrom, length)

    intervals = discover_intervals(plus_cov, minus_cov, annotated_mask, cfg)

    candidates: List[Candidate] = []
    for idx, (start, end) in enumerate(intervals):
        p = plus_cov[start:end]
        m = minus_cov[start:end]
        mm = multi_cov[start:end]
        metrics = compute_region_metrics(p, m, cfg, multi_cov=mm)
        label, reason, flags = classify_region(metrics, cfg)

        nearest_gene, distance, overlapping = annotation.nearest_gene(chrom, start, end)
        context = _context_label(
            metrics.dominant_strand, overlapping, exon_mask, start, end, distance, cfg
        )

        candidates.append(
            Candidate(
                region_id=f"{chrom}:{start}-{end}",
                chrom=chrom,
                start=start,
                end=end,
                metrics=metrics,
                label=label,
                reason=reason,
                flags=flags,
                context_label=context,
                nearest_feature_id=nearest_gene.gene_id if nearest_gene else None,
                nearest_feature_distance=distance,
                kept=is_kept(label),
            )
        )
    return candidates
