"""Automatic library-strandedness inference.

We sample reads overlapping annotated exons whose gene sits on a single strand,
and ask how often the read's *forward-hypothesis* transcription strand matches
the annotated gene strand. This mirrors RSeQC's `infer_experiment.py` logic:

  * fraction ~0.5            -> unstranded
  * fraction high (>= thr)   -> forward (fr-secondstrand)
  * fraction low  (<= 1-thr) -> reverse (fr-firststrand / dUTP)

The result and its confidence are reported in the summary JSON.
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

from .config import Config
from .gtf_io import Annotation
from .utils import get_logger

# random.Random is seeded deterministically so a run is reproducible.
_RNG = random.Random(12345)


def _sample_exon_intervals(
    annotation: Annotation, max_intervals: int
) -> List[Tuple[str, int, int, str]]:
    """Pick exons that unambiguously belong to one strand at their locus."""
    candidates: List[Tuple[str, int, int, str]] = []
    for chrom, exons in annotation.stranded_exons.items():
        for start, end, strand in exons:
            if strand in ("+", "-") and (end - start) >= 50:
                candidates.append((chrom, start, end, strand))
    if len(candidates) > max_intervals:
        candidates = _RNG.sample(candidates, max_intervals)
    return candidates


def infer_strandedness(
    bam_path: str, annotation: Annotation, cfg: Config
) -> Dict:
    """Return a metrics dict including the inferred library strandedness.

    Keys: decision, p_forward, p_reverse, reads_sampled, confidence, note.
    """
    from .bam_io import iter_reads_over_intervals, transcription_strand

    logger = get_logger(cfg.verbose)

    # Sample a subset of exons; cap read count via early break.
    exon_ivals = _sample_exon_intervals(annotation, max_intervals=5000)
    # Build a strand lookup keyed by (chrom, start, end) for the sampled set.
    strand_by_ival = {(c, s, e): st for (c, s, e, st) in exon_ivals}
    intervals = [(c, s, e) for (c, s, e, _) in exon_ivals]

    forward_votes = 0
    total = 0
    # We re-fetch per interval so we can attribute each read to that exon's strand.
    from .bam_io import open_bam, _passes_filters

    bam = open_bam(bam_path)
    for chrom, start, end, gene_strand in exon_ivals:
        try:
            reads = bam.fetch(chrom, start, end)
        except ValueError:
            continue
        for read in reads:
            if not _passes_filters(read, cfg):
                continue
            fwd_strand = transcription_strand(read, "forward")
            if fwd_strand == gene_strand:
                forward_votes += 1
            total += 1
            if total >= cfg.strand_infer_max_reads:
                break
        if total >= cfg.strand_infer_max_reads:
            break
    bam.close()

    if total == 0:
        logger.warning(
            "Strandedness inference found no reads over annotated exons; "
            "defaulting to unstranded."
        )
        return {
            "decision": "unstranded",
            "p_forward": None,
            "p_reverse": None,
            "reads_sampled": 0,
            "confidence": 0.0,
            "note": "no reads sampled over exons; defaulted to unstranded",
        }

    p_forward = forward_votes / total
    p_reverse = 1.0 - p_forward
    thr = cfg.strand_infer_min_confidence

    if p_forward >= thr:
        decision = "forward"
        confidence = p_forward
    elif p_reverse >= thr:
        decision = "reverse"
        confidence = p_reverse
    else:
        decision = "unstranded"
        confidence = 1.0 - abs(p_forward - 0.5) * 2  # ~1 when perfectly balanced

    note = (
        f"p_forward={p_forward:.3f} over {total} exon-overlapping reads "
        f"(threshold {thr})"
    )
    logger.info("Inferred strandedness: %s (%s)", decision, note)
    return {
        "decision": decision,
        "p_forward": p_forward,
        "p_reverse": p_reverse,
        "reads_sampled": total,
        "confidence": confidence,
        "note": note,
    }


def resolve_strandedness(
    bam_path: str, annotation: Annotation, cfg: Config
) -> Tuple[str, Dict]:
    """Resolve the effective strandedness, honouring an explicit user choice."""
    logger = get_logger(cfg.verbose)
    if cfg.library_strandedness != "auto":
        metrics = {
            "decision": cfg.library_strandedness,
            "note": "set explicitly by user (--library-strandedness)",
            "confidence": 1.0,
            "reads_sampled": 0,
        }
        return cfg.library_strandedness, metrics

    metrics = infer_strandedness(bam_path, annotation, cfg)
    if metrics["decision"] == "unstranded":
        logger.warning(
            "Library appears UNSTRANDED. gDNA detection relies on strand "
            "symmetry, which is unreliable without strand information: even "
            "genuine single-strand RNA maps to both strands in an unstranded "
            "library. Interpret likely_gDNA calls with caution."
        )
    return metrics["decision"], metrics
