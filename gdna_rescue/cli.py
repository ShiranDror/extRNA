"""Command-line interface for the gDNA-vs-novel-transcript tool."""

from __future__ import annotations

import argparse
import sys

from .config import Config, STRANDEDNESS_CHOICES, ANNOTATION_MODES
from .utils import get_logger


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="detect_gdna_vs_novel.py",
        description=(
            "Detect likely genomic-DNA contamination and rescue candidate "
            "unannotated transcripts from RNA-seq alignments, using strand "
            "symmetry and continuous-coverage evidence (not splice junctions)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    req = p.add_argument_group("required inputs")
    req.add_argument("--bam", required=True, help="Coordinate-sorted, indexed BAM.")
    req.add_argument("--gtf", required=True, help="Reference GTF annotation.")
    req.add_argument("--out-prefix", required=True, help="Output file prefix.")

    opt = p.add_argument_group("optional inputs")
    opt.add_argument("--fai", default=None,
                     help="FASTA .fai index or 2-column genome-sizes file for "
                          "coordinate sanity checks.")

    filt = p.add_argument_group("read filtering")
    filt.add_argument("--min-mapq", type=int, default=20)
    filt.add_argument("--min-baseq", type=int, default=0)
    filt.add_argument("--keep-duplicates", action="store_true",
                      help="Keep reads flagged as duplicates (default: drop).")
    filt.add_argument("--min-unique-fraction", type=float, default=0.50,
                      help="Region kept only if uniquely-mapped reads are >= this "
                           "fraction of total coverage; else flagged as a "
                           "multimapper/repeat artifact.")
    filt.add_argument("--no-count-secondary", dest="count_secondary",
                      action="store_false",
                      help="Do not count secondary alignments in the multimapper "
                           "track.")

    disc = p.add_argument_group("region discovery")
    disc.add_argument("--min-depth", type=int, default=10,
                      help="Min combined (both-strand) per-base depth for a base "
                           "to count as covered.")
    disc.add_argument("--strand-min-depth", type=int, default=3,
                      help="Per-strand per-base depth for a base to count as supported.")
    disc.add_argument("--max-gap", type=int, default=50)
    disc.add_argument("--min-region-length", type=int, default=200)
    disc.add_argument("--min-covered-bases", type=int, default=100)
    disc.add_argument("--min-covered-fraction", type=float, default=0.7)

    ann = p.add_argument_group("annotation")
    ann.add_argument("--annotation-mode", choices=ANNOTATION_MODES, default="exon")
    ann.add_argument("--nearest-feature-window", type=int, default=10000)

    strand = p.add_argument_group("library strandedness")
    strand.add_argument("--library-strandedness", choices=STRANDEDNESS_CHOICES,
                        default="auto")
    strand.add_argument("--strand-infer-max-reads", type=int, default=200000)
    strand.add_argument("--strand-infer-min-confidence", type=float, default=0.80)

    gd = p.add_argument_group("gDNA classification thresholds")
    gd.add_argument("--gdna-min-dual-strand-fraction", type=float, default=0.60)
    gd.add_argument("--gdna-max-strand-length-ratio-diff", type=float, default=0.25)
    gd.add_argument("--gdna-min-profile-correlation", type=float, default=0.70)
    gd.add_argument("--gdna-min-covered-fraction", type=float, default=0.50)
    gd.add_argument("--gdna-min-depth-balance", type=float, default=0.50)
    gd.add_argument("--gdna-flat-cv-threshold", type=float, default=0.40)

    nv = p.add_argument_group("novel-transcript threshold")
    nv.add_argument("--novel-min-dominant-strand-fraction", type=float, default=0.80)

    misc = p.add_argument_group("performance / outputs")
    misc.add_argument("--threads", type=int, default=4)
    misc.add_argument("--no-bed", dest="emit_bed", action="store_false",
                      help="Do not write candidate_regions.bed.")
    misc.add_argument("--emit-bedgraph", action="store_true",
                      help="Write per-strand candidate bedGraph files.")
    misc.add_argument("--verbose", action="store_true")
    return p


def args_to_config(args: argparse.Namespace) -> Config:
    return Config(
        bam=args.bam,
        gtf=args.gtf,
        out_prefix=args.out_prefix,
        fai=args.fai,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        keep_duplicates=args.keep_duplicates,
        count_secondary=args.count_secondary,
        min_unique_fraction=args.min_unique_fraction,
        min_depth=args.min_depth,
        strand_min_depth=args.strand_min_depth,
        max_gap=args.max_gap,
        min_region_length=args.min_region_length,
        min_covered_bases=args.min_covered_bases,
        min_covered_fraction=args.min_covered_fraction,
        annotation_mode=args.annotation_mode,
        nearest_feature_window=args.nearest_feature_window,
        library_strandedness=args.library_strandedness,
        strand_infer_max_reads=args.strand_infer_max_reads,
        strand_infer_min_confidence=args.strand_infer_min_confidence,
        gdna_min_dual_strand_fraction=args.gdna_min_dual_strand_fraction,
        gdna_max_strand_length_ratio_diff=args.gdna_max_strand_length_ratio_diff,
        gdna_min_profile_correlation=args.gdna_min_profile_correlation,
        gdna_min_covered_fraction=args.gdna_min_covered_fraction,
        gdna_min_depth_balance=args.gdna_min_depth_balance,
        gdna_flat_cv_threshold=args.gdna_flat_cv_threshold,
        novel_min_dominant_strand_fraction=args.novel_min_dominant_strand_fraction,
        threads=args.threads,
        emit_bed=args.emit_bed,
        emit_bedgraph=args.emit_bedgraph,
        verbose=args.verbose,
    )


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    cfg = args_to_config(args)
    logger = get_logger(cfg.verbose)

    # Import the pipeline lazily so `--help` works even without pysam installed.
    try:
        from .pipeline import run
    except ImportError as exc:  # pragma: no cover
        logger.error("Failed to import pipeline: %s", exc)
        return 2

    try:
        run(cfg)
    except (FileNotFoundError, ValueError, ImportError) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
