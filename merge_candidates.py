#!/usr/bin/env python3
"""Cross-sample consensus of candidate regions.

Takes the per-sample ``*.candidate_regions.tsv`` files produced by
``detect_gdna_vs_novel.py`` and keeps loci reproduced in at least --min-samples
samples, writing a consensus table + a consensus GTF of reproducible novel
transcripts. Pure polars/Python — no pysam, runs anywhere.

Example:
    python merge_candidates.py \
      --tsv A.candidate_regions.tsv B.candidate_regions.tsv \
            C.candidate_regions.tsv D.candidate_regions.tsv \
      --out-prefix cohort --min-samples 2
"""

from __future__ import annotations

import argparse
import sys

from gdna_rescue.crosssample import ConsensusConfig, run
from gdna_rescue.utils import get_logger


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="merge_candidates.py",
        description=(
            "Cross-sample consensus / reproducibility filter for candidate "
            "regions. Genuine novel transcripts should recur across samples; "
            "recurrent gDNA and multimapper loci are reported separately and are "
            "NOT added to the consensus GTF."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tsv", nargs="+", required=True,
                   help="Per-sample *.candidate_regions.tsv files (>=1).")
    p.add_argument("--sample-names", nargs="+", default=None,
                   help="Optional names matching --tsv order (default: from filename).")
    p.add_argument("--out-prefix", required=True, help="Output file prefix.")
    p.add_argument("--reference-gtf", default=None,
                   help="Reference GTF. If given, also writes "
                        "<prefix>.reference_plus_consensus.gtf — the "
                        "analysis-ready annotation (reference + reproducible novel "
                        "transcripts) to run featureCounts on the original STAR BAMs.")
    p.add_argument("--min-samples", type=int, default=2,
                   help="Keep loci present in at least this many samples.")
    p.add_argument("--min-reciprocal-overlap", type=float, default=0.85,
                   help="Two candidates match only if each covers >= this "
                        "fraction of the other.")
    p.add_argument("--ignore-strand", dest="strand_aware", action="store_false",
                   help="Match candidates regardless of dominant strand "
                        "(default: strand-aware).")
    p.add_argument("--include-bidirectional", action="store_true",
                   help="Also write reproducible bidirectional loci to the "
                        "consensus GTF (default: novel only).")

    ann = p.add_argument_group(
        "external annotation overlay",
        "Annotate the consensus coordinates with overlapping features from "
        "external GFF3/GFF/GTF files (e.g. the Ensembl Regulatory Build). "
        "Purely additive: adds columns to the consensus table and attributes to "
        "the GTF, and never changes a classification.",
    )
    ann.add_argument("--annotate", nargs="+", default=[], metavar="[LABEL=]FILE",
                     help="Annotation file(s), GFF3/GFF/GTF, optionally gzipped. "
                          "Prefix with 'label=' to set the output column prefix, "
                          "otherwise it is derived from the filename.")
    ann.add_argument("--annotate-labels", nargs="+", default=None,
                     help="Explicit column prefixes matching --annotate order.")
    ann.add_argument("--annotate-feature-types", nargs="+", default=None,
                     metavar="TYPE",
                     help="Keep only these column-3 feature types, e.g. "
                          "enhancer promoter CTCF_binding_site "
                          "open_chromatin_region (default: all types).")
    ann.add_argument("--annotate-nearest-window", type=int, default=10000,
                     help="When a region overlaps nothing, report the nearest "
                          "feature within this many bp (0 disables).")
    ann.add_argument("--no-annotate-names", dest="annotate_names",
                     action="store_false",
                     help="Do NOT append the overlapping feature type to "
                          "consensus transcript names. By default a locus on an "
                          "enhancer is named consensus_transcript_7-enhancer, so "
                          "the GTF and downstream count tables are "
                          "self-describing. Disable to keep bare IDs comparable "
                          "with an earlier un-annotated run.")
    ann.add_argument("--annotate-stranded", action="store_true",
                     help="Require the feature strand to match the region's "
                          "strand. Off by default: Ensembl regulatory features "
                          "are mostly strandless, and CTCF strand is motif "
                          "orientation, not transcription. Use for sources where "
                          "strand does mean transcription (e.g. CAGE peaks).")

    plots = p.add_argument_group(
        "per-transcript coverage plots",
        "Emit one IGV-like HTML per surviving transcript showing every sample's "
        "per-base coverage (split by strand and by unique/duplicate/multimapper) "
        "over the locus, with feature and reference-gene tracks and the full "
        "evidence panel. Needs the per-sample *.coverage.h5 written by "
        "detect_gdna_vs_novel.py --emit-coverage-store, plus h5py + plotly.",
    )
    plots.add_argument("--emit-plots", action="store_true",
                       help="Write plots into {out-prefix}_plots/.")
    plots.add_argument("--plot-shoulder", type=int, default=1000,
                       help="bp of context drawn either side of each locus.")
    plots.add_argument("--plot-all-passing", action="store_true",
                       help="Plot every locus passing --min-samples (incl. recurrent "
                            "gDNA / bidirectional / multimapper) as controls, not "
                            "just the reproducible novel transcripts.")
    plots.add_argument("--coverage-store", nargs="+", default=None,
                       metavar="H5",
                       help="Coverage-store paths matching --tsv order. Default: "
                            "auto-derived as {sample}.coverage.h5 next to each TSV.")

    p.add_argument("--verbose", action="store_true")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = ConsensusConfig(
        tsvs=args.tsv,
        sample_names=args.sample_names,
        out_prefix=args.out_prefix,
        min_samples=args.min_samples,
        min_reciprocal_overlap=args.min_reciprocal_overlap,
        strand_aware=args.strand_aware,
        include_bidirectional=args.include_bidirectional,
        reference_gtf=args.reference_gtf,
        annotate=args.annotate,
        annotate_labels=args.annotate_labels,
        annotate_feature_types=args.annotate_feature_types,
        annotate_nearest_window=args.annotate_nearest_window,
        annotate_stranded=args.annotate_stranded,
        annotate_names=args.annotate_names,
        emit_plots=args.emit_plots,
        plot_shoulder=args.plot_shoulder,
        plot_all_passing=args.plot_all_passing,
        coverage_stores=args.coverage_store,
        verbose=args.verbose,
    )
    logger = get_logger(cfg.verbose)
    try:
        run(cfg)
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
