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
