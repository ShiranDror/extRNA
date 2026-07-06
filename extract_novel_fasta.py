#!/usr/bin/env python3
"""Extract FASTA sequences for novel transcripts from a genome + GTF.

Works on any GTF the tool emits — a per-sample `*.unknown_transcripts.gtf`, the
cross-sample `*.consensus_transcripts.gtf`, or a merged annotation (use
--source to pull only the rescued features). Sequences are concatenated over
exons and reverse-complemented for '-' strand transcripts (5'->3' orientation),
ready for BLAST / ORF / homology searches. Pure Python (+ optional pysam) — runs
anywhere.

Examples:
    python extract_novel_fasta.py --genome genome.fa \
        --gtf sample.unknown_transcripts.gtf --out sample.novel.fa

    python extract_novel_fasta.py --genome genome.fa \
        --gtf cohort.reference_plus_consensus.gtf \
        --source gdna_rescue_consensus --out cohort.novel.fa
"""

from __future__ import annotations

import argparse
import sys

from gdna_rescue.fasta import extract_transcript_fasta
from gdna_rescue.utils import get_logger


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="extract_novel_fasta.py",
        description="Extract novel-transcript sequences to FASTA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--genome", required=True, help="Genome FASTA (.fai built if absent).")
    p.add_argument("--gtf", required=True,
                   help="GTF with the transcripts to extract (e.g. "
                        "*.unknown_transcripts.gtf or *.consensus_transcripts.gtf).")
    p.add_argument("--out", required=True, help="Output FASTA path.")
    p.add_argument("--source", default=None,
                   help="Only extract transcripts from this GTF source column "
                        "(e.g. gdna_rescue or gdna_rescue_consensus). Useful when "
                        "the GTF also contains reference annotation.")
    p.add_argument("--line-width", type=int, default=60,
                   help="Wrap sequence lines at this width.")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    logger = get_logger(args.verbose)
    try:
        n = extract_transcript_fasta(
            args.gtf, args.genome, args.out,
            source_filter=args.source, line_width=args.line_width,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1
    logger.info("Wrote %d transcript sequences to %s", n, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
