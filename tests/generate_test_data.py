"""Synthetic test-data generator.

Two levels:

1. ``make_archetype_coverage`` builds per-strand coverage arrays for the three
   biological archetypes the classifier must distinguish. These need only numpy
   and are used by the unit tests (runnable on any platform, no pysam).

2. ``write_synthetic_bam_gtf`` writes a tiny BAM + GTF exercising the full
   pipeline. It requires pysam and is therefore skipped automatically where
   pysam is unavailable (e.g. Windows).

Run directly to print the classifier's decision on each archetype:

    python -m tests.generate_test_data
"""

from __future__ import annotations

from typing import Dict, Tuple

import numpy as np


def make_archetype_coverage(
    seed: int = 7, length: int = 1200, depth: int = 12
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Return {name: (plus_cov, minus_cov)} for three archetypes.

    * gdna_symmetric  : both strands broad, similar length, correlated profiles.
    * novel_single    : one dominant strand, other strand ~empty.
    * bidir_asymmetric: both strands present but offset/partial, uncorrelated.
    """
    rng = np.random.default_rng(seed)

    # A smooth shared shape so gDNA strands are genuinely correlated.
    x = np.linspace(0, 6 * np.pi, length)
    shape = depth * (1.2 + 0.6 * np.sin(x))  # always positive, undulating

    # 1) gDNA: both strands follow the same shape (dsDNA fragments both ways).
    gdna_plus = np.clip(shape + rng.normal(0, 0.8, length), 0, None).astype(np.int32)
    gdna_minus = np.clip(shape + rng.normal(0, 0.8, length), 0, None).astype(np.int32)

    # 2) Novel single-strand transcript: plus carries essentially all signal.
    novel_plus = np.clip(shape + rng.normal(0, 0.8, length), 0, None).astype(np.int32)
    novel_minus = (rng.random(length) < 0.02).astype(np.int32)  # rare noise

    # 3) Bidirectional but asymmetric: plus on the left, minus on the right,
    #    only partial overlap and different shapes -> low correlation.
    bidir_plus = np.zeros(length, dtype=np.int32)
    bidir_minus = np.zeros(length, dtype=np.int32)
    left = slice(0, int(length * 0.55))
    right = slice(int(length * 0.45), length)
    bidir_plus[left] = np.clip(
        depth * 1.0 + rng.normal(0, 1.0, left.stop - left.start), 0, None
    ).astype(np.int32)
    bidir_minus[right] = np.clip(
        depth * 0.7 + rng.normal(0, 1.0, right.stop - right.start), 0, None
    ).astype(np.int32)

    return {
        "gdna_symmetric": (gdna_plus, gdna_minus),
        "novel_single": (novel_plus, novel_minus),
        "bidir_asymmetric": (bidir_plus, bidir_minus),
    }


def write_synthetic_bam_gtf(out_dir: str) -> Tuple[str, str]:
    """Write a small BAM + GTF for full-pipeline testing (requires pysam).

    Returns (bam_path, gtf_path). Raises ImportError if pysam is unavailable.
    """
    import os
    import pysam  # noqa: F401  (import error surfaces to caller)

    os.makedirs(out_dir, exist_ok=True)
    chrom = "chr_test"
    chrom_len = 20000
    gtf_path = os.path.join(out_dir, "synthetic.gtf")
    bam_path = os.path.join(out_dir, "synthetic.bam")

    # One annotated gene so discovery has something to mask and strandedness
    # inference has exons to sample.
    with open(gtf_path, "w") as fh:
        g_s, g_e = 1000, 3000  # 1-based inclusive
        attrs = 'gene_id "geneA"; transcript_id "txA"; gene_name "geneA";'
        fh.write(f"{chrom}\tsim\tgene\t{g_s}\t{g_e}\t.\t+\t.\t{attrs}\n")
        fh.write(f"{chrom}\tsim\ttranscript\t{g_s}\t{g_e}\t.\t+\t.\t{attrs}\n")
        fh.write(f"{chrom}\tsim\texon\t{g_s}\t{g_e}\t.\t+\t.\t{attrs}\n")

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": chrom, "LN": chrom_len}]}

    rng = np.random.default_rng(3)
    read_len = 75

    def add_reads(recs, start, end, strand, n, mapq=255, secondary=False):
        for i in range(n):
            pos = int(rng.integers(start, max(start + 1, end - read_len)))
            a = pysam.AlignedSegment()
            a.query_name = f"r{strand}_{start}_{i}_{mapq}"
            a.query_sequence = "A" * read_len
            flag = 16 if strand == "-" else 0
            if secondary:
                flag |= 256
            a.flag = flag
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = mapq  # 255 = STAR unique; low = multimapper
            a.cigar = [(0, read_len)]
            a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            recs.append(a)

    recs = []
    # Annotated gene body: forward reads (sense), should be masked out.
    add_reads(recs, 1000, 3000, "+", 400)
    # gDNA region: both strands, symmetric, broad, 6000-8000.
    add_reads(recs, 6000, 8000, "+", 500)
    add_reads(recs, 6000, 8000, "-", 500)
    # Novel single-strand region: plus only, 11000-12500.
    add_reads(recs, 11000, 12500, "+", 400)
    # Bidirectional asymmetric: plus left, minus right, partial overlap.
    add_reads(recs, 15000, 16200, "+", 300)
    add_reads(recs, 15800, 17000, "-", 220)
    # Multimapper artifact: dense-enough unique coverage to be discovered as one
    # region, but swamped by multimapped reads (MAPQ 1) -> flagged, not rescued.
    add_reads(recs, 18000, 18800, "+", 400, mapq=255)
    add_reads(recs, 18000, 18800, "+", 1000, mapq=1)

    recs.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for r in recs:
            out.write(r)
    pysam.index(bam_path)
    return bam_path, gtf_path


if __name__ == "__main__":
    from gdna_rescue.config import Config
    from gdna_rescue.classify import compute_region_metrics, classify_region

    cfg = Config()
    for name, (p, m) in make_archetype_coverage().items():
        metrics = compute_region_metrics(p, m, cfg)
        label, reason, _ = classify_region(metrics, cfg)
        print(f"{name:18s} -> {label}")
        print(f"    {reason}")
