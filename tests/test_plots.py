"""Smoke tests for the per-transcript coverage plots (plotly + h5py, no pysam)."""

from __future__ import annotations

import os

import numpy as np
import polars as pl
import pytest

pytest.importorskip("h5py")
pytest.importorskip("plotly")

from gdna_rescue import plots  # noqa: E402
from gdna_rescue.coverage_store import sparsify, write_coverage_store  # noqa: E402
from gdna_rescue.crosssample import ConsensusConfig  # noqa: E402
from gdna_rescue.crosssample import ConsensusRegion  # noqa: E402
from gdna_rescue.overlay import parse_overlay  # noqa: E402


def _make_store(path, chrom, length, plus_span, sample):
    up = np.zeros(length, dtype=np.int32)
    s, e = plus_span
    up[s:e] = 25
    dp = np.zeros(length, dtype=np.int32)
    dp[s:e] = 6  # some duplicate inflation to render a second channel
    pos, sparse = sparsify({"unique_plus": up, "dup_plus": dp})
    write_coverage_store(path, {chrom: (pos, sparse)}, sample_name=sample,
                         chrom_sizes={chrom: length}, strandedness="reverse")


def _region():
    return ConsensusRegion(
        consensus_id="1:1000-1400:+",
        chrom="1", start=1000, end=1400, strand="+",
        n_samples=2, n_members=2,
        samples=["A", "B"],
        member_region_ids=["1:1000-1400", "1:1010-1395"],
        member_classes=["likely_novel_transcript", "likely_novel_transcript"],
        majority_class="likely_novel_transcript",
        consensus_class="reproducible_novel",
        class_agreement=1.0,
        mean_unique_fraction=0.98,
        mean_dual_strand_fraction=0.05,
        mean_profile_correlation=0.1,
        mean_avg_depth=25.0,
        passes_min_samples=True,
        in_consensus_gtf=True,
        consensus_transcript_name="consensus_transcript_1",
        member_samples=["A", "B"],
    )


def _df():
    return pl.DataFrame({
        "sample": ["A", "B"],
        "region_id": ["1:1000-1400", "1:1010-1395"],
        "class": ["likely_novel_transcript", "likely_novel_transcript"],
        "avg_depth": [25.0, 24.0],
        "max_depth": [31, 30],
        "covered_fraction": [0.99, 0.98],
        "unique_fraction": [0.99, 0.97],
        "dual_strand_fraction": [0.04, 0.06],
        "profile_correlation": [0.1, 0.12],
        "dominant_strand": ["+", "+"],
        "context_label": ["intergenic", "intergenic"],
        "nearest_feature_id": ["ENSG1", "ENSG1"],
        "nearest_feature_distance": [500, 480],
        "reason_for_classification": ["single strand dominant", "single strand dominant"],
    })


def test_build_plots_smoke(tmp_path):
    length = 3000
    # C is a non-contributing sample with no coverage at the locus -> still drawn.
    store_by_sample = {}
    for sample, span in (("A", (1000, 1400)), ("B", (1010, 1395)), ("C", None)):
        p = str(tmp_path / f"{sample}.coverage.h5")
        if span is None:
            _make_store(p, "1", length, (5, 6), sample)  # coverage far from locus
        else:
            _make_store(p, "1", length, span, sample)
        store_by_sample[sample] = p

    gff = tmp_path / "reg.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "1\tt\tenhancer\t1100\t1300\t.\t.\t.\tID=E1\n",
        encoding="utf-8",
    )
    idx = parse_overlay(str(gff), label="reg")

    cfg = ConsensusConfig(tsvs=["A.candidate_regions.tsv"], out_prefix="cohort",
                          plot_shoulder=200, verbose=False)
    out_dir = str(tmp_path / "cohort_plots")

    n = plots.build_plots(
        [_region()], _df(), [idx], ["reg"], cfg,
        store_by_sample=store_by_sample, sample_order=["A", "B", "C"],
        out_dir=out_dir,
    )
    assert n == 1

    html_path = os.path.join(out_dir, "consensus_transcript_1.html")
    assert os.path.exists(html_path)
    assert os.path.exists(os.path.join(out_dir, "index.html"))
    # plotly.min.js is written alongside (include_plotlyjs="directory").
    assert os.path.exists(os.path.join(out_dir, "plotly.min.js"))

    page = open(html_path, encoding="utf-8").read()
    assert "consensus_transcript_1" in page
    for token in ("unique", "dup", "multi", "Evidence", "enhancer"):
        assert token in page
    # Every sample, including the non-contributing C, has a panel.
    for sample in ("A", "B", "C"):
        assert sample in page
