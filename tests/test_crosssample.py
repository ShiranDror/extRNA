"""Unit tests for the cross-sample consensus (polars, no pysam)."""

import polars as pl
import pytest

from gdna_rescue.crosssample import (
    ConsensusConfig,
    build_consensus,
    load_candidates,
    cluster_rows,
    _reciprocal_overlap,
    _vote_class,
    CONS_NOVEL,
    CONS_GDNA,
    CONS_BIDIR,
    run,
)
from gdna_rescue.classify import (
    LIKELY_NOVEL,
    LIKELY_GDNA,
    POSSIBLE_BIDIRECTIONAL,
    LIKELY_MULTIMAPPER,
)


def _write_tsv(path, rows):
    """rows: list of (chrom, start, end, cls, strand)."""
    df = pl.DataFrame(
        {
            "region_id": [f"{r[0]}:{r[1]}-{r[2]}" for r in rows],
            "chrom": [r[0] for r in rows],
            "start": [r[1] for r in rows],
            "end": [r[2] for r in rows],
            "class": [r[3] for r in rows],
            "dominant_strand": [r[4] for r in rows],
            "unique_fraction": [1.0 for _ in rows],
            "dual_strand_fraction": [0.0 for _ in rows],
            "profile_correlation": [0.0 for _ in rows],
            "avg_depth": [20.0 for _ in rows],
        }
    )
    df.write_csv(path, separator="\t")
    return str(path)


def test_reciprocal_overlap():
    assert _reciprocal_overlap(0, 100, 0, 100) == 1.0
    assert _reciprocal_overlap(0, 100, 200, 300) == 0.0
    # 90 / 100 each way.
    assert _reciprocal_overlap(0, 100, 10, 100) == pytest.approx(0.9)
    # Asymmetric: small inside big -> min ratio is small.
    assert _reciprocal_overlap(0, 1000, 0, 100) == pytest.approx(0.1)


def test_vote_class_majority_and_tiebreak():
    assert _vote_class([LIKELY_NOVEL, LIKELY_NOVEL, LIKELY_GDNA]) == LIKELY_NOVEL
    # Tie -> most conservative (gDNA) wins.
    assert _vote_class([LIKELY_NOVEL, LIKELY_GDNA]) == LIKELY_GDNA
    assert _vote_class(
        [POSSIBLE_BIDIRECTIONAL, LIKELY_MULTIMAPPER]
    ) == LIKELY_MULTIMAPPER


def test_high_overlap_clusters_together():
    cfg = ConsensusConfig(min_reciprocal_overlap=0.85)
    rows = [
        {"chrom": "c", "_s": 100, "_e": 1100, "_group": ("c", "+")},
        {"chrom": "c", "_s": 120, "_e": 1120, "_group": ("c", "+")},  # ~0.96 overlap
    ]
    clusters = cluster_rows(rows, cfg)
    assert len(clusters) == 1


def test_low_overlap_does_not_cluster():
    cfg = ConsensusConfig(min_reciprocal_overlap=0.85)
    rows = [
        {"chrom": "c", "_s": 0, "_e": 1000, "_group": ("c", "+")},
        {"chrom": "c", "_s": 500, "_e": 1500, "_group": ("c", "+")},  # 0.5 overlap
    ]
    clusters = cluster_rows(rows, cfg)
    assert len(clusters) == 2


def test_strand_aware_separates_strands():
    cfg = ConsensusConfig(min_reciprocal_overlap=0.85, strand_aware=True)
    rows = [
        {"chrom": "c", "_s": 0, "_e": 1000, "_group": ("c", "+")},
        {"chrom": "c", "_s": 0, "_e": 1000, "_group": ("c", "-")},
    ]
    assert len(cluster_rows(rows, cfg)) == 2


def test_min_samples_filter_and_gtf(tmp_path):
    # Region X: novel in 3 samples (reproducible). Region Y: novel in 1 sample.
    a = _write_tsv(tmp_path / "A.candidate_regions.tsv",
                   [("chr1", 1000, 2000, LIKELY_NOVEL, "+"),
                    ("chr1", 5000, 6000, LIKELY_NOVEL, "+")])
    b = _write_tsv(tmp_path / "B.candidate_regions.tsv",
                   [("chr1", 1010, 2010, LIKELY_NOVEL, "+")])
    c = _write_tsv(tmp_path / "C.candidate_regions.tsv",
                   [("chr1", 990, 1990, LIKELY_NOVEL, "+")])

    cfg = ConsensusConfig(tsvs=[a, b, c], out_prefix=str(tmp_path / "cohort"),
                          min_samples=2)
    df = load_candidates(cfg)
    regions = build_consensus(df, cfg)

    by_pass = {r.consensus_id: r for r in regions}
    # The 1000-2000 locus is in 3 samples -> passes; the 5000-6000 in 1 -> fails.
    reproducible = [r for r in regions if r.passes_min_samples]
    assert len(reproducible) == 1
    assert reproducible[0].n_samples == 3
    assert reproducible[0].consensus_class == CONS_NOVEL
    assert reproducible[0].in_consensus_gtf
    assert reproducible[0].consensus_transcript_name == "consensus_transcript_1"

    singletons = [r for r in regions if not r.passes_min_samples]
    assert all(not r.in_consensus_gtf for r in singletons)


def test_recurrent_gdna_reported_but_not_in_gtf(tmp_path):
    # Same locus called gDNA in 3 samples -> recurrent_gDNA, not in GTF.
    a = _write_tsv(tmp_path / "A.candidate_regions.tsv",
                   [("chr2", 1000, 3000, LIKELY_GDNA, "+")])
    b = _write_tsv(tmp_path / "B.candidate_regions.tsv",
                   [("chr2", 1005, 3005, LIKELY_GDNA, "+")])
    c = _write_tsv(tmp_path / "C.candidate_regions.tsv",
                   [("chr2", 995, 2995, LIKELY_GDNA, "+")])
    cfg = ConsensusConfig(tsvs=[a, b, c], out_prefix=str(tmp_path / "co"),
                          min_samples=2)
    summary = run(cfg)
    assert summary["n_recurrent_gDNA"] == 1
    assert summary["n_reproducible_novel"] == 0
    assert summary["n_written_to_consensus_gtf"] == 0
    # GTF exists but contains no feature lines.
    with open(str(tmp_path / "co.consensus_transcripts.gtf")) as fh:
        body = [l for l in fh if not l.startswith("#")]
    assert body == []


def test_reference_plus_consensus_gtf(tmp_path):
    # A tiny reference GTF with one gene.
    ref = tmp_path / "ref.gtf"
    ref.write_text(
        'chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "geneA"; gene_name "geneA";\n'
        'chr1\tsrc\texon\t100\t500\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";\n'
    )
    a = _write_tsv(tmp_path / "A.candidate_regions.tsv",
                   [("chr1", 1000, 2000, LIKELY_NOVEL, "+")])
    b = _write_tsv(tmp_path / "B.candidate_regions.tsv",
                   [("chr1", 1010, 2010, LIKELY_NOVEL, "+")])
    cfg = ConsensusConfig(tsvs=[a, b], out_prefix=str(tmp_path / "cohort"),
                          min_samples=2, reference_gtf=str(ref))
    run(cfg)
    merged = str(tmp_path / "cohort.reference_plus_consensus.gtf")
    import os
    assert os.path.exists(merged)
    text = open(merged).read()
    # Reference gene preserved AND consensus transcript appended.
    assert 'gene_id "geneA"' in text
    assert 'transcript_id "consensus_transcript_1"' in text
    # Consensus feature carries genomic union-span coords (1000..2010).
    assert "\ttranscript\t1000\t2010\t" in text
    assert "\texon\t1000\t2010\t" in text


def test_chrom_column_forced_to_string(tmp_path):
    # Regression: many numeric chromosomes followed by "MT" must not make polars
    # infer an integer chrom column and then fail to parse "MT".
    rows = [(str(i), 1000 + i, 2000 + i, LIKELY_NOVEL, "+") for i in range(1, 23)]
    rows.append(("MT", 500, 1500, LIKELY_NOVEL, "+"))
    path = _write_tsv(tmp_path / "A.candidate_regions.tsv", rows)
    cfg = ConsensusConfig(tsvs=[path], out_prefix=str(tmp_path / "co"), min_samples=1)
    df = load_candidates(cfg)
    assert str(df.schema["chrom"]) in ("String", "Utf8")
    # "MT" is present and grouped correctly.
    regions = build_consensus(df, cfg)
    assert any(r.chrom == "MT" for r in regions)


def test_numeric_only_chrom_still_string(tmp_path):
    # Even when every chromosome value is numeric, chrom must load as string so
    # a later sample with "X"/"MT" merges consistently.
    path = _write_tsv(tmp_path / "A.candidate_regions.tsv",
                      [("1", 1000, 2000, LIKELY_NOVEL, "+")])
    cfg = ConsensusConfig(tsvs=[path], out_prefix=str(tmp_path / "co"), min_samples=1)
    df = load_candidates(cfg)
    assert str(df.schema["chrom"]) in ("String", "Utf8")


def test_end_to_end_files_written(tmp_path):
    a = _write_tsv(tmp_path / "A.candidate_regions.tsv",
                   [("chr1", 1000, 2000, LIKELY_NOVEL, "+")])
    b = _write_tsv(tmp_path / "B.candidate_regions.tsv",
                   [("chr1", 1010, 2010, LIKELY_NOVEL, "+")])
    cfg = ConsensusConfig(tsvs=[a, b], out_prefix=str(tmp_path / "cohort"),
                          min_samples=2)
    summary = run(cfg)
    assert summary["n_reproducible_novel"] == 1
    import os
    for suffix in (".consensus_regions.tsv", ".consensus_transcripts.gtf",
                   ".consensus_summary.json"):
        assert os.path.exists(str(tmp_path / "cohort") + suffix)
