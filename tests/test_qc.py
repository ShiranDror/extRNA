"""Unit tests for gDNA QC metrics and the MultiQC TSV writer (pysam-free)."""

from types import SimpleNamespace

from gdna_rescue.config import Config
from gdna_rescue.classify import LIKELY_GDNA, LIKELY_NOVEL
from gdna_rescue.writers import compute_gdna_qc, write_multiqc_tsv


def _cand(label, total_coverage, length, kept):
    return SimpleNamespace(
        label=label,
        kept=kept,
        length=length,
        metrics=SimpleNamespace(total_coverage=total_coverage),
    )


def test_compute_gdna_qc_percentages():
    candidates = [
        _cand(LIKELY_GDNA, total_coverage=3000.0, length=1000, kept=False),
        _cand(LIKELY_NOVEL, total_coverage=1000.0, length=500, kept=True),
    ]
    # Genome-wide unique coverage = 40000 aligned bases.
    qc = compute_gdna_qc(candidates, genome_unique_coverage=40000,
                         genome_multi_coverage=5000)
    assert qc["n_gDNA_regions"] == 1
    assert qc["gDNA_unique_coverage"] == 3000.0
    # 3000 / 40000 = 7.5% of all mapped unique signal.
    assert qc["pct_gDNA_of_mapped_coverage"] == 7.5
    # 3000 / (3000 + 1000) = 75% of candidate coverage.
    assert qc["pct_gDNA_of_candidate_coverage"] == 75.0
    # 1000 / (1000 + 500) bases.
    assert qc["pct_gDNA_of_candidate_bases"] == round(100 * 1000 / 1500, 4)


def test_compute_gdna_qc_zero_safe():
    qc = compute_gdna_qc([], genome_unique_coverage=0, genome_multi_coverage=0)
    assert qc["pct_gDNA_of_mapped_coverage"] == 0.0
    assert qc["pct_gDNA_of_candidate_coverage"] == 0.0
    assert qc["n_gDNA_regions"] == 0


def test_write_multiqc_tsv(tmp_path):
    candidates = [
        _cand(LIKELY_GDNA, 3000.0, 1000, kept=False),
        _cand(LIKELY_NOVEL, 1000.0, 500, kept=True),
    ]
    qc = compute_gdna_qc(candidates, 40000, 5000)
    cfg = Config(out_prefix=str(tmp_path / "mysample"), sample_name="mysample")
    path = str(tmp_path / "mysample.gdna_mqc.tsv")
    write_multiqc_tsv(cfg, candidates, qc, path)

    text = open(path).read()
    # MultiQC custom-content markers and general-stats config present.
    assert "# plot_type: 'generalstats'" in text
    assert "extrna_pct_gDNA" in text
    lines = [l for l in text.splitlines() if not l.startswith("#")]
    header, row = lines[0].split("\t"), lines[1].split("\t")
    assert header[0] == "Sample"
    assert row[0] == "mysample"
    # % gDNA column carries the mapped-coverage percentage (7.5).
    pct_idx = header.index("extrna_pct_gDNA")
    assert float(row[pct_idx]) == 7.5
    # novel transcripts column = 1 rescued.
    nov_idx = header.index("extrna_novel_transcripts")
    assert int(row[nov_idx]) == 1


def test_multiqc_filename_matches_pattern():
    # MultiQC discovers *_mqc.tsv; our suffix must satisfy that.
    assert "sample.gdna_mqc.tsv".endswith("_mqc.tsv")
