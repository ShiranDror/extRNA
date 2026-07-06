"""Unit tests for gDNA QC metrics, read assignment, and the MultiQC TSV writer.

Pysam-free: candidates are mocked with SimpleNamespace.
"""

from types import SimpleNamespace

from gdna_rescue.config import Config
from gdna_rescue.classify import (
    LIKELY_GDNA,
    LIKELY_NOVEL,
    LIKELY_MULTIMAPPER,
    POSSIBLE_BIDIRECTIONAL,
)
from gdna_rescue.writers import (
    compute_gdna_qc,
    compute_read_assignment,
    write_multiqc_tsv,
)


def _cand(label, total_coverage, length, kept):
    return SimpleNamespace(
        label=label,
        kept=kept,
        length=length,
        metrics=SimpleNamespace(total_coverage=total_coverage),
    )


# --- coverage-based gDNA QC (summary.json) --------------------------------

def test_compute_gdna_qc_percentages():
    candidates = [
        _cand(LIKELY_GDNA, total_coverage=3000.0, length=1000, kept=False),
        _cand(LIKELY_NOVEL, total_coverage=1000.0, length=500, kept=True),
    ]
    qc = compute_gdna_qc(candidates, genome_unique_coverage=40000,
                         genome_multi_coverage=5000)
    assert qc["n_gDNA_regions"] == 1
    assert qc["pct_gDNA_of_mapped_coverage"] == 7.5
    assert qc["pct_gDNA_of_candidate_coverage"] == 75.0


def test_compute_gdna_qc_zero_safe():
    qc = compute_gdna_qc([], genome_unique_coverage=0, genome_multi_coverage=0)
    assert qc["pct_gDNA_of_mapped_coverage"] == 0.0
    assert qc["n_gDNA_regions"] == 0


# --- read-count assignment (MultiQC bargraph) ----------------------------

def test_compute_read_assignment_partitions_reads():
    candidates = [
        _cand(LIKELY_GDNA, 0, 1000, False),
        _cand(LIKELY_NOVEL, 0, 500, True),
        _cand(POSSIBLE_BIDIRECTIONAL, 0, 400, True),
        _cand(LIKELY_MULTIMAPPER, 0, 300, False),
    ]
    region_reads = [800, 120, 60, 40]  # per-candidate unique read counts
    cats = compute_read_assignment(
        candidates, region_reads,
        total_unique_reads=1_000_000, total_annotated_reads=950_000,
    )
    assert cats["annotated"] == 950_000
    assert cats[LIKELY_GDNA] == 800
    assert cats[LIKELY_NOVEL] == 120
    assert cats[POSSIBLE_BIDIRECTIONAL] == 60
    assert cats[LIKELY_MULTIMAPPER] == 40
    # remainder = 1,000,000 - 950,000 - (800+120+60+40)
    assert cats["other_unannotated"] == 1_000_000 - 950_000 - 1020
    # Categories should sum back to the total unique reads.
    assert sum(v for v in cats.values()) == 1_000_000


def test_read_assignment_clamps_negative_remainder():
    # If annotated + candidate reads exceed the total (edge/midpoint effects),
    # other_unannotated must not go negative.
    candidates = [_cand(LIKELY_NOVEL, 0, 500, True)]
    cats = compute_read_assignment(candidates, [100],
                                   total_unique_reads=100, total_annotated_reads=100)
    assert cats["other_unannotated"] == 0


def test_write_multiqc_tsv_bargraph(tmp_path):
    candidates = [
        _cand(LIKELY_GDNA, 0, 1000, False),
        _cand(LIKELY_NOVEL, 0, 500, True),
    ]
    cats = compute_read_assignment(candidates, [800, 120],
                                   total_unique_reads=1000, total_annotated_reads=50)
    cfg = Config(out_prefix=str(tmp_path / "mysample"), sample_name="mysample")
    path = str(tmp_path / "mysample.gdna_mqc.tsv")
    write_multiqc_tsv(cfg, cats, path)

    text = open(path).read()
    assert "# plot_type: 'bargraph'" in text
    lines = [l for l in text.splitlines() if not l.startswith("#")]
    header, row = lines[0].split("\t"), lines[1].split("\t")
    assert header[0] == "Sample"
    assert row[0] == "mysample"
    # gDNA column = 800 reads.
    assert int(row[header.index("gDNA")]) == 800
    assert int(row[header.index("novel_transcript")]) == 120
    assert int(row[header.index("annotated")]) == 50


def test_multiqc_filename_matches_pattern():
    assert "sample.gdna_mqc.tsv".endswith("_mqc.tsv")
