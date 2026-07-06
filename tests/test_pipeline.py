"""End-to-end integration test. Requires pysam; auto-skipped where unavailable
(e.g. Windows). Run on Linux/macOS/WSL with pysam installed."""

import json
import os

import pytest

pysam = pytest.importorskip("pysam")

from gdna_rescue.config import Config
from gdna_rescue.pipeline import run
from gdna_rescue.classify import (
    LIKELY_GDNA,
    LIKELY_NOVEL,
    LIKELY_MULTIMAPPER,
    POSSIBLE_BIDIRECTIONAL,
)
from tests.generate_test_data import write_synthetic_bam_gtf


@pytest.fixture
def synthetic(tmp_path):
    bam, gtf = write_synthetic_bam_gtf(str(tmp_path / "data"))
    return bam, gtf, str(tmp_path / "out")


def test_full_pipeline_classifies_and_rescues(synthetic):
    bam, gtf, prefix = synthetic
    cfg = Config(
        bam=bam, gtf=gtf, out_prefix=prefix,
        library_strandedness="forward", threads=1, min_region_length=200,
    )
    summary = run(cfg)

    # Synthetic data contains: 1 symmetric gDNA region, 1 single-strand novel,
    # 1 asymmetric bidirectional region, 1 multimapper artifact.
    assert summary["n_likely_gDNA"] == 1
    assert summary["n_likely_novel_transcript"] == 1
    assert summary["n_possible_bidirectional_RNA"] == 1
    assert summary["n_likely_multimapper_artifact"] == 1
    # gDNA and the multimapper artifact are not rescued; the other two are.
    assert summary["n_rescued_unknown_transcripts"] == 2

    for suffix in (
        ".candidate_regions.tsv",
        ".unknown_transcripts.gtf",
        ".annotation_plus_unknowns.gtf",
        ".summary.json",
        ".candidate_regions.bed",
        ".gdna_mqc.tsv",
    ):
        assert os.path.exists(prefix + suffix), suffix

    # MultiQC file: bargraph of read assignment, keyed by sample.
    with open(prefix + ".gdna_mqc.tsv") as fh:
        mqc = fh.read()
    assert "# plot_type: 'bargraph'" in mqc
    mqc_lines = [l for l in mqc.splitlines() if not l.startswith("#")]
    mqc_header = mqc_lines[0].split("\t")
    for col in ("annotated", "novel_transcript", "gDNA"):
        assert col in mqc_header

    # Read assignment + gDNA QC are recorded in the summary JSON too.
    assert summary["gdna_contamination_qc"]["n_gDNA_regions"] == 1
    ra = summary["read_assignment_counts"]
    assert ra["annotated"] > 0            # gene-body reads
    assert ra[LIKELY_GDNA] > 0            # gDNA region carries reads
    # Categories sum to total unique reads.
    assert sum(ra.values()) == (
        ra["annotated"] + ra[LIKELY_NOVEL] + ra[POSSIBLE_BIDIRECTIONAL]
        + ra[LIKELY_GDNA] + ra[LIKELY_MULTIMAPPER] + ra["other_unannotated"]
    )

    # The rescued GTF must contain sequential unknown_transcript names.
    with open(prefix + ".unknown_transcripts.gtf") as fh:
        text = fh.read()
    assert 'transcript_id "unknown_transcript_1"' in text
    assert 'transcript_id "unknown_transcript_2"' in text


def test_auto_strandedness_detects_forward(synthetic):
    bam, gtf, prefix = synthetic
    cfg = Config(
        bam=bam, gtf=gtf, out_prefix=prefix + "_auto",
        library_strandedness="auto", threads=1, min_region_length=200,
    )
    run(cfg)
    with open(prefix + "_auto.summary.json") as fh:
        summary = json.load(fh)
    assert summary["inferred_library_strandedness"] == "forward"
    assert summary["strandedness_metrics"]["p_forward"] > 0.9
