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


def test_pair_geometry_filter_flags_half_mapped_region(tmp_path):
    """Region-level integration test for the pair-geometry filter (Config.pair_filter).

    Adds a region (chr_test:4000-4800) with 400 unique single-end reads and 1000
    paired plus-strand reads whose mate is unmapped. With pair_filter="concordant"
    (default) the paired reads are reclassified as noise regardless of MAPQ, so
    unique coverage (400) is swamped by multi coverage (1000) and the region is
    flagged likely_multimapper_artifact. With pair_filter="off" all 1400 reads
    count as unique on one dominant strand, so the same region is instead
    classified likely_novel_transcript.
    """
    bam, gtf = write_synthetic_bam_gtf(
        str(tmp_path / "data"), include_half_mapped=True
    )

    def _region_row(prefix):
        with open(prefix + ".candidate_regions.tsv") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row["chrom"] != "chr_test":
                    continue
                # TSV start/end are 1-based inclusive; overlap test against
                # the 0-based half-open region we generated (4000, 4800).
                if int(row["start"]) - 1 < 4800 and int(row["end"]) > 4000:
                    return row
        return None

    # (a) Default config: pair_filter="concordant".
    prefix_on = str(tmp_path / "out_concordant")
    cfg_on = Config(
        bam=bam, gtf=gtf, out_prefix=prefix_on,
        library_strandedness="forward", threads=1, min_region_length=200,
    )
    assert cfg_on.pair_filter == "concordant"
    summary_on = run(cfg_on)

    row_on = _region_row(prefix_on)
    assert row_on is not None, "no candidate region found over chr_test:4000-4800"
    assert row_on["class"] == LIKELY_MULTIMAPPER

    read_totals_on = summary_on["read_totals"]
    assert read_totals_on["n_half_mapped_reads"] == 1000
    assert read_totals_on["n_discordant_reads"] == 0

    # (b) Same data, pair_filter="off": paired mate-unmapped reads count as
    # unique again, so the region is instead a clean novel-transcript call.
    prefix_off = str(tmp_path / "out_off")
    cfg_off = Config(
        bam=bam, gtf=gtf, out_prefix=prefix_off,
        library_strandedness="forward", threads=1, min_region_length=200,
        pair_filter="off",
    )
    summary_off = run(cfg_off)

    row_off = _region_row(prefix_off)
    assert row_off is not None, "no candidate region found over chr_test:4000-4800"
    assert row_off["class"] == LIKELY_NOVEL

    read_totals_off = summary_off["read_totals"]
    assert read_totals_off["n_half_mapped_reads"] == 0
    assert read_totals_off["n_discordant_reads"] == 0


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
