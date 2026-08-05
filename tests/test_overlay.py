"""Tests for the external annotation overlay (pure stdlib, no pysam)."""

from __future__ import annotations

import gzip
import json
import os

import pytest

from gdna_rescue.crosssample import (
    CONS_GDNA,
    CONS_NOVEL,
    ConsensusConfig,
    annotate_consensus,
    annotation_crosstab,
    apply_annotation_names,
    build_consensus,
    consensus_to_dataframe,
    write_consensus_gtf,
    write_consensus_summary,
)
from gdna_rescue.overlay import (
    annotation_name_token,
    derive_label,
    load_sources,
    normalise_chrom,
    overlay_columns,
    parse_attributes,
    parse_overlay,
    sanitise_name_token,
)

# Real Ensembl Regulatory Build lines (release 116 formatting), 1-based inclusive.
ENSEMBL_GFF3 = """\
##gff-version 3
# Gene IDs in promoter attributes are from Ensembl release 114
1\tEnsembl\tpromoter\t10936\t11436\t.\t.\t.\tID=ENSR1_958;extended_start=9953;extended_end=11436;gene_id=ENSG00000290825;gene_name=DDX11L16;gene_biotype=lncRNA;color=#ff0000
1\tEnsembl\tCTCF_binding_site\t11222\t11243\t.\t+\t.\tID=ENSR1_53X2;color=#40e0d0
1\tEnsembl\tenhancer\t11437\t11649\t.\t.\t.\tID=ENSR1_88N;color=#faca00
1\tEnsembl\tenhancer\t12686\t14671\t.\t.\t.\tID=ENSR1_CZ;color=#faca00
1\tEnsembl\topen_chromatin_region\t50000\t50500\t.\t.\t.\tID=ENSR1_OCR1;color=#000000
2\tEnsembl\tenhancer\t900\t1100\t.\t-\t.\tID=ENSR2_E1;color=#faca00
"""


def _write(tmp_path, name, text, gz=False):
    path = os.path.join(str(tmp_path), name)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return path


# --------------------------------------------------------------------------- #
# Attribute parsing / labels / chrom naming
# --------------------------------------------------------------------------- #

def test_parse_gff3_attributes():
    attrs = parse_attributes("ID=ENSR1_958;gene_name=DDX11L16;color=#ff0000")
    assert attrs["ID"] == "ENSR1_958"
    assert attrs["gene_name"] == "DDX11L16"
    assert attrs["color"] == "#ff0000"


def test_parse_gtf_attributes_still_work():
    """The overlay must accept GTF too, not just GFF3."""
    attrs = parse_attributes('gene_id "ENSG001"; transcript_id "ENST001";')
    assert attrs == {"gene_id": "ENSG001", "transcript_id": "ENST001"}


def test_gff3_values_are_percent_decoded():
    assert parse_attributes("Name=foo%2Cbar%20baz")["Name"] == "foo,bar baz"


@pytest.mark.parametrize(
    "raw,expected",
    [("chr1", "1"), ("1", "1"), ("chrX", "X"), ("X", "X"),
     ("chrM", "MT"), ("M", "MT"), ("MT", "MT"), ("chrMT", "MT")],
)
def test_normalise_chrom(raw, expected):
    assert normalise_chrom(raw) == expected


def test_derive_label_from_ensembl_filename():
    assert derive_label("Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz") == \
        "regulatory"
    assert derive_label("Homo_sapiens.GRCh38.EMARs.v116.gff.gz") == "emar"


def test_derive_label_uniquifies():
    used = ["regulatory"]
    assert derive_label("x.regulatory_features.gff3", used) == "regulatory2"


# --------------------------------------------------------------------------- #
# Parsing and overlap
# --------------------------------------------------------------------------- #

def test_parse_overlay_counts_and_types(tmp_path):
    idx = parse_overlay(_write(tmp_path, "reg.gff3", ENSEMBL_GFF3), label="reg")
    assert idx.n_features == 6
    assert idx.type_counts == {
        "promoter": 1, "CTCF_binding_site": 1,
        "enhancer": 3, "open_chromatin_region": 1,
    }
    # Ensembl names chromosomes '1'/'2' with no 'chr' prefix.
    assert idx.chroms == {"1", "2"}


def test_gzipped_input_is_supported(tmp_path):
    idx = parse_overlay(_write(tmp_path, "reg.gff3.gz", ENSEMBL_GFF3, gz=True))
    assert idx.n_features == 6


def test_coordinates_are_converted_to_half_open(tmp_path):
    """GFF3 1-based inclusive 11437..11649 -> 0-based half-open [11436, 11649)."""
    idx = parse_overlay(_write(tmp_path, "reg.gff3", ENSEMBL_GFF3))
    enh = [f for f in idx.features["1"] if f.fid == "ENSR1_88N"][0]
    assert (enh.start, enh.end) == (11436, 11649)

    # A region ending exactly at the feature start must NOT overlap.
    assert idx.overlap("1", 11000, 11436) == [] or all(
        f.fid != "ENSR1_88N" for f in idx.overlap("1", 11000, 11436)
    )
    # One base into it must overlap.
    assert any(f.fid == "ENSR1_88N" for f in idx.overlap("1", 11000, 11437))


def test_overlap_matches_across_chr_prefix(tmp_path):
    """A UCSC-style 'chr1' region must still hit an Ensembl-style '1' feature."""
    idx = parse_overlay(_write(tmp_path, "reg.gff3", ENSEMBL_GFF3))
    hits = idx.overlap("chr1", 12700, 12800)
    assert [f.fid for f in hits] == ["ENSR1_CZ"]


def test_long_feature_is_found_from_a_late_query(tmp_path):
    """The max-length lookback must not miss a long feature starting far left."""
    gff = (
        "1\tEnsembl\tenhancer\t100\t100000\t.\t.\t.\tID=LONG\n"
        "1\tEnsembl\tenhancer\t99000\t99100\t.\t.\t.\tID=SHORT\n"
    )
    idx = parse_overlay(_write(tmp_path, "long.gff3", gff))
    hits = {f.fid for f in idx.overlap("1", 99050, 99060)}
    assert hits == {"LONG", "SHORT"}


def test_feature_type_filter(tmp_path):
    path = _write(tmp_path, "reg.gff3", ENSEMBL_GFF3)
    idx = parse_overlay(path, feature_types=["enhancer"])
    assert idx.n_features == 3
    assert set(idx.type_counts) == {"enhancer"}


def test_stranded_matching_is_opt_in(tmp_path):
    path = _write(tmp_path, "reg.gff3", ENSEMBL_GFF3)
    # chr2 enhancer is on '-'. Strand-agnostic (default) hits from either strand.
    loose = parse_overlay(path)
    assert loose.overlap("2", 950, 1000, "+")
    # Stranded matching drops the mismatch.
    strict = parse_overlay(path, stranded=True)
    assert not strict.overlap("2", 950, 1000, "+")
    assert strict.overlap("2", 950, 1000, "-")


def test_footprint_is_union_not_sum(tmp_path):
    """Overlapping same-type features must not double-count covered bp."""
    gff = (
        "1\tEnsembl\tenhancer\t1\t100\t.\t.\t.\tID=A\n"
        "1\tEnsembl\tenhancer\t51\t150\t.\t.\t.\tID=B\n"
    )
    idx = parse_overlay(_write(tmp_path, "u.gff3", gff))
    assert idx.covered_bp_by_type == {"enhancer": 150}


def test_load_sources_label_syntax_and_duplicates(tmp_path):
    path = _write(tmp_path, "reg.gff3", ENSEMBL_GFF3)
    (idx,) = load_sources([f"myreg={path}"])
    assert idx.label == "myreg"

    with pytest.raises(ValueError, match="duplicate annotation label"):
        load_sources([f"a={path}", f"a={path}"])

    with pytest.raises(FileNotFoundError):
        load_sources([os.path.join(str(tmp_path), "nope.gff3")])


def test_windows_path_is_not_mistaken_for_a_label(tmp_path):
    """A drive-letter path has no '=', but guard the label split anyway."""
    path = _write(tmp_path, "reg.gff3", ENSEMBL_GFF3)
    (idx,) = load_sources([path])
    assert idx.n_features == 6


# --------------------------------------------------------------------------- #
# Integration with the consensus step
# --------------------------------------------------------------------------- #

def _consensus_rows():
    """Two loci reproduced in 2 samples: one on an enhancer, one on nothing."""
    import polars as pl

    rows = []
    for sample in ("A", "B"):
        rows.append({
            "region_id": f"1:12700-13000", "chrom": "1", "start": 12701,
            "end": 13000, "class": "likely_novel_transcript",
            "dominant_strand": "+", "unique_fraction": 0.99,
            "dual_strand_fraction": 0.05, "profile_correlation": 0.1,
            "avg_depth": 30.0, "sample": sample,
        })
        rows.append({
            "region_id": f"1:800000-800400", "chrom": "1", "start": 800001,
            "end": 800400, "class": "likely_gDNA",
            "dominant_strand": "+", "unique_fraction": 0.98,
            "dual_strand_fraction": 0.9, "profile_correlation": 0.9,
            "avg_depth": 12.0, "sample": sample,
        })
    return pl.DataFrame(rows)


def test_annotation_adds_columns_and_hits(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    assert len(regions) == 2

    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    stats = annotate_consensus(regions, indexes, cfg)

    assert stats[0]["n_features_loaded"] == 6
    assert stats[0]["n_regions_with_overlap"] == 1

    on_enhancer = [c for c in regions if c.start == 12700][0]
    assert on_enhancer.annotations["reg_n"] == 1
    assert on_enhancer.annotations["reg_types"] == "enhancer:1"
    assert on_enhancer.annotations["reg_ids"] == "ENSR1_CZ"
    assert on_enhancer.annotations["reg_overlap_bp"] == 300
    assert on_enhancer.annotations["reg_overlap_frac"] == 1.0

    empty = [c for c in regions if c.start == 800000][0]
    assert empty.annotations["reg_n"] == 0
    assert empty.annotations["reg_types"] == "NA"
    # Nothing within the default window out at 800 kb.
    assert empty.annotations["reg_nearest_distance"] == -1


def test_promoter_gene_names_are_reported(tmp_path):
    """Ensembl promoters carry gene_id/gene_name; surface them."""
    import polars as pl

    rows = [{
        "region_id": "1:11000-11200", "chrom": "1", "start": 11001, "end": 11200,
        "class": "likely_novel_transcript", "dominant_strand": "-",
        "unique_fraction": 0.99, "dual_strand_fraction": 0.1,
        "profile_correlation": 0.0, "avg_depth": 20.0, "sample": s,
    } for s in ("A", "B")]

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    assert regions[0].annotations["reg_genes"] == "DDX11L16"
    assert "promoter" in str(regions[0].annotations["reg_types"])


def test_nearest_feature_is_reported_within_window(tmp_path):
    import polars as pl

    # 2 kb downstream of the open_chromatin_region ending at 50500.
    rows = [{
        "region_id": "1:52000-52300", "chrom": "1", "start": 52001, "end": 52300,
        "class": "likely_novel_transcript", "dominant_strand": "+",
        "unique_fraction": 0.99, "dual_strand_fraction": 0.1,
        "profile_correlation": 0.0, "avg_depth": 20.0, "sample": s,
    } for s in ("A", "B")]

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2,
                          annotate_nearest_window=10000)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    assert regions[0].annotations["reg_n"] == 0
    assert regions[0].annotations["reg_nearest"] == "open_chromatin_region:ENSR1_OCR1"
    assert regions[0].annotations["reg_nearest_distance"] == 52000 - 50500


def test_chrom_naming_mismatch_raises_instead_of_silent_na(tmp_path):
    """A silent table of NA is the worst outcome; fail loudly instead."""
    import polars as pl

    rows = [{
        "region_id": "scaffold_9:100-400", "chrom": "scaffold_9", "start": 101,
        "end": 400, "class": "likely_novel_transcript", "dominant_strand": "+",
        "unique_fraction": 0.99, "dual_strand_fraction": 0.1,
        "profile_correlation": 0.0, "avg_depth": 20.0, "sample": s,
    } for s in ("A", "B")]

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])

    with pytest.raises(ValueError, match="shares no chromosome names"):
        annotate_consensus(regions, indexes, cfg)


def test_table_keeps_overlay_columns_when_nothing_overlaps(tmp_path):
    """Zero hits must still yield the columns, not a silently narrower table."""
    import polars as pl

    rows = [{
        "region_id": "1:700000-700300", "chrom": "1", "start": 700001,
        "end": 700300, "class": "likely_novel_transcript",
        "dominant_strand": "+", "unique_fraction": 0.99,
        "dual_strand_fraction": 0.1, "profile_correlation": 0.0,
        "avg_depth": 20.0, "sample": s,
    } for s in ("A", "B")]

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    df = consensus_to_dataframe(regions, overlay_labels=["reg"])
    for col in overlay_columns("reg"):
        assert col in df.columns


def test_classification_is_untouched_by_annotation(tmp_path):
    """The overlay is additive: classes before == classes after."""
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    before = [(c.consensus_class, c.in_consensus_gtf) for c in regions]

    indexes = load_sources([_write(tmp_path, "reg.gff3", ENSEMBL_GFF3)])
    annotate_consensus(regions, indexes, cfg)

    assert [(c.consensus_class, c.in_consensus_gtf) for c in regions] == before


def test_crosstab_separates_novel_from_gdna(tmp_path):
    """The class x feature-type table is the built-in control."""
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    tab = annotation_crosstab(regions, ["reg"])["reg"]
    assert tab["n_regions_with_overlap"] == 1
    assert tab["by_consensus_class"][CONS_NOVEL]["enhancer"] == 1
    assert "enhancer" not in tab["by_consensus_class"][CONS_GDNA]


def test_overlay_types_reach_the_consensus_gtf(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    out = os.path.join(str(tmp_path), "c.consensus_transcripts.gtf")
    assert write_consensus_gtf(regions, out) == 1
    with open(out) as fh:
        text = fh.read()
    assert 'reg_types "enhancer:1"' in text


# --------------------------------------------------------------------------- #
# Annotation-derived transcript names
# --------------------------------------------------------------------------- #

def test_name_token_is_first_type_in_the_types_column():
    """The suffix must be readable straight off <label>_types (sorted order)."""
    types = {"reg": {"enhancer": 1, "CTCF_binding_site": 1}}
    # _types renders as 'CTCF_binding_site:1,enhancer:1', so the first is CTCF.
    assert annotation_name_token(types, ["reg"]) == "CTCF_binding_site"


def test_name_token_needs_no_known_vocabulary():
    """An unheard-of feature type must work with no code change."""
    types = {"src": {"my_novel_element_v9": 3}}
    assert annotation_name_token(types, ["src"]) == "my_novel_element_v9"


def test_name_token_uses_annotate_source_order():
    """Source precedence is the user's --annotate order, not hardcoded."""
    types = {"first": {}, "second": {"EMAR": 1}, "third": {"enhancer": 1}}
    assert annotation_name_token(types, ["first", "second", "third"]) == "EMAR"
    # Reordering --annotate reorders precedence.
    assert annotation_name_token(types, ["third", "second"]) == "enhancer"


def test_name_token_is_none_without_overlap():
    assert annotation_name_token({"reg": {}}, ["reg"]) is None
    assert annotation_name_token({}, ["reg"]) is None


@pytest.mark.parametrize(
    "raw,expected",
    [("enhancer", "enhancer"), ("CTCF_binding_site", "CTCF_binding_site"),
     ("open chromatin region", "open_chromatin_region"),
     ('we"ird;type', "we_ird_type"), ("  spaced  ", "spaced"), ("!!!", "annotated")],
)
def test_sanitise_name_token(raw, expected):
    assert sanitise_name_token(raw) == expected


def test_names_get_feature_type_suffix(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    assert apply_annotation_names(regions, ["reg"]) == 1
    named = [c.consensus_transcript_name for c in regions
             if c.consensus_transcript_name]
    # chr1:12700-13000 hits enhancer ENSR1_CZ only.
    assert named == ["consensus_transcript_1-enhancer"]


def test_unannotated_loci_keep_bare_names(tmp_path):
    """A locus overlapping nothing must not gain a suffix."""
    import polars as pl

    rows = [{
        "region_id": "1:700000-700300", "chrom": "1", "start": 700001,
        "end": 700300, "class": "likely_novel_transcript",
        "dominant_strand": "+", "unique_fraction": 0.99,
        "dual_strand_fraction": 0.1, "profile_correlation": 0.0,
        "avg_depth": 20.0, "sample": s,
    } for s in ("A", "B")]

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    assert apply_annotation_names(regions, ["reg"]) == 0
    assert regions[0].consensus_transcript_name == "consensus_transcript_1"


def test_annotate_names_can_be_disabled(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)

    assert apply_annotation_names(regions, ["reg"], enabled=False) == 0
    assert [c.consensus_transcript_name for c in regions
            if c.consensus_transcript_name] == ["consensus_transcript_1"]


def test_suffix_keeps_the_consensus_transcript_prefix_matchable(tmp_path):
    """The point of a suffix: ^consensus_transcript_ still matches."""
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)
    apply_annotation_names(regions, ["reg"])

    names = [c.consensus_transcript_name for c in regions
             if c.consensus_transcript_name]
    assert names and all(n.startswith("consensus_transcript_") for n in names)


def test_decorated_names_stay_unique(tmp_path):
    """Two adjacent loci on the same feature type must not collide."""
    import polars as pl

    rows = []
    for s in ("A", "B"):
        for start in (12700, 13500):    # both inside enhancer ENSR1_CZ span
            rows.append({
                "region_id": f"1:{start}-{start + 300}", "chrom": "1",
                "start": start + 1, "end": start + 300,
                "class": "likely_novel_transcript", "dominant_strand": "+",
                "unique_fraction": 0.99, "dual_strand_fraction": 0.05,
                "profile_correlation": 0.0, "avg_depth": 20.0, "sample": s,
            })

    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(pl.DataFrame(rows), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)
    apply_annotation_names(regions, ["reg"])

    names = [c.consensus_transcript_name for c in regions
             if c.consensus_transcript_name]
    assert len(names) == 2
    assert len(set(names)) == 2
    assert all(n.endswith("-enhancer") for n in names)


def test_decorated_name_reaches_gtf_and_table(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2)
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    annotate_consensus(regions, indexes, cfg)
    apply_annotation_names(regions, ["reg"])

    out = os.path.join(str(tmp_path), "gtf.gtf")
    write_consensus_gtf(regions, out)
    with open(out) as fh:
        text = fh.read()
    assert 'transcript_id "consensus_transcript_1-enhancer"' in text
    assert 'gene_id "consensus_transcript_1-enhancer_gene"' in text

    df = consensus_to_dataframe(regions, overlay_labels=["reg"])
    assert "consensus_transcript_1-enhancer" in df["consensus_transcript_name"].to_list()


def test_summary_records_annotation_block(tmp_path):
    cfg = ConsensusConfig(tsvs=["dummy"], out_prefix="c", min_samples=2,
                          annotate=["reg=x.gff3"])
    regions = build_consensus(_consensus_rows(), cfg)
    indexes = load_sources([f"reg={_write(tmp_path, 'reg.gff3', ENSEMBL_GFF3)}"])
    stats = annotate_consensus(regions, indexes, cfg)

    out = os.path.join(str(tmp_path), "c.consensus_summary.json")
    summary = write_consensus_summary(
        cfg, regions, out, annotation_stats=stats, overlay_labels=["reg"]
    )
    with open(out) as fh:
        on_disk = json.load(fh)

    assert summary["annotation"]["sources"][0]["label"] == "reg"
    assert on_disk["annotation"]["sources"][0]["covered_bp_by_type"]["enhancer"] > 0
    assert CONS_NOVEL in on_disk["annotation"]["crosstab_passing"]["reg"][
        "by_consensus_class"
    ]
