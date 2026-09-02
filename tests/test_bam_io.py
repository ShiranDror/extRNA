"""Unit tests for gdna_rescue.bam_io.read_category (pysam-free).

read_category(read, cfg) is duck-typed: it never touches pysam, only reads
plain attributes off the ``read`` object. So we can exercise it on any
platform (including native Windows, where pysam has no wheels) using a
simple fake read stand-in instead of a real pysam.AlignedSegment.
"""

from dataclasses import dataclass

import pytest

from gdna_rescue.config import Config
from gdna_rescue.bam_io import read_category


@dataclass
class FakeRead:
    """Minimal stand-in for pysam.AlignedSegment covering the flags/fields
    read_category() (and _pair_disposition()) actually look at.

    Defaults describe a plain single-end primary read that would be "unique"
    under default settings: no flags set, MAPQ 255 (STAR's unique-mapping
    MAPQ), not part of a pair.
    """

    is_unmapped: bool = False
    is_qcfail: bool = False
    is_supplementary: bool = False
    is_duplicate: bool = False
    is_secondary: bool = False
    is_paired: bool = False
    mate_is_unmapped: bool = False
    is_proper_pair: bool = False
    mapping_quality: int = 255
    # Only used to document the TLEN trap case below; read_category() must
    # never consult this field.
    template_length: int = 0


@pytest.fixture
def cfg():
    return Config()


# --- Pair-geometry filtering (pair_filter="concordant", the default) ------


def test_half_mapped_pair_is_multi_despite_max_mapq(cfg):
    # Paired, mate unmapped, MAPQ 255 -> still "multi": a half-mapped mate
    # never counts as unique regardless of MAPQ.
    read = FakeRead(is_paired=True, mate_is_unmapped=True, mapping_quality=255)
    assert read_category(read, cfg) == "multi"


def test_discordant_pair_is_multi_despite_max_mapq(cfg):
    # Paired, mate mapped but not flagged as a proper pair, MAPQ 255 ->
    # "multi": discordant orientation/placement is noise regardless of MAPQ.
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=False,
        is_proper_pair=False,
        mapping_quality=255,
    )
    assert read_category(read, cfg) == "multi"


def test_proper_pair_with_huge_tlen_is_unique(cfg):
    # THE TRAP CASE: paired, proper-pair FLAG set, mate mapped, MAPQ 255, but
    # a huge genomic TLEN (80 kb) as would occur for an intron-spanning pair
    # from STAR. Must still be "unique".
    #
    # The filter has to key off the aligner's proper-pair FLAG, never a TLEN
    # cutoff: splice-aware aligners flag intron-spanning pairs as proper even
    # when genomic TLEN is 10-100 kb, so a TLEN cutoff would systematically
    # bias AGAINST long-intron novel transcripts -- exactly what this tool
    # exists to find.
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=False,
        is_proper_pair=True,
        mapping_quality=255,
        template_length=80000,
    )
    assert read_category(read, cfg) == "unique"


def test_single_end_high_mapq_is_unique_unchanged(cfg):
    read = FakeRead(is_paired=False, mapping_quality=255)
    assert read_category(read, cfg) == "unique"


def test_single_end_ignores_nonsense_mate_unmapped_flag(cfg):
    # is_paired=False with mate_is_unmapped=True is a nonsense combination
    # that can appear in malformed flags; pair geometry is only consulted
    # when is_paired is True, so this is still "unique".
    read = FakeRead(is_paired=False, mate_is_unmapped=True, mapping_quality=255)
    assert read_category(read, cfg) == "unique"


# --- pair_filter="off": recovers the old, pre-pair-filter behavior --------


def test_pair_filter_off_half_mapped_becomes_unique():
    cfg = Config(pair_filter="off")
    read = FakeRead(is_paired=True, mate_is_unmapped=True, mapping_quality=255)
    assert read_category(read, cfg) == "unique"


def test_pair_filter_off_discordant_becomes_unique():
    cfg = Config(pair_filter="off")
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=False,
        is_proper_pair=False,
        mapping_quality=255,
    )
    assert read_category(read, cfg) == "unique"


def test_pair_filter_off_low_mapq_primary_is_still_multi():
    cfg = Config(pair_filter="off")
    read = FakeRead(is_paired=True, is_proper_pair=True, mapping_quality=3)
    assert read_category(read, cfg) == "multi"


def test_pair_filter_off_duplicate_is_still_none():
    cfg = Config(pair_filter="off")
    read = FakeRead(is_duplicate=True, mapping_quality=255)
    assert read_category(read, cfg) is None


# --- Additional coverage ----------------------------------------------------


def test_half_mapped_duplicate_is_dropped_not_multi(cfg):
    # Duplicate drop happens before pair-geometry classification, so a
    # duplicate half-mapped read is None, not "multi".
    read = FakeRead(
        is_paired=True, mate_is_unmapped=True, is_duplicate=True, mapping_quality=255
    )
    assert read_category(read, cfg) is None


def test_half_mapped_secondary_is_multi_when_counted(cfg):
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=True,
        is_secondary=True,
        mapping_quality=255,
    )
    assert cfg.count_secondary is True
    assert read_category(read, cfg) == "multi"


def test_half_mapped_secondary_is_none_when_not_counted():
    cfg = Config(count_secondary=False)
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=True,
        is_secondary=True,
        mapping_quality=255,
    )
    assert read_category(read, cfg) is None


def test_low_mapq_discordant_pair_is_still_multi(cfg):
    # Already "multi" via MAPQ alone; pair-geometry classification doesn't
    # change the outcome.
    read = FakeRead(
        is_paired=True,
        mate_is_unmapped=False,
        is_proper_pair=False,
        mapping_quality=3,
    )
    assert read_category(read, cfg) == "multi"


def test_bogus_pair_filter_rejected_by_validate():
    cfg = Config(pair_filter="bogus")
    with pytest.raises(ValueError):
        cfg.validate()
