"""Unit tests for region metrics and the rule-based classifier (pysam-free)."""

import numpy as np
import pytest

from gdna_rescue.config import Config
from gdna_rescue.classify import (
    compute_region_metrics,
    classify_region,
    LIKELY_GDNA,
    LIKELY_NOVEL,
    LIKELY_MULTIMAPPER,
    POSSIBLE_BIDIRECTIONAL,
    is_kept,
)
from tests.generate_test_data import make_archetype_coverage


@pytest.fixture
def cfg():
    return Config()


def _label(plus, minus, cfg):
    m = compute_region_metrics(plus, minus, cfg)
    label, reason, flags = classify_region(m, cfg)
    return label, m, reason


def test_gdna_symmetric_is_flagged(cfg):
    plus, minus = make_archetype_coverage()["gdna_symmetric"]
    label, m, reason = _label(plus, minus, cfg)
    assert label == LIKELY_GDNA, reason
    assert not is_kept(label)
    assert m.profile_correlation >= cfg.gdna_min_profile_correlation
    assert m.dual_strand_fraction >= cfg.gdna_min_dual_strand_fraction


def test_flat_symmetric_gdna_is_flagged(cfg):
    # Randomly-fragmented gDNA -> flat, uncorrelated coverage on BOTH strands.
    # Correlation is ~0 here, so this exercises the flatness branch.
    rng = np.random.default_rng(11)
    n = 1500
    plus = np.clip(30 + rng.normal(0, 3, n), 0, None).astype(np.int32)
    minus = np.clip(30 + rng.normal(0, 3, n), 0, None).astype(np.int32)
    label, m, reason = _label(plus, minus, cfg)
    assert m.profile_correlation < cfg.gdna_min_profile_correlation  # not correlated
    assert label == LIKELY_GDNA, reason
    assert not is_kept(label)


def test_novel_single_strand_is_rescued(cfg):
    plus, minus = make_archetype_coverage()["novel_single"]
    label, m, reason = _label(plus, minus, cfg)
    assert label == LIKELY_NOVEL, reason
    assert is_kept(label)
    assert m.dominant_strand_fraction >= cfg.novel_min_dominant_strand_fraction
    assert m.dominant_strand == "+"


def test_bidirectional_asymmetric_is_ambiguous(cfg):
    plus, minus = make_archetype_coverage()["bidir_asymmetric"]
    label, m, reason = _label(plus, minus, cfg)
    assert label == POSSIBLE_BIDIRECTIONAL, reason
    assert is_kept(label)  # not gDNA -> rescued


def test_reverse_of_novel_picks_minus_strand(cfg):
    plus, minus = make_archetype_coverage()["novel_single"]
    # Swap strands: now the minus strand should dominate.
    label, m, reason = _label(minus, plus, cfg)
    assert label == LIKELY_NOVEL
    assert m.dominant_strand == "-"


def test_multimapper_artifact_is_flagged_and_not_rescued(cfg):
    n = 800
    # Unique coverage alone would look like a clean single-strand novel locus...
    plus = np.full(n, 12, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    # ...but the locus is swamped by multimapped reads -> artifact.
    multi = np.full(n, 100, dtype=np.int32)
    m = compute_region_metrics(plus, minus, cfg, multi_cov=multi)
    label, reason, flags = classify_region(m, cfg)
    assert label == LIKELY_MULTIMAPPER, reason
    assert not is_kept(label)
    assert m.unique_fraction < cfg.min_unique_fraction

    # Without the multimapper pile-up the same region is a novel transcript.
    m2 = compute_region_metrics(plus, minus, cfg)
    assert classify_region(m2, cfg)[0] == LIKELY_NOVEL
    assert m2.unique_fraction == 1.0


def test_empty_and_flat_profiles_do_not_crash(cfg):
    z = np.zeros(300, dtype=np.int32)
    m = compute_region_metrics(z, z, cfg)
    label, reason, _ = classify_region(m, cfg)
    # All-zero -> dominant fraction defaults to 1.0 -> treated as single strand.
    assert label in (LIKELY_NOVEL, POSSIBLE_BIDIRECTIONAL)
    assert m.profile_correlation == 0.0


def test_non_flat_uncorrelated_symmetric_is_not_gdna(cfg):
    # Symmetric, broad, dual-strand, but structured AND anti-correlated:
    # plus ramps up while minus ramps down. High CV (not flat) and negative r
    # -> pattern is not gDNA-consistent, so it must NOT be called gDNA.
    n = 800
    plus = np.linspace(1, 40, n).astype(np.int32)
    minus = np.linspace(40, 1, n).astype(np.int32)
    m = compute_region_metrics(plus, minus, cfg)
    label, reason, flags = classify_region(m, cfg)
    assert not flags["gdna_profile_correlated"]
    assert not flags["gdna_both_strands_flat"]
    assert label != LIKELY_GDNA, reason


def test_metric_ranges(cfg):
    plus, minus = make_archetype_coverage()["gdna_symmetric"]
    m = compute_region_metrics(plus, minus, cfg)
    assert 0.0 <= m.dual_strand_fraction <= 1.0
    assert 0.0 <= m.covered_fraction <= 1.0
    assert 0.5 <= m.dominant_strand_fraction <= 1.0
    assert -1.0 <= m.profile_correlation <= 1.0
    assert 0.0 <= m.strand_entropy <= 1.0
