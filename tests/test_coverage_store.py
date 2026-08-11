"""Tests for the sparse coverage store (numpy + h5py, no pysam)."""

from __future__ import annotations

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from gdna_rescue.coverage_store import (  # noqa: E402
    CHANNELS,
    CoverageStore,
    sparsify,
    write_coverage_store,
)


def _dense(length, spikes):
    """Build a dense int32 array of ``length`` with {pos: value} spikes."""
    arr = np.zeros(length, dtype=np.int32)
    for pos, val in spikes.items():
        arr[pos] = val
    return arr


def test_sparsify_drops_zero_and_masked_bases():
    length = 20
    channels = {
        "unique_plus": _dense(length, {2: 5, 10: 7}),
        "unique_minus": _dense(length, {2: 3}),
        "multi_plus": _dense(length, {15: 4}),
    }
    mask_plus = np.zeros(length, dtype=bool)
    mask_plus[10] = True  # annotated -> must be dropped from *_plus only

    pos, sparse = sparsify(channels, mask_plus=mask_plus, mask_minus=None)

    # Non-zero, non-masked union of positions: 2 and 15 (10 was masked on +).
    assert list(pos) == [2, 15]
    # Every channel is present and aligned to pos.
    for name in CHANNELS:
        assert sparse[name].shape[0] == pos.shape[0]
    # Values line up; masked base 10 removed, base 2 kept on both strands.
    idx = {int(p): i for i, p in enumerate(pos)}
    assert sparse["unique_plus"][idx[2]] == 5
    assert sparse["unique_minus"][idx[2]] == 3
    assert sparse["multi_plus"][idx[15]] == 4


def test_write_and_slice_roundtrip(tmp_path):
    length = 100
    up = _dense(length, {10: 8, 11: 9, 50: 4})
    um = _dense(length, {10: 2})
    dp = _dense(length, {11: 5})
    channels = {"unique_plus": up, "unique_minus": um, "dup_plus": dp}
    pos, sparse = sparsify(channels)

    path = str(tmp_path / "s.coverage.h5")
    n = write_coverage_store(
        path, {"1": (pos, sparse)},
        sample_name="S1", chrom_sizes={"1": length}, strandedness="reverse",
        mask_mode="stranded",
    )
    assert n == 1

    with CoverageStore(path) as store:
        assert store.sample == "S1"
        assert store.strandedness == "reverse"
        assert "1" in store.chroms
        win = store.window("1", 8, 14)  # bases 8..13
        # window is dense incl. zeros; base 10 -> offset 2, base 11 -> offset 3.
        assert win["unique_plus"][2] == 8
        assert win["unique_plus"][3] == 9
        assert win["unique_minus"][2] == 2
        assert win["dup_plus"][3] == 5
        assert win["unique_plus"][0] == 0  # uncovered base is zero
        # A chromosome not in the store reads back as all zeros.
        assert store.window("MT", 0, 5)["unique_plus"].sum() == 0


def test_empty_store(tmp_path):
    path = str(tmp_path / "e.coverage.h5")
    empty_pos = np.zeros(0, dtype=np.int64)
    empty_sparse = {c: np.zeros(0, dtype=np.int32) for c in CHANNELS}
    n = write_coverage_store(path, {"1": (empty_pos, empty_sparse)},
                             sample_name="E")
    assert n == 0  # nothing written for a chrom with no coverage
    with CoverageStore(path) as store:
        assert store.window("1", 0, 10)["unique_plus"].sum() == 0
