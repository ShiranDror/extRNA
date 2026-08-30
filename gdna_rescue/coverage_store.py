"""Sparse, multi-channel per-base coverage store (HDF5).

The detect pipeline computes dense per-strand coverage per chromosome anyway;
instead of throwing it away it hands the arrays here to be *sparsified* (only the
non-zero, non-masked bases are kept) and written to ``{prefix}.coverage.h5``.
``merge_candidates`` later reads small windows back out to draw per-transcript
coverage plots. numpy + h5py only -- NO pysam -- so both the write side (Linux/
WSL, where the BAMs are read) and the read side (anywhere, including native
Windows) work.

Channels are strand-resolved so strand is first-class evidence on the plots, and
split by read category so a reviewer can see how much of a pileup is genuine
unique signal versus PCR-duplicate or low-quality/multimapper coverage:

    unique_plus / unique_minus : dedup, MAPQ>=min_mapq, primary (defines regions)
    dup_plus    / dup_minus    : reads that would be unique but are duplicate-flagged
    multi_plus  / multi_minus  : low-MAPQ / secondary ("low quality") coverage

Coordinates are 0-based. The store holds only *unannotated* coverage: positions
masked as annotated during discovery are dropped, matching the "all non-masked
coverage" contract, and appear as zero when a window is read back.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

# Canonical channel order. Kept in one place so writer, reader and plots agree.
CHANNELS = (
    "unique_plus", "unique_minus",
    "dup_plus", "dup_minus",
    "multi_plus", "multi_minus",
)

FORMAT_TAG = "gdna_rescue_coverage_store_v1"


def _require_h5py():
    try:
        import h5py  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without h5py
        raise ImportError(
            "h5py is required for the coverage store but is not installed. "
            "Install it with `pip install h5py` or `conda install -c conda-forge "
            "h5py`."
        ) from exc
    return h5py


def sparsify(
    channels: Dict[str, np.ndarray],
    mask_plus: Optional[np.ndarray] = None,
    mask_minus: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Reduce dense per-base channel arrays to their non-zero, non-masked bases.

    ``channels`` maps names in :data:`CHANNELS` to dense int arrays over one
    chromosome (missing channels are treated as all-zero). ``mask_plus`` /
    ``mask_minus`` (bool arrays, same length) mark annotated positions to remove
    from the ``*_plus`` / ``*_minus`` channels respectively, so the store keeps
    only unannotated coverage.

    Returns ``(pos, sparse)`` where ``pos`` is the sorted 0-based array of bases
    that are non-zero in at least one channel, and ``sparse`` maps every channel
    name to an int32 array aligned to ``pos``. Empty bases are never stored.
    """
    length = 0
    masked: Dict[str, np.ndarray] = {}
    for name in CHANNELS:
        arr = channels.get(name)
        if arr is None:
            continue
        length = max(length, int(arr.shape[0]))
        m = mask_plus if name.endswith("_plus") else mask_minus
        if m is not None:
            m = m[: arr.shape[0]]
            arr = np.where(m, 0, arr)
        masked[name] = arr

    if length == 0:
        empty = np.zeros(0, dtype=np.int64)
        return empty, {n: np.zeros(0, dtype=np.int32) for n in CHANNELS}

    nz = np.zeros(length, dtype=bool)
    for arr in masked.values():
        nz[: arr.shape[0]] |= arr != 0
    pos = np.nonzero(nz)[0].astype(np.int64)

    sparse: Dict[str, np.ndarray] = {}
    for name in CHANNELS:
        arr = masked.get(name)
        if arr is None:
            sparse[name] = np.zeros(pos.shape[0], dtype=np.int32)
        else:
            vals = np.zeros(pos.shape[0], dtype=np.int32)
            in_range = pos < arr.shape[0]
            vals[in_range] = arr[pos[in_range]].astype(np.int32)
            sparse[name] = vals
    return pos, sparse


def write_coverage_store(
    path: str,
    per_chrom: Dict[str, Tuple[np.ndarray, Dict[str, np.ndarray]]],
    *,
    sample_name: Optional[str] = None,
    chrom_sizes: Optional[Dict[str, int]] = None,
    strandedness: Optional[str] = None,
    mask_mode: str = "positional",
    compression: str = "gzip",
) -> int:
    """Write sparsified per-chrom channels to an HDF5 store.

    ``per_chrom`` maps chrom -> ``(pos, sparse)`` as returned by :func:`sparsify`.
    Returns the number of chromosomes actually written (those with any coverage).
    """
    h5py = _require_h5py()
    chrom_sizes = chrom_sizes or {}
    n_written = 0
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = FORMAT_TAG
        h5.attrs["sample"] = sample_name or ""
        h5.attrs["channels"] = list(CHANNELS)
        h5.attrs["strandedness"] = strandedness or "unknown"
        h5.attrs["mask_mode"] = mask_mode
        root = h5.create_group("chroms")
        for chrom, (pos, sparse) in per_chrom.items():
            if pos is None or pos.shape[0] == 0:
                continue
            g = root.create_group(str(chrom))
            g.attrs["length"] = int(chrom_sizes.get(chrom, 0))
            g.create_dataset("pos", data=pos.astype(np.int64), compression=compression)
            for name in CHANNELS:
                g.create_dataset(
                    name, data=sparse[name].astype(np.int32), compression=compression
                )
            n_written += 1
    return n_written


class CoverageStore:
    """Read-only accessor for a coverage store; slices dense windows on demand.

    Per-chromosome sparse arrays are cached on first access, so drawing many
    transcripts on the same chromosome reads the file once.
    """

    def __init__(self, path: str):
        h5py = _require_h5py()
        self.path = path
        self._h5 = h5py.File(path, "r")
        self.sample = str(self._h5.attrs.get("sample", "") or "")
        self.strandedness = str(self._h5.attrs.get("strandedness", "unknown") or "unknown")
        self.mask_mode = str(self._h5.attrs.get("mask_mode", "") or "")
        self._cache: Dict[str, Optional[Tuple[np.ndarray, Dict[str, np.ndarray]]]] = {}

    def evict_all(self) -> None:
        """Drop all cached chromosome arrays (frees the resident window data)."""
        self._cache.clear()

    def close(self) -> None:
        self._h5.close()

    def __enter__(self) -> "CoverageStore":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    @property
    def chroms(self):
        grp = self._h5.get("chroms")
        return set(grp.keys()) if grp is not None else set()

    def _load(self, chrom: str):
        if chrom not in self._cache:
            key = f"chroms/{chrom}"
            if key not in self._h5:
                self._cache[chrom] = None
            else:
                g = self._h5[key]
                pos = g["pos"][:]
                data = {n: g[n][:] for n in CHANNELS}
                self._cache[chrom] = (pos, data)
        return self._cache[chrom]

    def window(self, chrom: str, start: int, end: int) -> Dict[str, np.ndarray]:
        """Dense per-channel coverage over ``[start, end)`` (0-based half-open).

        Bases with no stored coverage (uncovered or masked as annotated) are 0.
        """
        n = max(0, int(end) - int(start))
        dense = {name: np.zeros(n, dtype=np.int32) for name in CHANNELS}
        loaded = self._load(chrom)
        if loaded is None or n == 0:
            return dense
        pos, data = loaded
        lo = int(np.searchsorted(pos, start, side="left"))
        hi = int(np.searchsorted(pos, end, side="left"))
        if hi <= lo:
            return dense
        idx = (pos[lo:hi] - start).astype(np.int64)
        for name in CHANNELS:
            dense[name][idx] = data[name][lo:hi]
        return dense
