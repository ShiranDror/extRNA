"""Small dependency-light helpers: logging, interval maths, coordinate utils.

Deliberately free of pysam so it can be imported by the pure-numpy modules and
by the test suite on any platform.
"""

from __future__ import annotations

import logging
import sys
from typing import List, Tuple

import numpy as np


def get_logger(verbose: bool = False) -> logging.Logger:
    """Return a configured module logger writing to stderr."""
    logger = logging.getLogger("gdna_rescue")
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("[%(asctime)s] %(levelname)s %(message)s", "%H:%M:%S")
        )
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger


def flip_strand(strand: str) -> str:
    """Return the opposite strand character."""
    if strand == "+":
        return "-"
    if strand == "-":
        return "+"
    return strand  # '.' stays '.'


def find_true_runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Return contiguous [start, end) runs where a boolean mask is True.

    Uses a diff of the padded mask so it is O(n) and vectorised.
    """
    if mask.size == 0:
        return []
    m = mask.astype(np.int8)
    # Pad with zeros so runs at the very ends are detected.
    diff = np.diff(np.concatenate(([0], m, [0])))
    starts = np.flatnonzero(diff == 1)
    ends = np.flatnonzero(diff == -1)
    return list(zip(starts.tolist(), ends.tolist()))


def merge_runs_with_gap(
    runs: List[Tuple[int, int]],
    max_gap: int,
    blocked_mask: np.ndarray | None = None,
) -> List[Tuple[int, int]]:
    """Merge adjacent [start, end) runs separated by <= ``max_gap`` bases.

    If ``blocked_mask`` is given, a gap is only bridged when it contains no
    blocked position. This is how we avoid merging two candidate intervals
    across an annotated feature that happens to sit between them.
    """
    if not runs:
        return []
    runs = sorted(runs)
    merged = [list(runs[0])]
    for start, end in runs[1:]:
        prev = merged[-1]
        gap = start - prev[1]
        bridgeable = gap <= max_gap
        if bridgeable and blocked_mask is not None and gap > 0:
            if blocked_mask[prev[1]:start].any():
                bridgeable = False
        if bridgeable:
            prev[1] = max(prev[1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson correlation that returns 0.0 when either vector has no variance.

    scipy.stats.pearsonr raises / warns on constant input; for coverage
    profiles a flat strand simply means "no correlation information", which we
    encode as 0.0.
    """
    if a.size < 2 or b.size < 2:
        return 0.0
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    sa = a.std()
    sb = b.std()
    if sa == 0.0 or sb == 0.0:
        return 0.0
    r = float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))
    # Guard against tiny floating point excursions outside [-1, 1].
    return max(-1.0, min(1.0, r))


def binary_entropy(p: float) -> float:
    """Shannon entropy (bits) of a Bernoulli(p); 1.0 = balanced, 0.0 = pure."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return float(-p * np.log2(p) - (1 - p) * np.log2(1 - p))
