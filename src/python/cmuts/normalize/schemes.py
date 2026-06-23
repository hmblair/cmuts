"""Normalization schemes for reactivity profiles.

This module provides different normalization methods for converting raw
mutation rates into normalized reactivity values.

Normalization Schemes:
    RAW: No normalization, raw mutation rates
    UBR: Upper-bound reference percentile normalization (default)
    OUTLIER: 2-8% outlier-based normalization
"""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..internal import Opts, ProbingData

__all__ = [
    "NormScheme",
    "get_norm",
    "pooled_norm",
]


class NormScheme(Enum):
    """Normalization scheme for reactivity values."""

    RAW = 0
    UBR = 1
    OUTLIER = 2


def _get_norm_scheme(norm: str) -> NormScheme:
    """Determine normalization scheme from string."""
    norm_lower = norm.lower()
    if norm_lower == "raw":
        return NormScheme.RAW
    if norm_lower == "ubr":
        return NormScheme.UBR
    if norm_lower == "outlier":
        return NormScheme.OUTLIER

    raise ValueError(f"Unknown normalization scheme: {norm}. Use 'ubr', 'raw', or 'outlier'.")


def _get_norm_raw(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """No normalization."""
    return np.ones(1, dtype=data.reactivity.dtype)


# UBR (upper-bound reference) normalization parameters.
_UBR_READ_CUTOFF = 500  # per-reference read floor for "high-coverage" positions
_UBR_PERCENTILE = 90  # reactivity percentile used as the normalization factor


def _ubr_high_coverage(data: ProbingData) -> np.ndarray:
    """Reactivity values at masked, high-coverage positions of one dataset."""
    mask = np.asarray(data.mask)
    reads = np.asarray(data.reads)
    reactivity = np.asarray(data.reactivity)
    good_pos = mask & (reads > _UBR_READ_CUTOFF)[:, None]
    return reactivity[good_pos]


def _ubr_norm(datasets: list[ProbingData], opts: Opts) -> np.ndarray:
    """UBR factor: the percentile of high-coverage reactivity pooled across one
    or more datasets. Falls back to no normalization when no position clears the
    coverage floor.
    """
    vals = [v for v in (_ubr_high_coverage(d) for d in datasets) if v.size]
    if not vals:
        return _get_norm_raw(datasets[0], opts)
    pooled = np.concatenate([v.flatten() for v in vals])
    return np.asarray(np.percentile(pooled, _UBR_PERCENTILE))


def _get_norm_percentile(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """UBR percentile normalization: the 90th percentile of high-coverage
    positions, computed as a single dataset pooled through :func:`_ubr_norm`.
    """
    return _ubr_norm([data], opts)


def _get_norm_outlier(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """2-8% Normalization.

    From top 10%, ignore top 2% and divide by average of remaining 8%.
    """

    def norm_1d(
        arr: np.ndarray,
        low: float = 0.02,
        high: float = 0.1,
    ) -> float:
        arr = arr[~np.isnan(arr)]

        p2 = max(1, round(len(arr) * low) - 1)
        p10 = max(1, round(len(arr) * high) - 1)

        sarr = np.sort(arr)[::-1]
        if sarr.size == 0:
            m: float = 1.0
        else:
            m = float(np.mean(sarr[p2 : p10 + 1]))

        return m

    return np.apply_along_axis(norm_1d, axis=1, arr=data.reactivity)[:, None]


def get_norm(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """Get normalization factor based on options.

    Args:
        data: ProbingData containing reactivity values
        opts: Options specifying normalization scheme

    Returns:
        Normalization factor(s) as numpy array
    """
    scheme = _get_norm_scheme(opts.norm)

    if scheme == NormScheme.RAW:
        return _get_norm_raw(data, opts)
    if scheme == NormScheme.UBR:
        return _get_norm_percentile(data, opts)
    if scheme == NormScheme.OUTLIER:
        return _get_norm_outlier(data, opts)

    raise ValueError(f"Invalid normalization scheme {scheme}.")


def pooled_norm(
    datasets: list[ProbingData],
    opts: Opts,
) -> np.ndarray:
    """Compute a single normalization factor pooled across multiple datasets.

    Used to make reactivities directly comparable across experimental groups.
    For UBR, pools high-coverage reactivity values across all datasets and takes
    the percentile. For OUTLIER, averages the per-reference factors across
    datasets. For RAW, returns ones.

    Args:
        datasets: List of computed ProbingData objects (one per group).
        opts: Options specifying normalization scheme.

    Returns:
        Shared normalization factor as a numpy array.
    """
    if not datasets:
        return np.ones(1, dtype=np.float64)
    if len(datasets) == 1:
        return get_norm(datasets[0], opts)

    scheme = _get_norm_scheme(opts.norm)

    if scheme == NormScheme.RAW:
        return _get_norm_raw(datasets[0], opts)

    if scheme == NormScheme.UBR:
        return _ubr_norm(datasets, opts)

    if scheme == NormScheme.OUTLIER:
        norms = [_get_norm_outlier(data, opts) for data in datasets]
        return np.mean(norms, axis=0)

    raise ValueError(f"Invalid normalization scheme {scheme}.")
