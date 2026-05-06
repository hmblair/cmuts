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


def _get_norm_percentile(
    data: ProbingData,
    opts: Opts,
) -> np.ndarray:
    """UBR percentile normalization.

    Uses the 90th percentile of high-coverage positions as the
    normalization factor.
    """
    READ_CUTOFF = 100
    PERCENTILE = 90

    _good_pos = data.mask & (data.reads > READ_CUTOFF)[:, None]
    high_reactivity = data.reactivity[_good_pos]

    if not high_reactivity.size:
        return _get_norm_raw(data, opts)

    return np.asarray(np.percentile(high_reactivity.flatten(), PERCENTILE))


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
        READ_CUTOFF = 100
        PERCENTILE = 90

        all_vals: list[np.ndarray] = []
        for data in datasets:
            mask = np.asarray(data.mask)
            reads = np.asarray(data.reads)
            reactivity = np.asarray(data.reactivity)
            good_pos = mask & (reads > READ_CUTOFF)[:, None]
            vals = reactivity[good_pos]
            if vals.size:
                all_vals.append(vals)

        if not all_vals:
            return _get_norm_raw(datasets[0], opts)

        pooled = np.concatenate(all_vals).flatten()
        return np.asarray(np.percentile(pooled, PERCENTILE))

    if scheme == NormScheme.OUTLIER:
        norms = [_get_norm_outlier(data, opts) for data in datasets]
        return np.mean(norms, axis=0)

    raise ValueError(f"Invalid normalization scheme {scheme}.")
