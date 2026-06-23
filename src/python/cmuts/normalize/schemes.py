"""Normalization schemes for reactivity profiles.

A scheme converts raw mutation rates into a reactivity scale factor. Schemes are
small classes registered in a ``name -> instance`` table; both the CLI ``--norm``
choices and the dispatch read from it, so adding a scheme is a single localized
change (subclass :class:`Scheme`, decorate with :func:`register`).

Built-in schemes:
    raw     : no normalization
    ubr     : upper-bound reference -- 90th percentile of high-coverage positions
    outlier : 2-8% outlier-based, per-reference normalization
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..internal import Opts, ProbingData

__all__ = [
    "Scheme",
    "get_norm",
    "pooled_norm",
    "register",
    "scheme_names",
]


# UBR (upper-bound reference) normalization parameters.
_UBR_READ_CUTOFF = 500  # per-reference read floor for "high-coverage" positions
_UBR_PERCENTILE = 90  # reactivity percentile used as the normalization factor


# ---------------------------------------------------------------------------
# Scheme base class and registry
# ---------------------------------------------------------------------------


class Scheme(ABC):
    """A normalization scheme: maps probing data to a reactivity scale factor.

    Subclasses set ``name`` and implement :meth:`single`. The default
    :meth:`pooled` averages each group's factor, which is correct for
    per-reference and constant schemes; override it when groups must share a
    factor computed jointly (e.g. UBR pools raw values before the percentile).
    """

    name: str

    @abstractmethod
    def single(self, data: ProbingData, opts: Opts) -> np.ndarray:
        """Normalization factor for one dataset."""

    def pooled(self, datasets: list[ProbingData], opts: Opts) -> np.ndarray:
        """Shared factor across groups (default: average the per-group factors)."""
        return np.mean([self.single(d, opts) for d in datasets], axis=0)


_SCHEMES: dict[str, Scheme] = {}


def register(scheme_cls: type[Scheme]) -> type[Scheme]:
    """Class decorator: register an instance of ``scheme_cls`` under its name."""
    instance = scheme_cls()
    _SCHEMES[instance.name] = instance
    return scheme_cls


def scheme_names() -> list[str]:
    """Names of all registered schemes (drives the CLI ``--norm`` choices)."""
    return list(_SCHEMES)


def _resolve(name: str) -> Scheme:
    try:
        return _SCHEMES[name.lower()]
    except KeyError:
        raise ValueError(
            f"Unknown normalization scheme: {name}. Use one of {scheme_names()}."
        ) from None


# ---------------------------------------------------------------------------
# Built-in schemes
# ---------------------------------------------------------------------------


@register
class RawScheme(Scheme):
    """No normalization."""

    name = "raw"

    def single(self, data: ProbingData, opts: Opts) -> np.ndarray:
        return np.ones(1, dtype=np.asarray(data.reactivity).dtype)

    def pooled(self, datasets: list[ProbingData], opts: Opts) -> np.ndarray:
        return self.single(datasets[0], opts)


def _ubr_high_coverage(data: ProbingData) -> np.ndarray:
    """Reactivity values at masked, high-coverage positions of one dataset."""
    mask = np.asarray(data.mask)
    reads = np.asarray(data.reads)
    reactivity = np.asarray(data.reactivity)
    good_pos = mask & (reads > _UBR_READ_CUTOFF)[:, None]
    return reactivity[good_pos]


@register
class UBRScheme(Scheme):
    """Upper-bound reference: the percentile of high-coverage reactivity pooled
    across datasets. Falls back to no normalization when no position clears the
    coverage floor.
    """

    name = "ubr"

    def _norm(self, datasets: list[ProbingData]) -> np.ndarray:
        vals = [v for v in (_ubr_high_coverage(d) for d in datasets) if v.size]
        if not vals:
            return np.ones(1, dtype=np.asarray(datasets[0].reactivity).dtype)
        pooled = np.concatenate([v.flatten() for v in vals])
        return np.asarray(np.percentile(pooled, _UBR_PERCENTILE))

    def single(self, data: ProbingData, opts: Opts) -> np.ndarray:
        return self._norm([data])

    def pooled(self, datasets: list[ProbingData], opts: Opts) -> np.ndarray:
        return self._norm(datasets)


@register
class OutlierScheme(Scheme):
    """2-8% normalization: from the top 10%, drop the top 2% and divide by the
    mean of the remaining 8%. Produces a per-reference factor.
    """

    name = "outlier"

    def single(self, data: ProbingData, opts: Opts) -> np.ndarray:
        def norm_1d(arr: np.ndarray, low: float = 0.02, high: float = 0.1) -> float:
            arr = arr[~np.isnan(arr)]
            p2 = max(1, round(len(arr) * low) - 1)
            p10 = max(1, round(len(arr) * high) - 1)
            sarr = np.sort(arr)[::-1]
            if sarr.size == 0:
                return 1.0
            return float(np.mean(sarr[p2 : p10 + 1]))

        return np.apply_along_axis(norm_1d, axis=1, arr=data.reactivity)[:, None]


# ---------------------------------------------------------------------------
# Public dispatch API
# ---------------------------------------------------------------------------


def get_norm(data: ProbingData, opts: Opts) -> np.ndarray:
    """Normalization factor for a single dataset under ``opts.norm``."""
    return _resolve(opts.norm).single(data, opts)


def pooled_norm(datasets: list[ProbingData], opts: Opts) -> np.ndarray:
    """Shared normalization factor across datasets under ``opts.norm``.

    Used to make reactivities directly comparable across experimental groups.
    Empty input yields ones; a single dataset is equivalent to :func:`get_norm`.
    """
    if not datasets:
        return np.ones(1, dtype=np.float64)
    scheme = _resolve(opts.norm)
    if len(datasets) == 1:
        return scheme.single(datasets[0], opts)
    return scheme.pooled(datasets, opts)
