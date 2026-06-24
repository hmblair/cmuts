"""Unit tests for the normalization granularity dispatch (`schemes.normalization`).

These check the *semantics* of the two axes: the per-position factor must be
constant along the axis a flag pools over, and vary along the axis it splits on.
A 2-experiment x 2-reference grid is built so that references and experiments
have distinct value ranges, making their pooled factors differ.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from cmuts.normalize.schemes import _resolve, normalization


class _Data:
    """Minimal stand-in for ProbingData (the dispatch reads reactivity/reads)."""

    def __init__(self, reactivity: list[list[float]], reads: list[float]) -> None:
        self.reactivity = np.asarray(reactivity, dtype=float)
        self.reads = np.asarray(reads, dtype=float)
        self.sequences = None


def _factors(per_experiment: bool, per_reference: bool) -> list[np.ndarray]:
    # 2 experiments x 2 references x 4 positions. Reference 1 is an order of
    # magnitude hotter than reference 0, and experiment 1 hotter than 0, so each
    # pool yields a distinct factor. High coverage so UBR keeps every position.
    e0 = _Data([[1, 2, 3, 4], [10, 20, 30, 40]], [1000, 1000])
    e1 = _Data([[5, 6, 7, 8], [50, 60, 70, 80]], [1000, 1000])
    opts = SimpleNamespace(norm="ubr", per_experiment=per_experiment, per_reference=per_reference)
    return normalization([e0, e1], opts)  # one (R, P) factor array per experiment


def _constant_along_positions(factors: list[np.ndarray]) -> bool:
    """Every reference's factor is constant across positions (uniform scheme)."""
    return all(np.allclose(f, f[:, [0]]) for f in factors)


def test_global_norm_is_constant_everywhere():
    """No flags: a single factor over all experiments, references, and positions."""
    factors = _factors(per_experiment=False, per_reference=False)
    stacked = np.stack(factors)  # (E, R, P)
    assert np.allclose(stacked, stacked.flat[0])


def test_per_reference_constant_across_experiments_varies_across_references():
    factors = _factors(per_experiment=False, per_reference=True)
    assert _constant_along_positions(factors)
    # pooled across experiments -> reference r has the same factor in every experiment
    assert np.allclose(factors[0], factors[1])
    # split on reference -> different references get different factors
    assert not np.isclose(factors[0][0, 0], factors[0][1, 0])


def test_per_experiment_constant_across_references_varies_across_experiments():
    factors = _factors(per_experiment=True, per_reference=False)
    assert _constant_along_positions(factors)
    # pooled across references -> within an experiment every reference shares a factor
    for f in factors:
        assert np.allclose(f[0], f[1])
    # split on experiment -> different experiments get different factors
    assert not np.isclose(factors[0][0, 0], factors[1][0, 0])


def test_both_axes_give_independent_factors():
    factors = _factors(per_experiment=True, per_reference=True)
    assert _constant_along_positions(factors)
    # each (experiment, reference) cell is normalized on its own -> all four differ
    cells = [factors[e][r, 0] for e in range(2) for r in range(2)]
    assert len({round(float(c), 6) for c in cells}) == 4


def _shapemapper_boxplot(x: np.ndarray) -> float:
    """Reference port of ShapeMapper2's find_boxplot_factor (calc_quartile type 7
    == numpy's default linear percentile)."""
    x = np.sort(x[np.isfinite(x)])
    n = len(x)
    ten = n // 10
    q = 1.5 * abs(np.percentile(x, 25) - np.percentile(x, 75))
    cut = max(q, x[n - 1 - ten]) if n >= 100 else max(q, x[n - 1 - n // 20])
    kept = x[x < cut]
    return float(np.mean(kept[-ten:]))


def test_sm_shape_matches_shapemapper_boxplot():
    """The sm-shape scheme reproduces ShapeMapper2's SHAPE boxplot factor exactly."""
    scheme = _resolve("sm-shape")
    rng = np.random.default_rng(0)
    for n in (50, 207, 1000):
        arr = np.abs(rng.standard_exponential(n)) * 0.1
        arr[rng.integers(0, n, n // 20)] = np.nan  # scattered missing positions
        block = arr.reshape(1, n)  # one reference of n positions, high coverage
        mine = float(scheme.block_factor(block, np.array([1e6]), None)[0, 0])
        assert np.isclose(mine, _shapemapper_boxplot(arr))
