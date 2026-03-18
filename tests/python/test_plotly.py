"""Tests for the Plotly plotting backend (cmuts.visualize.plotly)."""

from __future__ import annotations

import numpy as np
import pytest

go = pytest.importorskip("plotly.graph_objects")

from cmuts.internal import ProbingData  # noqa: E402
from cmuts.visualize.plotly import (  # noqa: E402
    plot_correlation,
    plot_coverage,
    plot_cumulative_reads,
    plot_examples,
    plot_heatmap,
    plot_mi,
    plot_multiple_examples,
    plot_profile,
    plot_profiles,
    plot_read_hist,
    plot_snr_scaling,
    plot_termination,
)

pytestmark = pytest.mark.no_external_dependencies


def _make_probing_data(
    nseq: int = 1,
    length: int = 12,
    with_pairwise: bool = False,
) -> ProbingData:
    x = np.linspace(0.05, 0.95, length, dtype=float)
    reactivity = np.tile(x, (nseq, 1))
    reads = np.arange(1, nseq + 1, dtype=float) * 100
    error = np.full((nseq, length), 0.05, dtype=float)
    snr = np.divide(reactivity, error, out=np.zeros_like(reactivity), where=error > 0).mean(-1)
    mask = np.ones((nseq, length), dtype=bool)
    heatmap = np.full((4, 7), 0.1, dtype=float)
    coverage = np.full(length, reads.mean(), dtype=float)
    terminations = np.tile(np.linspace(1, length, length, dtype=float), (nseq, 1))

    if with_pairwise:
        base = np.eye(length, dtype=float)[None, :, :]
        covariance = np.repeat(base, nseq, axis=0)
        mi = np.repeat(base * 0.1, nseq, axis=0)
    else:
        covariance = None
        mi = None

    data = ProbingData(
        sequences=None,
        reactivity=reactivity,
        reads=reads,
        error=error,
        snr=snr,
        mask=mask,
        heatmap=heatmap,
        coverage=coverage,
        terminations=terminations,
        pairs=None,
        probability=None,
        covariance=covariance,
        mi=mi,
    )
    data.norm = np.ones(nseq, dtype=float)
    return data


class TestPlotlyFigures:
    """Each plot function should return a go.Figure."""

    def test_plot_heatmap(self) -> None:
        data = _make_probing_data()
        fig = plot_heatmap(data.heatmap, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_profile_single(self) -> None:
        data = _make_probing_data()
        fig = plot_profile(data.reactivity[0], data.error[0], "sample")
        assert isinstance(fig, go.Figure)
        assert len(fig.data) == 2  # reactivity + error

    def test_plot_profile_no_error(self) -> None:
        data = _make_probing_data()
        fig = plot_profile(data.reactivity[0], name="sample")
        assert isinstance(fig, go.Figure)
        assert len(fig.data) == 1

    def test_plot_examples_single(self) -> None:
        data = _make_probing_data(nseq=1)
        fig = plot_examples(data.reactivity, data.error, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_examples_multi(self) -> None:
        data = _make_probing_data(nseq=10)
        fig = plot_examples(data.reactivity, data.error, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_multiple_examples(self) -> None:
        data = _make_probing_data(nseq=10)
        fig = plot_multiple_examples(data.reactivity)
        assert isinstance(fig, go.Figure)

    def test_plot_termination(self) -> None:
        data = _make_probing_data()
        fig = plot_termination(data.terminations, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_coverage(self) -> None:
        data = _make_probing_data()
        fig = plot_coverage(data.coverage, data.reads, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_read_hist(self) -> None:
        data = _make_probing_data(nseq=50)
        fig = plot_read_hist(data.reads, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_cumulative_reads(self) -> None:
        data = _make_probing_data(nseq=200)
        fig = plot_cumulative_reads(data.reads, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_profiles(self) -> None:
        r1 = np.linspace(0.1, 0.8, 10)
        r2 = np.linspace(0.2, 0.9, 10)
        fig = plot_profiles([r1, r2], ["A", "B"])
        assert isinstance(fig, go.Figure)
        assert len(fig.data) == 2

    def test_plot_correlation(self) -> None:
        data = _make_probing_data(nseq=1, with_pairwise=True)
        fig = plot_correlation(data.covariance[0], "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_mi(self) -> None:
        data = _make_probing_data(nseq=1, with_pairwise=True)
        fig = plot_mi(data.mi[0], "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_snr_scaling_with_nomod(self) -> None:
        mod = _make_probing_data()
        nomod = _make_probing_data()
        combined = _make_probing_data()
        fig = plot_snr_scaling(mod, nomod, combined, "sample")
        assert isinstance(fig, go.Figure)

    def test_plot_snr_scaling_without_nomod(self) -> None:
        mod = _make_probing_data()
        combined = _make_probing_data()
        fig = plot_snr_scaling(mod, None, combined, "sample")
        assert isinstance(fig, go.Figure)
