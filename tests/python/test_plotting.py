"""
General integration tests for cmuts plotting entry points.

These tests use small synthetic ProbingData objects and a temporary output
directory to verify that the plotting APIs run headlessly and produce files.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from cmuts.internal import ProbingData
from cmuts.visualize.plotting import main, plot_all, plot_profiles, plot_snr_scaling

pytestmark = pytest.mark.no_external_dependencies


def _make_probing_data(
    nseq: int = 1,
    length: int = 12,
    with_pairwise: bool = False,
) -> ProbingData:
    """Build a small but representative ProbingData object for plotting."""
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


class TestPlotAll:
    """Smoke tests for the main plot generator."""

    def test_plot_all_single_sequence_creates_standard_figures(
        self,
        tmp_path: Path,
    ) -> None:
        data = _make_probing_data(nseq=1)

        plot_all(data, "sample", str(tmp_path))

        created = {path.name for path in tmp_path.glob("*.png")}
        assert "sample-heatmap.png" in created
        assert "sample-profile.png" in created
        assert "sample-termination.png" in created
        assert "sample-coverage.png" in created

    def test_plot_all_multi_sequence_creates_summary_figures(self, tmp_path: Path) -> None:
        data = _make_probing_data(nseq=120)

        plot_all(data, "sample", str(tmp_path))

        created = {path.name for path in tmp_path.glob("*.png")}
        assert "sample-heatmap.png" in created
        assert "sample-examples.png" in created
        assert "sample-reads-hist.png" in created
        assert "sample-cumulative-reads.png" in created

    def test_plot_all_handles_pairwise_data(
        self,
        tmp_path: Path,
    ) -> None:
        data = _make_probing_data(nseq=1, with_pairwise=True)

        plot_all(data, "sample", str(tmp_path))

        created = {path.name for path in tmp_path.glob("*.png")}
        assert "sample-mutual-information.png" in created
        assert "sample-correlation.png" in created

    def test_plot_all_indexes_pairwise_outputs_for_multiple_sequences(self, tmp_path: Path) -> None:
        data = _make_probing_data(nseq=3, with_pairwise=True)

        plot_all(data, "sample", str(tmp_path))

        created = {path.name for path in tmp_path.glob("*.png")}
        assert "sample-1-mutual-information.png" in created
        assert "sample-2-mutual-information.png" in created
        assert "sample-3-mutual-information.png" in created
        assert "sample-1-correlation.png" in created
        assert "sample-2-correlation.png" in created
        assert "sample-3-correlation.png" in created


class TestStandalonePlotting:
    """Smoke tests for additional public plotting entry points."""

    def test_plot_profiles_overlays_multiple_profiles(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.chdir(tmp_path)

        reactivities = [
            np.linspace(0.1, 0.8, 10, dtype=float),
            np.linspace(0.2, 0.9, 10, dtype=float),
        ]
        plot_profiles(reactivities, ["A", "B"], str(tmp_path))

        assert (tmp_path / "profile.png").exists()

    def test_plot_snr_scaling_writes_output(self, tmp_path: Path) -> None:
        mod = _make_probing_data(nseq=1)
        nomod = _make_probing_data(nseq=1)
        combined = _make_probing_data(nseq=1)

        plot_snr_scaling(mod, nomod, combined, "sample", str(tmp_path))

        assert (tmp_path / "sample-snr-scaling.png").exists()


class TestPlotCli:
    """Smoke test for the cmuts-plot CLI entry point."""

    def test_main_reads_hdf5_and_creates_plots(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        data = _make_probing_data(nseq=1, with_pairwise=True)
        input_path = tmp_path / "input.h5"
        out_dir = tmp_path / "figures"
        data.save("sample", str(input_path))

        monkeypatch.setattr(
            sys,
            "argv",
            ["cmuts-plot", str(input_path), "--group", "sample", "--out", str(out_dir)],
        )

        main()

        created = {path.name for path in out_dir.glob("*.png")}
        assert "sample-heatmap.png" in created
        assert "sample-profile.png" in created
        assert "sample-correlation.png" in created
        assert "sample-mutual-information.png" in created
        assert "Plots saved to" in capsys.readouterr().out
