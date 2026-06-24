"""
Integration tests for the cmuts-normalize CLI script.

Tests that the wrapper script correctly constructs Opts from CLI arguments
and invokes compute_reactivity. Uses synthetic HDF5 count data generated
by `cmuts generate` + `cmuts core`.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import h5py
import pytest
from helpers import TestParams, TestRunner


@pytest.fixture
def test_data(tmp_path: Path) -> tuple[Path, Path]:
    """Generate a small counts HDF5 + FASTA via cmuts generate + cmuts core."""
    params = TestParams(
        length=50,
        references=2,
        queries=100,
        seed=99,
    )
    runner = TestRunner(params)
    runner._work_dir = tmp_path
    runner.generate_test_data()
    runner.run_cmuts()
    return runner.cmuts_h5_path, runner.fasta_path


_SCRIPT = Path(__file__).resolve().parents[2] / "src" / "python" / "cmuts-normalize"


def _run_normalize(
    counts_h5: Path,
    fasta: Path,
    out: Path,
    mod: list[str],
    *,
    name: str = "exp",
    nomod: list[str] | None = None,
    extra_args: list[str] | None = None,
) -> subprocess.CompletedProcess:
    """Run cmuts-normalize for a single experiment named ``name``."""
    spec = [name, "mod=" + ",".join(mod)]
    if nomod:
        spec.append("nomod=" + ",".join(nomod))
    cmd = [
        sys.executable,
        str(_SCRIPT),
        str(counts_h5),
        "--fasta",
        str(fasta),
        "--experiment",
        *spec,
        "-o",
        str(out),
        "--overwrite",
    ]
    if extra_args:
        cmd.extend(extra_args)
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60, cwd=out.parent)


def _find_mod_groups(h5_path: Path) -> list[str]:
    """Find group paths containing counts-1d datasets."""
    groups = []
    with h5py.File(h5_path, "r") as f:

        def visitor(name: str, obj: object) -> None:
            if name.endswith("counts-1d") and hasattr(obj, "shape"):
                # Parent group is everything before /counts-1d
                groups.append(name.rsplit("/counts-1d", 1)[0])

        f.visititems(visitor)
    return groups


class TestNormalizeSmoke:
    """Basic smoke tests: does the script run and produce output?"""

    def test_runs_successfully(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Script exits 0 with valid inputs."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        assert groups, "No count groups found in test HDF5"

        out = tmp_path / "out.h5"
        result = _run_normalize(counts_h5, fasta, out, mod=groups)
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_produces_output_file(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Script creates an HDF5 output file with expected datasets."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(counts_h5, fasta, out, mod=groups)
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        assert out.exists()

        with h5py.File(out, "r") as f:
            assert "reactivity" in f["exp"]
            assert "reads" in f["exp"]
            assert "error" in f["exp"]


class TestNormalizeOpts:
    """Tests that CLI flags are correctly wired to Opts fields."""

    def test_outlier_normalization(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--norm outlier produces output."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(counts_h5, fasta, out, mod=groups, extra_args=["--norm", "outlier"])
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        assert out.exists()

    def test_sm_dms_normalization(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--norm sm-dms (ShapeMapper-style per-nucleotide) produces output.

        Exercises the sequence-aware path: the scheme needs the per-position base
        identity, which is attached from the FASTA when the counts file carries
        no tokenized sequence.
        """
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(counts_h5, fasta, out, mod=groups, extra_args=["--norm", "sm-dms"])
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        assert out.exists()
        with h5py.File(out, "r") as f:
            assert "reactivity" in f["exp"]

    def test_clip_flags(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--clip-below and --clip-above flags work with explicit thresholds."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            fasta,
            out,
            mod=groups,
            extra_args=["--clip-below", "0", "--clip-above", "1"],
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_blank_flags(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--blank-5p and --blank-3p flags work."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            fasta,
            out,
            mod=groups,
            extra_args=["--blank-5p", "3", "--blank-3p", "3"],
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_no_insertions(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--no-insertions flag works."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            fasta,
            out,
            mod=groups,
            extra_args=["--no-insertions"],
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_no_deletions(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--no-deletions flag works."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            fasta,
            out,
            mod=groups,
            extra_args=["--no-deletions"],
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_experiment_name_is_output_group(self, test_data: tuple[Path, Path], tmp_path: Path):
        """The experiment name becomes the output HDF5 group."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(counts_h5, fasta, out, mod=groups, name="mygroup")
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

        with h5py.File(out, "r") as f:
            assert "mygroup" in f
            assert "reactivity" in f["mygroup"]

    def test_per_reference_norm(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--per-reference-norm runs and produces output."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5, fasta, out, mod=groups, extra_args=["--per-reference-norm"]
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        with h5py.File(out, "r") as f:
            assert "reactivity" in f["exp"]


class TestNormalizeMultiExperiment:
    """Tests for multiple --experiment flags."""

    def _run(
        self,
        counts_h5: Path,
        fasta: Path,
        out: Path,
        experiments: list[tuple[str, list[str]]],
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess:
        cmd = [sys.executable, str(_SCRIPT), str(counts_h5), "--fasta", str(fasta)]
        for name, mod in experiments:
            cmd += ["--experiment", name, "mod=" + ",".join(mod)]
        cmd += ["-o", str(out), "--overwrite"]
        if extra_args:
            cmd.extend(extra_args)
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60, cwd=out.parent)

    def test_two_experiments(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Two experiments produce a single HDF5 with both as top-level groups."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        assert groups, "No count groups found in test HDF5"
        # Reuse the single available dataset for both experiments -- enough to
        # exercise the multi-experiment path without distinct counts.
        ds = groups[0]

        out = tmp_path / "out.h5"
        result = self._run(counts_h5, fasta, out, [("g1", [ds]), ("g2", [ds])])
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        assert out.exists()

        with h5py.File(out, "r") as f:
            assert "g1" in f and "g2" in f
            assert "reactivity" in f["g1"]
            assert "reactivity" in f["g2"]

    def test_per_experiment_norm(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--per-experiment-norm runs without error."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        ds = groups[0]

        out = tmp_path / "out.h5"
        result = self._run(
            counts_h5, fasta, out, [("g1", [ds]), ("g2", [ds])], ["--per-experiment-norm"]
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_duplicate_names_rejected(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Two experiments with the same name are rejected."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        ds = groups[0]

        out = tmp_path / "out.h5"
        result = self._run(counts_h5, fasta, out, [("g1", [ds]), ("g1", [ds])])
        assert result.returncode != 0


class TestNormalizeErrors:
    """Tests that the script fails gracefully on bad inputs."""

    def test_missing_fasta(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Script fails when FASTA file doesn't exist."""
        counts_h5, _ = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            tmp_path / "nonexistent.fasta",
            out,
            mod=groups,
        )
        assert result.returncode != 0

    def test_missing_input_file(self, tmp_path: Path):
        """Script fails when input HDF5 doesn't exist."""
        out = tmp_path / "out.h5"
        result = _run_normalize(
            tmp_path / "nonexistent.h5",
            tmp_path / "ref.fasta",
            out,
            mod=["sample"],
        )
        assert result.returncode != 0
