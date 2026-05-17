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


def _run_normalize(
    counts_h5: Path,
    fasta: Path,
    out: Path,
    mod: list[str],
    extra_args: list[str] | None = None,
) -> subprocess.CompletedProcess:
    """Run the cmuts-normalize script as a subprocess."""
    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parents[2] / "src" / "python" / "cmuts-normalize"),
        str(counts_h5),
        "--fasta",
        str(fasta),
        "--mod",
        *mod,
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
            assert "reactivity" in f
            assert "reads" in f
            assert "error" in f


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

    def test_group_flag(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--group flag places output under a group."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        out = tmp_path / "out.h5"
        result = _run_normalize(
            counts_h5,
            fasta,
            out,
            mod=groups,
            extra_args=["--group", "mygroup"],
        )
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

        with h5py.File(out, "r") as f:
            assert "mygroup" in f
            assert "reactivity" in f["mygroup"]


class TestNormalizeMultiGroup:
    """Tests for --groups multi-group mode."""

    def _write_toml(self, path: Path, entries: list[dict]) -> None:
        lines: list[str] = []
        for e in entries:
            lines.append("[[group]]")
            lines.append(f'name = "{e["name"]}"')
            mods = ", ".join(f'"{m}"' for m in e["mod"])
            lines.append(f"mod = [{mods}]")
            if e.get("nomod"):
                nomods = ", ".join(f'"{m}"' for m in e["nomod"])
                lines.append(f"nomod = [{nomods}]")
            lines.append("")
        path.write_text("\n".join(lines))

    def _run_groups(
        self,
        counts_h5: Path,
        fasta: Path,
        out: Path,
        groups_toml: Path,
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess:
        cmd = [
            sys.executable,
            str(Path(__file__).resolve().parents[2] / "src" / "python" / "cmuts-normalize"),
            str(counts_h5),
            "--fasta",
            str(fasta),
            "--groups",
            str(groups_toml),
            "-o",
            str(out),
            "--overwrite",
        ]
        if extra_args:
            cmd.extend(extra_args)
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60, cwd=out.parent)

    def test_two_groups_shared_norm(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Two groups produce a single HDF5 with both as top-level groups."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        assert groups, "No count groups found in test HDF5"
        # Reuse the single available dataset for both groups — sufficient to
        # exercise the multi-group code path without needing distinct counts.
        ds = groups[0]

        toml_path = tmp_path / "groups.toml"
        self._write_toml(
            toml_path,
            [
                {"name": "g1", "mod": [ds]},
                {"name": "g2", "mod": [ds]},
            ],
        )

        out = tmp_path / "out.h5"
        result = self._run_groups(counts_h5, fasta, out, toml_path)
        assert result.returncode == 0, f"stderr:\n{result.stderr}"
        assert out.exists()

        with h5py.File(out, "r") as f:
            assert "g1" in f and "g2" in f
            assert "reactivity" in f["g1"]
            assert "reactivity" in f["g2"]

    def test_independent_norm_flag(self, test_data: tuple[Path, Path], tmp_path: Path):
        """--independent-norm runs without error."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)
        assert groups
        ds = groups[0]

        toml_path = tmp_path / "groups.toml"
        self._write_toml(
            toml_path,
            [
                {"name": "g1", "mod": [ds]},
                {"name": "g2", "mod": [ds]},
            ],
        )

        out = tmp_path / "out.h5"
        result = self._run_groups(counts_h5, fasta, out, toml_path, ["--independent-norm"])
        assert result.returncode == 0, f"stderr:\n{result.stderr}"

    def test_groups_and_mod_mutually_exclusive(self, test_data: tuple[Path, Path], tmp_path: Path):
        """Passing both --groups and --mod is rejected."""
        counts_h5, fasta = test_data
        groups = _find_mod_groups(counts_h5)

        toml_path = tmp_path / "groups.toml"
        self._write_toml(toml_path, [{"name": "g1", "mod": [groups[0]]}])

        out = tmp_path / "out.h5"
        cmd = [
            sys.executable,
            str(Path(__file__).resolve().parents[2] / "src" / "python" / "cmuts-normalize"),
            str(counts_h5),
            "--fasta",
            str(fasta),
            "--groups",
            str(toml_path),
            "--mod",
            groups[0],
            "-o",
            str(out),
            "--overwrite",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
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
