"""
Integration tests for cmuts core --tokenize and the --token-map option.

These tests write a small FASTA and run `cmuts core --tokenize` with no
alignment files, which writes the per-base `sequence` dataset directly. They
verify the integer tokens emitted for each base, including custom mappings
supplied via --token-map.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

# Tokenization needs only the cmuts CLI, not samtools.
pytestmark = pytest.mark.no_external_dependencies

_CMUTS = shutil.which("cmuts")


@pytest.fixture(autouse=True)
def _require_cmuts() -> None:
    if _CMUTS is None:
        pytest.skip("cmuts not found in PATH")


def _tokenize(
    tmp_path: Path,
    sequences: list[str],
    token_map: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Write a FASTA and run `cmuts core --tokenize`, returning the process."""
    fasta = tmp_path / "ref.fa"
    fasta.write_text("".join(f">r{i}\n{s}\n" for i, s in enumerate(sequences)))
    out = tmp_path / "out.h5"
    args = [
        "cmuts",
        "core",
        "--tokenize",
        "--fasta",
        str(fasta),
        "--output",
        str(out),
    ]
    if token_map is not None:
        args += ["--token-map", token_map]
    return subprocess.run(args, capture_output=True, text=True)


def _read_sequence(tmp_path: Path) -> np.ndarray:
    with h5py.File(tmp_path / "out.h5", "r") as f:
        return f["sequence"][:]


def test_default_token_map_is_identity(tmp_path: Path) -> None:
    """With no --token-map, bases map to A=0, C=1, G=2, U=3."""
    proc = _tokenize(tmp_path, ["ACGU"])
    assert proc.returncode == 0, proc.stderr
    seq = _read_sequence(tmp_path)
    assert seq[0].tolist() == [0, 1, 2, 3]


def test_custom_token_map(tmp_path: Path) -> None:
    """--token-map remaps each base to an arbitrary integer."""
    proc = _tokenize(tmp_path, ["ACGU"], token_map="4,5,6,7")
    assert proc.returncode == 0, proc.stderr
    seq = _read_sequence(tmp_path)
    assert seq[0].tolist() == [4, 5, 6, 7]


def test_t_and_u_share_the_fourth_token(tmp_path: Path) -> None:
    """T and U both map to the fourth value (they collapse internally)."""
    proc = _tokenize(tmp_path, ["ACGU", "ACGT"], token_map="4,5,6,7")
    assert proc.returncode == 0, proc.stderr
    seq = _read_sequence(tmp_path)
    assert seq[0].tolist() == [4, 5, 6, 7]
    assert seq[1].tolist() == [4, 5, 6, 7]


def test_token_map_permutation(tmp_path: Path) -> None:
    """A reordering map relabels bases without changing the alphabet size."""
    proc = _tokenize(tmp_path, ["ACGUT"], token_map="3,2,1,0")
    assert proc.returncode == 0, proc.stderr
    seq = _read_sequence(tmp_path)
    assert seq[0].tolist() == [3, 2, 1, 0, 0]


def test_custom_token_map_with_uneven_lengths(tmp_path: Path) -> None:
    """Mapping is applied per base across references of differing lengths."""
    sequences = ["ACGU", "AC"]
    proc = _tokenize(tmp_path, sequences, token_map="4,5,6,7")
    assert proc.returncode == 0, proc.stderr
    seq = _read_sequence(tmp_path)
    assert seq[0, : len(sequences[0])].tolist() == [4, 5, 6, 7]
    assert seq[1, : len(sequences[1])].tolist() == [4, 5]


def test_too_few_tokens_is_rejected(tmp_path: Path) -> None:
    """A map without exactly four values fails and writes no dataset."""
    proc = _tokenize(tmp_path, ["ACGU"], token_map="4,5,6")
    assert proc.returncode != 0
    assert "token-map" in proc.stderr.lower()
    with h5py.File(tmp_path / "out.h5", "r") as f:
        assert "sequence" not in f


def test_negative_token_is_rejected(tmp_path: Path) -> None:
    """A negative token fails with a clear message and writes no dataset."""
    proc = _tokenize(tmp_path, ["ACGU"], token_map="4,5,6,-7")
    assert proc.returncode != 0
    assert "non-negative" in proc.stderr.lower()
    with h5py.File(tmp_path / "out.h5", "r") as f:
        assert "sequence" not in f


def test_out_of_range_token_is_rejected(tmp_path: Path) -> None:
    """A token outside the int8 output range fails and writes no dataset."""
    proc = _tokenize(tmp_path, ["ACGU"], token_map="4,5,6,9999")
    assert proc.returncode != 0
    assert "range" in proc.stderr.lower()
    with h5py.File(tmp_path / "out.h5", "r") as f:
        assert "sequence" not in f
