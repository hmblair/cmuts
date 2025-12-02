"""
Test utilities for cmuts integration tests.

This module provides fixtures for generating test data and running the cmuts
pipeline, enabling both fixed edge-case tests and random fuzzing.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np


@lru_cache(maxsize=1)
def has_mpi_support() -> bool:
    """Check if cmuts was built with MPI support."""
    try:
        result = subprocess.run(
            ["cmuts", "core", "--version"],
            capture_output=True,
            text=True,
        )
        return "MPI" in result.stdout
    except Exception:
        return False


@dataclass
class TestParams:  # noqa: pytest collection warning - not a test class
    """Parameters for test data generation and cmuts execution."""

    # Sequence parameters
    length: int = 100
    references: int = 1
    queries: int = 100

    # Threading
    threads: int = 1
    chunk_size: int = 128

    # Quality filters
    min_mapq: int = 0
    min_phred: int = 0
    quality_window: int = 1

    # Length filters
    min_length: int = 10
    max_length: int = 200
    max_indel_length: int = 10

    # Collapse distance
    collapse: int = 0

    # Modification flags
    no_mismatches: bool = False
    no_insertions: bool = False
    no_deletions: bool = False

    # File format
    use_cram: bool = False

    # Seed for reproducibility (None = random)
    seed: int | None = None

    def to_generate_args(self) -> list[str]:
        """Convert to cmuts-generate-tests arguments."""
        args = [
            "--length", str(self.length),
            "--queries", str(self.queries),
            "--references", str(self.references),
            "--min-mapq", str(self.min_mapq),
            "--min-phred", str(self.min_phred),
            "--min-length", str(self.min_length),
            "--max-length", str(self.max_length),
            "--max-indel-length", str(self.max_indel_length),
            "--quality-window", str(self.quality_window),
            "--collapse", str(self.collapse),
        ]
        if self.no_mismatches:
            args.append("--no-mismatches")
        if self.no_insertions:
            args.append("--no-insertions")
        if self.no_deletions:
            args.append("--no-deletions")
        if self.seed is not None:
            args.extend(["--seed", str(self.seed)])
        return args

    def to_cmuts_args(self) -> list[str]:
        """Convert to cmuts core arguments."""
        # Limit threads to 1 if MPI support is not available
        threads = self.threads if has_mpi_support() else 1
        args = [
            "--threads", str(threads),
            "--min-mapq", str(self.min_mapq),
            "--min-phred", str(self.min_phred),
            "--quality-window", str(self.quality_window),
            "--min-length", str(self.min_length),
            "--max-length", str(self.max_length),
            "--max-indel-length", str(self.max_indel_length),
            "--collapse", str(self.collapse),
            "--chunk-size", str(self.chunk_size),
            "--disable-ambiguous",
        ]
        if self.no_mismatches:
            args.append("--no-mismatches")
        if self.no_insertions:
            args.append("--no-insertions")
        if self.no_deletions:
            args.append("--no-deletions")
        return args


@dataclass
class TestResult:
    """Results from running a cmuts test."""

    cmuts_counts: np.ndarray
    expected_counts: np.ndarray
    params: TestParams
    work_dir: Path

    @property
    def difference(self) -> np.ndarray:
        """Element-wise difference between cmuts and expected."""
        return self.cmuts_counts - self.expected_counts

    @property
    def max_absolute_error(self) -> float:
        """Maximum absolute difference in any count."""
        return float(np.max(np.abs(self.difference)))

    @property
    def total_absolute_error(self) -> float:
        """Sum of all absolute differences."""
        return float(np.sum(np.abs(self.difference)))

    @property
    def mismatched_positions(self) -> int:
        """Number of positions where counts differ."""
        return int(np.sum(self.difference != 0))

    @property
    def total_coverage(self) -> float:
        """Total coverage in expected data."""
        return float(self.expected_counts.sum())

    @property
    def counts_match(self) -> bool:
        """Check if all counts match exactly."""
        return np.allclose(self.cmuts_counts, self.expected_counts, rtol=0, atol=1e-6)

    def summary(self) -> str:
        """Human-readable summary of test results."""
        lines = [
            f"Counts match: {self.counts_match}",
            f"Max absolute error: {self.max_absolute_error:.6f}",
            f"Total absolute error: {self.total_absolute_error:.6f}",
            f"Mismatched positions: {self.mismatched_positions}",
            f"Total coverage: {self.total_coverage:.0f}",
        ]
        return "\n".join(lines)


class TestRunner:
    """Manages test execution in a temporary directory."""

    def __init__(self, params: TestParams, keep_files: bool = False):
        self.params = params
        self.keep_files = keep_files
        self._work_dir: Path | None = None

    def __enter__(self) -> "TestRunner":
        self._work_dir = Path(tempfile.mkdtemp(prefix="cmuts_test_"))
        return self

    def __exit__(self, *args) -> None:
        if not self.keep_files and self._work_dir:
            shutil.rmtree(self._work_dir, ignore_errors=True)

    @property
    def work_dir(self) -> Path:
        if self._work_dir is None:
            raise RuntimeError("TestRunner must be used as context manager")
        return self._work_dir

    @property
    def fasta_path(self) -> Path:
        return self.work_dir / "seq.fasta"

    @property
    def sam_path(self) -> Path:
        return self.work_dir / "aln.sam"

    @property
    def alignment_path(self) -> Path:
        ext = "cram" if self.params.use_cram else "bam"
        return self.work_dir / f"aln.{ext}"

    @property
    def expected_h5_path(self) -> Path:
        return self.work_dir / "expected.h5"

    @property
    def cmuts_h5_path(self) -> Path:
        return self.work_dir / "cmuts.h5"

    def generate_test_data(self) -> None:
        """Generate synthetic test data."""
        args = ["cmuts-generate-tests"] + self.params.to_generate_args() + [
            "--out-fasta", str(self.fasta_path),
            "--out-sam", str(self.sam_path),
            "--out-h5", str(self.expected_h5_path),
        ]
        subprocess.run(args, check=True, capture_output=True)

    def convert_alignment(self) -> None:
        """Convert SAM to BAM/CRAM using samtools."""
        bam_tmp = self.work_dir / "aln_tmp.bam"
        threads = str(self.params.threads)

        if self.params.use_cram:
            # Sort then convert to CRAM
            subprocess.run(
                ["samtools", "sort", f"-@{threads}", "-o", str(bam_tmp), str(self.sam_path)],
                check=True, capture_output=True
            )
            subprocess.run(
                ["samtools", "view", "-T", str(self.fasta_path), "-C",
                 "--output-fmt-option", "version=3.0", "-o", str(self.alignment_path), str(bam_tmp)],
                check=True, capture_output=True
            )
        else:
            # Add MD tags and sort
            calmd_result = subprocess.run(
                ["samtools", "calmd", "-b", f"-@{threads}", str(self.sam_path), str(self.fasta_path)],
                check=True, capture_output=True
            )
            bam_tmp.write_bytes(calmd_result.stdout)
            subprocess.run(
                ["samtools", "sort", f"-@{threads}", "-o", str(self.alignment_path), str(bam_tmp)],
                check=True, capture_output=True
            )

    def run_cmuts(self) -> None:
        """Run cmuts core on the test data."""
        args = ["cmuts", "core"] + self.params.to_cmuts_args() + [
            "--fasta", str(self.fasta_path),
            "--output", str(self.cmuts_h5_path),
            str(self.alignment_path),
        ]
        subprocess.run(args, check=True, capture_output=True)

    def run(self) -> TestResult:
        """Execute full test pipeline and return results."""
        self.generate_test_data()
        self.convert_alignment()
        self.run_cmuts()

        def find_counts_dataset(f: h5py.File) -> str:
            """Find the counts-1d dataset path in an HDF5 file."""
            result = []
            f.visititems(lambda n, o: result.append(n) if n.endswith("counts-1d") else None)
            if not result:
                datasets = []
                f.visititems(lambda n, o: datasets.append(f"{n}: {o.shape}") if hasattr(o, "shape") else None)
                raise KeyError(f"No counts-1d dataset found. Available: {datasets}")
            return result[0]

        with h5py.File(self.cmuts_h5_path, "r") as f:
            dataset_path = find_counts_dataset(f)
            cmuts_counts = f[dataset_path][..., :-1]

        with h5py.File(self.expected_h5_path, "r") as f:
            dataset_path = find_counts_dataset(f)
            expected_counts = f[dataset_path][:]

        return TestResult(
            cmuts_counts=cmuts_counts,
            expected_counts=expected_counts,
            params=self.params,
            work_dir=self.work_dir,
        )


def run_test(params: TestParams, keep_files: bool = False) -> TestResult:
    """Run a single test with given parameters."""
    runner = TestRunner(params, keep_files=keep_files)
    runner.__enter__()
    try:
        result = runner.run()
        # Keep files for debugging failed CRAM tests
        if not result.counts_match and params.use_cram:
            runner.keep_files = True
            print(f"\nDEBUG: CRAM test failed, files kept at: {runner.work_dir}")
        return result
    finally:
        runner.__exit__(None, None, None)
