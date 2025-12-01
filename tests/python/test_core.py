"""
Integration tests for cmuts core.

This module contains both fixed edge-case tests and random fuzzing tests
to verify cmuts produces correct mutation counts.
"""
from __future__ import annotations

import random

import pytest

from helpers import TestParams, run_test


# =============================================================================
# Fixed Edge Case Tests
# =============================================================================


class TestEdgeCases:
    """Fixed, deterministic test cases for edge conditions."""

    def test_minimal_single_read(self):
        """Single read on single reference - simplest case."""
        params = TestParams(
            length=50,
            references=1,
            queries=1,
            seed=42,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_short_reads_at_min_length(self):
        """Reads exactly at minimum length threshold."""
        params = TestParams(
            length=20,
            references=1,
            queries=50,
            min_length=15,
            max_length=25,
            seed=123,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_high_mapq_filter(self):
        """High MAPQ threshold should filter most reads."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            min_mapq=50,
            seed=456,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_high_phred_filter(self):
        """High PHRED threshold should mask low-quality bases."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            min_phred=30,
            seed=789,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_quality_window_averaging(self):
        """Quality window averaging over multiple bases."""
        params = TestParams(
            length=100,
            references=1,
            queries=100,
            quality_window=5,
            min_phred=20,
            seed=101,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_no_mismatches(self):
        """Counting only insertions and deletions."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            no_mismatches=True,
            seed=202,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_no_insertions(self):
        """Counting only mismatches and deletions."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            no_insertions=True,
            seed=303,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_no_deletions(self):
        """Counting only mismatches and insertions."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            no_deletions=True,
            seed=404,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_mismatches_only(self):
        """Counting only mismatches (no indels)."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            no_insertions=True,
            no_deletions=True,
            seed=505,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_collapse_adjacent_modifications(self):
        """Collapse modifications within distance threshold."""
        params = TestParams(
            length=100,
            references=2,
            queries=200,
            collapse=5,
            seed=606,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_large_collapse_distance(self):
        """Large collapse distance (most modifications collapsed)."""
        params = TestParams(
            length=50,
            references=1,
            queries=100,
            collapse=20,
            seed=707,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_max_indel_length_filter(self):
        """Filter out long indels."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            max_indel_length=3,
            seed=808,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_multiple_references(self):
        """Multiple reference sequences."""
        params = TestParams(
            length=80,
            references=10,
            queries=50,
            seed=909,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_many_reads_per_reference(self):
        """High coverage scenario."""
        params = TestParams(
            length=50,
            references=2,
            queries=500,
            seed=1010,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_combined_filters(self):
        """Multiple filters applied simultaneously."""
        params = TestParams(
            length=100,
            references=3,
            queries=200,
            min_mapq=10,
            min_phred=15,
            min_length=30,
            max_length=150,
            max_indel_length=5,
            quality_window=3,
            collapse=3,
            seed=1111,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()


class TestCramFormat:
    """Tests specifically for CRAM file format handling."""

    def test_cram_basic(self):
        """Basic CRAM file processing."""
        params = TestParams(
            length=100,
            references=2,
            queries=100,
            use_cram=True,
            seed=2001,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_cram_with_filters(self):
        """CRAM with quality filters."""
        params = TestParams(
            length=100,
            references=3,
            queries=150,
            min_mapq=15,
            min_phred=20,
            use_cram=True,
            seed=2002,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()


# =============================================================================
# Random Fuzzing Tests
# =============================================================================


def random_params(seed: int | None = None) -> TestParams:
    """Generate random test parameters for fuzzing."""
    rng = random.Random(seed)

    length = rng.randint(10, 250)
    return TestParams(
        length=length,
        references=rng.randint(1, 32),
        queries=rng.randint(1, 512),
        threads=rng.randint(1, 4),
        chunk_size=rng.randint(16, 256),
        min_mapq=rng.randint(0, 30),
        min_phred=rng.randint(0, 20),
        quality_window=rng.randint(0, min(length, 20)),
        min_length=rng.randint(2, length),
        max_length=rng.randint(length, 2 * length),
        max_indel_length=rng.randint(0, length),
        collapse=rng.randint(0, length),
        no_mismatches=rng.random() < 0.2,
        no_insertions=rng.random() < 0.2,
        no_deletions=rng.random() < 0.2,
        use_cram=rng.random() < 0.3,
        seed=seed,
    )


class TestRandomFuzzing:
    """Random parameter fuzzing tests."""

    @pytest.mark.parametrize("iteration", range(10))
    def test_random_bam(self, iteration: int):
        """Random parameters with BAM format."""
        # Use iteration as part of seed for reproducibility
        seed = 10000 + iteration
        params = random_params(seed)
        params.use_cram = False

        result = run_test(params)
        assert result.counts_match, f"seed={seed}\n{result.summary()}"

    @pytest.mark.parametrize("iteration", range(5))
    def test_random_cram(self, iteration: int):
        """Random parameters with CRAM format."""
        seed = 20000 + iteration
        params = random_params(seed)
        params.use_cram = True

        result = run_test(params)
        assert result.counts_match, f"seed={seed}\n{result.summary()}"


# =============================================================================
# Stress Tests (marked slow)
# =============================================================================


@pytest.mark.slow
class TestStress:
    """Stress tests with larger data (marked slow)."""

    def test_large_reference(self):
        """Very long reference sequence."""
        params = TestParams(
            length=1000,
            references=1,
            queries=100,
            seed=30001,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_many_references(self):
        """Many reference sequences."""
        params = TestParams(
            length=50,
            references=100,
            queries=20,
            seed=30002,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    def test_high_coverage(self):
        """Very high coverage."""
        params = TestParams(
            length=100,
            references=5,
            queries=2000,
            seed=30003,
        )
        result = run_test(params)
        assert result.counts_match, result.summary()

    @pytest.mark.parametrize("iteration", range(50))
    def test_extended_fuzzing(self, iteration: int):
        """Extended fuzzing with more iterations."""
        seed = 40000 + iteration
        params = random_params(seed)
        result = run_test(params)
        assert result.counts_match, f"seed={seed}\n{result.summary()}"
