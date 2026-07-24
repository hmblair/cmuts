"""Focused arithmetic tests for the deletion/coverage denominator fix in
`_data_from_counts`.

Base-call coverage (`coverage = sum(counts[..., :IX_DEL], axis=(2, 3))`)
excludes deletion outcomes entirely -- a deletion only ever increments the
`(rbase, IX_DEL)` cell, never a match/mismatch cell. When `--deletions` is
on, deletion counts are added to the modification numerator but, before this
fix, never to the denominator -- so `modifications / coverage` could exceed
1 at deletion-heavy positions, silently absorbed by `_reactivity`'s
unconditional `clip(0, 1)`.

Each test here hand-builds a single-reference, single-position counts-1d
tensor (shape (1, 1, 4, 7): true base x [A, C, G, T, del, ins, term]) with a
known true reference base, so the expected reactivity/error/mask/coverage
values can be computed by hand and checked exactly -- deliberately narrower
than an end-to-end integration test, to pin down the arithmetic itself.
"""

from __future__ import annotations

import dask.array as da
import numpy as np
import pytest

from cmuts.internal import DataGroups, Opts, _data_from_counts

pytestmark = pytest.mark.no_external_dependencies

IX_DEL = 4
IX_INS = 5
TRUE_BASE_A = 0  # counts[..., TRUE_BASE_A, :] is populated; all other true-base slices are zero


def _opts(cutoff: int, ins: bool, dels: bool) -> Opts:
    return Opts(
        mod=DataGroups(None),
        nomod=DataGroups(None),
        cutoff=cutoff,
        ins=ins,
        dels=dels,
        norm="raw",
        blank=(0, 0),
        clip=(None, None),
        sig=0.05,
    )


def _single_position_counts(match: int, mismatch: int, deletions: int, insertions: int) -> da.Array:
    """(1, 1, 4, 7) tensor: one reference, one position, true base = A.
    `mismatch` all goes to the A->G outcome for simplicity."""
    row = np.zeros(7)
    row[TRUE_BASE_A] = match  # match count (diagonal)
    row[2] = mismatch  # A -> G
    row[IX_DEL] = deletions
    row[IX_INS] = insertions
    counts = np.zeros((1, 1, 4, 7))
    counts[0, 0, TRUE_BASE_A, :] = row
    return da.from_array(counts)


def _run(counts: da.Array, cutoff: int, ins: bool, dels: bool):
    opts = _opts(cutoff, ins, dels)
    lengths = da.array([1])
    reactivity, reads, err, snr, mask, heatmap, mc, term, *_ = _data_from_counts(counts, None, opts, lengths)
    return {
        "reactivity": float(np.asarray(reactivity)[0, 0]),
        "error": float(np.asarray(err)[0, 0]),
        "mask": bool(np.asarray(mask)[0, 0]),
        "reads": float(np.asarray(reads)[0]),
        "mean_coverage": float(np.asarray(mc)[0]),
    }


class TestUnambiguousDeletionArithmetic:
    """80 matches, 10 mismatches, 10 deletions -- the correct event rate
    pairing mismatches+deletions in both numerator and denominator is
    (10+10)/(80+10+10) = 0.20, not the legacy (10+10)/(80+10)."""

    def test_deletions_on_matches_the_paired_event_rate(self):
        counts = _single_position_counts(match=80, mismatch=10, deletions=10, insertions=0)
        result = _run(counts, cutoff=0, ins=False, dels=True)
        assert result["reactivity"] == pytest.approx(0.20)
        assert result["reads"] == pytest.approx(100)
        assert result["mean_coverage"] == pytest.approx(100)
        assert result["error"] == pytest.approx(np.sqrt(0.20 * 0.80 / 100))

    def test_deletions_off_excludes_deletions_from_both_sides(self):
        """--no-deletions: deletions contribute to neither numerator nor
        denominator -- legacy behavior, unaffected by this patch."""
        counts = _single_position_counts(match=80, mismatch=10, deletions=10, insertions=0)
        result = _run(counts, cutoff=0, ins=False, dels=False)
        assert result["reactivity"] == pytest.approx(10 / 90)
        assert result["reads"] == pytest.approx(90)
        assert result["mean_coverage"] == pytest.approx(90)


class TestNoRegressionWithoutDeletions:
    def test_zero_deletions_is_a_noop(self):
        """With no deletion mass present, --deletions on/off must agree --
        there's nothing for the fix to change."""
        counts = _single_position_counts(match=80, mismatch=10, deletions=0, insertions=0)
        with_dels = _run(counts, cutoff=0, ins=False, dels=True)
        without_dels = _run(counts, cutoff=0, ins=False, dels=False)
        assert with_dels == pytest.approx(without_dels)
        assert with_dels["reactivity"] == pytest.approx(10 / 90)


class TestCutoffBoundary:
    def test_position_admitted_only_after_deletion_inclusive_denominator(self):
        """90 base calls (below a 95 cutoff) + 10 deletions = 100 (above
        it). This position is masked out under the legacy base-call-only
        denominator and admitted under the fix -- a real, if usually small,
        behavior change worth covering explicitly."""
        counts = _single_position_counts(match=80, mismatch=10, deletions=10, insertions=0)

        admitted_with_fix = _run(counts, cutoff=95, ins=False, dels=True)
        assert admitted_with_fix["mask"] is True

        # Legacy-equivalent: base-call-only coverage is 90, below the cutoff.
        legacy_equivalent = _run(counts, cutoff=95, ins=False, dels=False)
        assert legacy_equivalent["mask"] is False


class TestInsertionPlusDeletionNotBoundedByOne:
    """The fix corrects the deletion-specific gap; it does NOT make the
    overall ratio bounded by 1 in general. With insertions also enabled,
    the numerator can still exceed the (deletion-inclusive) denominator, and
    the clip in _reactivity remains load-bearing for that case. This test
    documents that limitation rather than hiding it."""

    def test_insertions_can_still_push_reactivity_to_the_clip(self):
        counts = _single_position_counts(match=80, mismatch=10, deletions=10, insertions=100)
        result = _run(counts, cutoff=0, ins=True, dels=True)
        # modifications = 10 (mismatch) + 100 (insertions) + 10 (deletions) = 120
        # coverage_denom = 90 + 10 = 100 -- insertions are not added to it
        # raw ratio = 1.2, clipped to 1.0 with error 0.
        assert result["reactivity"] == pytest.approx(1.0)
        assert result["error"] == pytest.approx(0.0)
