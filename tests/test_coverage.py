"""Coverage at a band of zero against samtools depth.

A band of zero pins the marginal to the CIGAR as written. A read without
run-length tags then covers exactly the positions its aligned bases pair, which
samtools depth counts, once the reads cmuts hmm always rejects are removed,
along with the reads whose CIGAR the model cannot follow as written.

A read with the tags admits cDNAs with more or fewer bases than the read shows,
and a base the read does not show can pair with a base the CIGAR deletes. Such
a read covers at most the positions its aligned bases pair or delete, which
samtools depth counts when told to count deletions.
"""

import numpy as np

from alignments import drop_adjacent_indels, strip_tags
from filters import UNFILTERED
from oracle import rows_by_name, samtools_depth
from outputs import COVERAGE, field_of
from programs import run_cmuts

# How far a stored coverage may stand from the integer depth. Each read's
# contribution is a posterior of one, accumulated in doubles and stored as a
# 32-bit float, so the gap is storage rounding and not model spread.
RELATIVE = 1e-6
ABSOLUTE = 1e-4


def expected_row(depths: dict, columns: int) -> np.ndarray:
    """Builds the coverage row a reference should hold: its depth at every
    covered position, and the fill of zero everywhere else."""
    row = np.zeros(columns)

    for position, depth in depths.items():
        row[position] = depth

    return row


def test_coverage_at_band_zero_matches_samtools_depth_without_the_tags(
        data, falsifiable, tmp_path):
    pinned = strip_tags(drop_adjacent_indels(data, tmp_path), tmp_path)
    depths = samtools_depth(pinned, tmp_path)

    falsifiable(any(depths.values()))

    output = run_cmuts(pinned, tmp_path / "out.h5", band=0, **UNFILTERED)
    coverage = field_of(output, COVERAGE)

    for name, row in rows_by_name(pinned.fasta).items():
        expected = expected_row(depths[name], coverage.shape[1])

        assert np.allclose(coverage[row], expected,
                           rtol=RELATIVE, atol=ABSOLUTE), name


def test_coverage_at_band_zero_stays_within_the_depth_with_deletions(
        data, falsifiable, tmp_path):
    pinned = drop_adjacent_indels(data, tmp_path)
    depths = samtools_depth(pinned, tmp_path, deletions=True)

    falsifiable(any(depths.values()))

    output = run_cmuts(pinned, tmp_path / "out.h5", band=0, **UNFILTERED)
    coverage = field_of(output, COVERAGE)

    for name, row in rows_by_name(pinned.fasta).items():
        most = expected_row(depths[name], coverage.shape[1])

        assert np.all(coverage[row] <= most * (1 + RELATIVE) + ABSOLUTE), name
