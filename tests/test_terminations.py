"""Every counted read stops on exactly one base, or copies the reference to its end.

The reverse transcriptase copies the reference 3' to 5'. A read that stops at a
base has the base 3' of it as its 5'-most base, and the rate at a base is the
share of the reads pairing that 5'-most base which stopped there. The reads
which stop at a base are therefore the rate there times the coverage one base
3' of it. A read that pairs the first base copied the whole reference and
stopped at no base. The two kinds together are every read counted.

A band of zero pins each read's stop to the base 5' of its placed span. A wider
band may place the stop more than one base beyond the 5' end, where nothing is
accumulated.

The depth floor is lowered so that every position reports a rate. A position
of exactly one read can otherwise fall a rounding step below the floor, which
would drop that read from the sum.
"""

import numpy as np

from alignments import drop_adjacent_indels
from filters import UNFILTERED
from oracle import reference_lengths, rows_by_name
from outputs import PRIMARY_COUNTED, COVERAGE, TERMINATION_RATE, field_of
from programs import run_cmuts

# How far a summed posterior may stand from the whole number of reads. Each
# read contributes a posterior of one, accumulated in doubles and stored as a
# 32-bit float, so the gap is storage rounding.
RELATIVE = 1e-6
ABSOLUTE = 1e-4


def stopped_inside(rate, coverage):
    """Returns the reads of each reference that stopped at a base of it."""
    return np.nansum(rate[:, :-1] * coverage[:, 1:], axis=1)


def copied_to_the_end(coverage):
    """Returns the reads of each reference that pair its first base."""
    return coverage[:, 0]


def test_terminations_at_band_zero_sum_to_the_reads_counted(data, falsifiable, tmp_path):
    pinned = drop_adjacent_indels(data, tmp_path)

    output = run_cmuts(pinned, tmp_path / "out.h5", band=0, min_depth=0, **UNFILTERED)
    rate = np.asarray(field_of(output, TERMINATION_RATE))
    coverage = np.asarray(field_of(output, COVERAGE))
    counted = np.asarray(field_of(output, PRIMARY_COUNTED))

    falsifiable(np.any(counted > 0))

    stopped = stopped_inside(rate, coverage) + copied_to_the_end(coverage)

    assert np.allclose(stopped, counted, rtol=RELATIVE, atol=ABSOLUTE)


def test_no_termination_is_reported_at_the_last_base(data, falsifiable, tmp_path):
    """A termination at the last base would be decided from a base past the
    reference, which no read pairs."""
    output = run_cmuts(data, tmp_path / "out.h5", min_depth=0, **UNFILTERED)
    rate = np.asarray(field_of(output, TERMINATION_RATE))
    lengths = reference_lengths(data.fasta)
    row_of = rows_by_name(data.fasta)

    falsifiable(len(lengths) > 0)

    for name, length in lengths.items():
        assert np.isnan(rate[row_of[name], length - 1])
