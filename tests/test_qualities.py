"""The mapping-quality histogram, against what samtools makes of the same file.

There is one bin per mapping quality, from 0 to 254, and the same bins in
every row. A quality of 255 means "unavailable" and is rejected by the filter,
so it has no bin. Only the reads that passed the filter are counted, so
samtools is given the same criteria.
"""

import numpy as np
import pytest

from oracle import rows_by_name, samtools_quality_histogram
from outputs import PRIMARY_COUNTED, QUALITIES, QUALITY_BINS, field_of
from programs import run_cmuts

# Unfiltered, and above the mapping quality of most reads the generator
# writes, so the histogram is checked over two different sets of reads.
MAPPING_QUALITIES = (0, 30)


def expected_row(histogram: dict) -> np.ndarray:
    """Builds the row samtools implies. Bin i holds quality i."""
    row = np.zeros(QUALITY_BINS)

    for quality, count in histogram.items():
        row[quality] += count

    return row


def compare_against_samtools(output, data, min_mapq: int) -> int:
    """Asserts every reference's row against samtools, and returns how many
    references it ran over. The criterion is passed to both, since the two do
    not share one."""
    expected = samtools_quality_histogram(data, min_mapq=min_mapq)
    row_of = rows_by_name(data.fasta)

    written = field_of(output, QUALITIES)

    assert written.shape[1] == QUALITY_BINS

    for name, histogram in expected.items():
        assert np.array_equal(written[row_of[name]], expected_row(histogram)), name

    return len(expected)


@pytest.mark.parametrize("min_mapq", MAPPING_QUALITIES)
def test_histogram_matches_samtools(data, falsifiable, tmp_path, min_mapq):
    output = tmp_path / "out.h5"
    run_cmuts(data, output, min_mapq=min_mapq)

    falsifiable(compare_against_samtools(output, data, min_mapq) > 0)


def test_every_row_sums_to_the_primary_reads_counted(data, falsifiable, tmp_path):
    """Every counted primary read has a quality below 255, so every one lands in
    a bin."""
    output = tmp_path / "out.h5"
    run_cmuts(data, output, min_mapq=0)

    counted = field_of(output, PRIMARY_COUNTED)
    written = field_of(output, QUALITIES).sum(axis=1)

    falsifiable(counted.sum() > 0)

    assert np.array_equal(written, counted)
