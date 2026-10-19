"""The run-length tags of a sequencer that measures run lengths, which an
Ultima read carries.

The quality of a base is the chance that the run it belongs to is one base
longer or shorter than the read shows, and tp says which. t0 gives the chance
that a run of one base of another letter is missing beside the base, where the
flow order of the read group admits that letter. A base the read does not show
is a base of the cDNA all the same: it pairs with the reference as the base of
its run or of its flow, and it is counted. Each test writes a read with a tag
and an alignment together, and scores it with the tag at several chances, or
with and without the tag.
"""

import numpy as np
import pytest

from alignments import NATIVE, build_placements, strip_tags
from outputs import (COVERAGE, DELETION_RATE, INSERTION_RATE, MISMATCH_RATE,
                     PRIMARY_COUNTED, field_of, fields_of, outputs_agree)
from programs import run_cmuts

# The bases on each side of the run, as in the marginalization tests. A slide
# of the read by one base costs two mismatches.
LEFT = "GC"
RIGHT = "TG"

# The band of the run, and the depth below which a position is not reported.
BAND = 2
EVERY_POSITION = dict(min_depth=0)

# The tp values: the template run is one base longer, or one base shorter,
# than the read run.
LONGER = 1
SHORTER = -1
NONE = 0

# The offset that the t0 string and the quality string add to a score.
PHRED_OFFSET = 33

# The largest score the strings can hold. The chance of an error at this score
# is below one in a billion.
IMPOSSIBLE = ord("~") - PHRED_OFFSET

# Scores in order of rising chance. The chance at a score of zero is one,
# which the decoder caps so that the called length keeps a chance.
SCORES = (IMPOSSIBLE, 30, 20, 10, 0)

# The largest difference allowed between the gap count of an untagged read and
# that of a tagged read whose tags give no chance of an error. The rates are
# stored in single precision.
UNCHANGED = 1e-4

# The largest share of the untagged gap count allowed at the largest chance.
TAKEN = 0.5

# The largest difference allowed between the gap count of an untagged read and
# that of a read whose tag names a template length the read cannot show.
KEPT = 1e-2

# The flow order of the read group, which is the cycle of four bases the
# sequencer repeats. The flows between a run of A and a run of C are T and G.
# A missing run of G can therefore lie between those runs. A missing run of C
# cannot lie between a run of A and a run of G.
FLOW_ORDER = "TGCA"
READ_GROUP = "g"
READ_GROUP_TAG = f"RG:Z:{READ_GROUP}"


def read_group_line(flow_order=FLOW_ORDER) -> str:
    return f"@RG\tID:{READ_GROUP}\tFO:{flow_order}"


READ_GROUP_LINE = read_group_line()

# The flow order as each base of the cycle starts it.
ROTATIONS = [FLOW_ORDER[i:] + FLOW_ORDER[:i] for i in range(len(FLOW_ORDER))]

REVERSE_FLAG = 16

# The quality string of a read that stores no qualities.
NO_QUALITIES = "*"


def tp_tag(values) -> str:
    return "tp:B:c," + ",".join(str(value) for value in values)


def scores_string(scores) -> str:
    return "".join(chr(score + PHRED_OFFSET) for score in scores)


def t0_tag(scores) -> str:
    return "t0:Z:" + scores_string(scores)


def counts_of(output) -> dict:
    """Returns the insertion count and the deletion count, summed over the
    reference.

    A rate divides the events by the coverage of the position where the
    decision is made. That position is the position itself for an insertion
    and one base 3' of it for a deletion, so each rate is multiplied by the
    coverage there.
    """
    fields = fields_of(output, [COVERAGE, INSERTION_RATE, DELETION_RATE])
    coverage = fields[COVERAGE]

    return {
        INSERTION_RATE: float(np.nansum(fields[INSERTION_RATE] * coverage)),
        DELETION_RATE: float(np.nansum(fields[DELETION_RATE][:, :-1] * coverage[:, 1:])),
    }


def coverage_of(output, positions) -> float:
    """Returns the coverage summed over the given positions of the reference."""
    return float(np.nansum(field_of(output, COVERAGE)[:, list(positions)]))


def mismatches_at(output, position) -> float:
    """Returns the mismatch count at one position of the reference: the rate
    times the coverage there."""
    fields = fields_of(output, [COVERAGE, MISMATCH_RATE])

    return float(np.nansum(fields[MISMATCH_RATE][:, position] * fields[COVERAGE][:, position]))


def write_and_run(tmp_path, name, reference, read, cigar, **placement):
    """Scores one read written one way, and returns the output path. A read
    written without tags is written without a quality string, so that it is
    never misread, as a tagged read is not, and the two differ in the length
    errors alone."""
    if "tags" not in placement:
        placement.setdefault("qualities", NO_QUALITIES)

    data = build_placements(tmp_path, name, reference, read, [cigar], **placement)
    output = tmp_path / f"{name}.h5"

    return run_cmuts(data, output, band=BAND, **EVERY_POSITION)


def score(tmp_path, name, reference, read, cigar, **placement) -> dict:
    """Scores one read written one way, and returns its gap counts."""
    return counts_of(write_and_run(tmp_path, name, reference, read, cigar, **placement))


def run_length_tags(read, at, direction, quality):
    """Returns the tags of a read with one tp value, on base at, at the given
    quality. Every other base has no chance of a length error."""
    values = [NONE] * len(read)
    scores = [IMPOSSIBLE] * len(read)
    values[at] = direction
    scores[at] = quality

    return dict(
        qualities=scores_string(scores),
        tags=(tp_tag(values), t0_tag([IMPOSSIBLE] * len(read)), READ_GROUP_TAG),
        header=(READ_GROUP_LINE,),
    )


def missed_run_tags(read, before, quality, header=(READ_GROUP_LINE,)):
    """Returns the tags of a read with the given t0 score on base before and on
    the base after it. Every other base has no chance of a missing run."""
    scores = [IMPOSSIBLE] * len(read)
    scores[before] = quality
    scores[before + 1] = quality

    return dict(
        tags=(tp_tag([NONE] * len(read)), t0_tag(scores), READ_GROUP_TAG),
        header=header,
    )


def assert_falls_with_the_chance(counted, gap):
    """Asserts that the gap count at the impossible score equals the untagged
    count, that it falls at each larger chance, and that it is below the
    allowed share at the largest chance."""
    untagged, *tagged = counted

    assert abs(tagged[0][gap] - untagged[gap]) < UNCHANGED

    for lower, higher in zip(tagged, tagged[1:]):
        assert higher[gap] < lower[gap]

    assert tagged[-1][gap] < TAKEN * untagged[gap]


def assert_rises_with_the_chance(measured):
    """Asserts that a count at the impossible score equals the untagged count,
    that it rises at each larger chance, and that it exceeds the untagged count
    by at least the allowed share at the largest chance."""
    untagged, *tagged = measured

    assert abs(tagged[0] - untagged) < UNCHANGED

    for lower, higher in zip(tagged, tagged[1:]):
        assert higher > lower

    assert tagged[-1] > untagged + TAKEN


# ---------------------------------------------------------------------------
# Run-length errors
# ---------------------------------------------------------------------------

RUN = 5

# A read whose run is one base longer than the reference run, aligned with an
# insertion, and one whose run is one base shorter, aligned with a deletion.
# Each comes with the tp value of its template and the gap it is scored on.
LONGER_READ = (LEFT + "A" * (RUN + 1) + RIGHT, f"{len(LEFT) + RUN}M1I{len(RIGHT)}M",
               SHORTER, INSERTION_RATE)
SHORTER_READ = (LEFT + "A" * (RUN - 1) + RIGHT, f"{len(LEFT) + RUN - 1}M1D{len(RIGHT)}M",
                LONGER, DELETION_RATE)

REFERENCE = LEFT + "A" * RUN + RIGHT

# The first base of the run.
IN_RUN = len(LEFT)


@pytest.mark.parametrize("read, cigar, direction, gap", [LONGER_READ, SHORTER_READ],
                         ids=["longer", "shorter"])
def test_a_run_length_error_takes_the_gap_from_the_reverse_transcriptase(
        tmp_path, read, cigar, direction, gap):
    counted = [score(tmp_path, "untagged", REFERENCE, read, cigar)]

    for quality in SCORES:
        counted.append(score(tmp_path, f"q{quality}", REFERENCE, read, cigar,
                             **run_length_tags(read, IN_RUN, direction, quality)))

    assert_falls_with_the_chance(counted, gap)


def test_the_chances_of_the_bases_of_a_run_add(tmp_path):
    """Tags one base of the run at a score, and then two bases at the same
    score."""
    read, cigar, direction, gap = LONGER_READ
    quality = 20

    one = run_length_tags(read, IN_RUN, direction, quality)
    two = run_length_tags(read, IN_RUN, direction, quality)
    values = [NONE] * len(read)
    scores = [IMPOSSIBLE] * len(read)
    values[IN_RUN] = values[IN_RUN + 1] = direction
    scores[IN_RUN] = scores[IN_RUN + 1] = quality
    two["qualities"] = scores_string(scores)
    two["tags"] = (tp_tag(values), *two["tags"][1:])

    assert score(tmp_path, "two", REFERENCE, read, cigar, **two)[gap] \
        < score(tmp_path, "one", REFERENCE, read, cigar, **one)[gap]


def test_a_tag_pointing_away_from_the_template_leaves_the_gap(tmp_path):
    """Writes the tp value of a longer template on a read whose run is already
    one base longer than the reference."""
    read, cigar, _, gap = LONGER_READ
    untagged = score(tmp_path, "untagged", REFERENCE, read, cigar)
    tagged = score(tmp_path, "away", REFERENCE, read, cigar,
                   **run_length_tags(read, IN_RUN, LONGER, 0))

    assert abs(tagged[gap] - untagged[gap]) < KEPT


# ---------------------------------------------------------------------------
# Missed runs of one base
# ---------------------------------------------------------------------------

# A reference with a run of one base between two runs, and a read without that
# base. The flow order has G between A and C.
MISSED_REFERENCE = LEFT + "AAA" + "G" + "CCC" + RIGHT
MISSED_READ = LEFT + "AAA" + "CCC" + RIGHT
MISSED_CIGAR = f"{len(LEFT) + 3}M1D{3 + len(RIGHT)}M"

# The last base before the gap, and the position of the missing base.
BEFORE_GAP = len(LEFT) + 2
GAP = len(LEFT) + 3

# A reference whose base between the runs is a C, where the flow order admits
# only a T between a run of A and a run of G.
OTHER_BASE_REFERENCE = LEFT + "AAA" + "C" + "GGG" + RIGHT
OTHER_BASE_READ = LEFT + "AAA" + "GGG" + RIGHT

# A reference whose runs are adjacent in the flow order, since T follows A, so
# that no base can be missed between them.
NOT_ADMITTED_REFERENCE = LEFT + "AAA" + "G" + "TTT" + RIGHT
NOT_ADMITTED_READ = LEFT + "AAA" + "TTT" + RIGHT

# A reference whose run between the two read runs holds two bases, so that a
# missed base pairs one of them and the reverse transcriptase deleted the other.
LONGER_GAP_REFERENCE = LEFT + "AAA" + "GG" + "CCC" + RIGHT
LONGER_GAP_CIGAR = f"{len(LEFT) + 3}M2D{3 + len(RIGHT)}M"
LONGER_GAP = (len(LEFT) + 3, len(LEFT) + 4)


def test_a_missed_run_takes_the_deletion_from_the_reverse_transcriptase(tmp_path):
    counted = [score(tmp_path, "untagged", MISSED_REFERENCE, MISSED_READ, MISSED_CIGAR)]

    for quality in SCORES:
        counted.append(score(tmp_path, f"q{quality}", MISSED_REFERENCE, MISSED_READ,
                             MISSED_CIGAR,
                             **missed_run_tags(MISSED_READ, BEFORE_GAP, quality)))

    assert_falls_with_the_chance(counted, DELETION_RATE)


def test_a_missed_run_on_the_reverse_strand_is_read_through_the_complement(tmp_path):
    """Scores the missed run on a reverse-strand read, which was sequenced as
    the reverse complement of what is stored."""
    untagged = score(tmp_path, "untagged", MISSED_REFERENCE, MISSED_READ, MISSED_CIGAR,
                     flag=REVERSE_FLAG)
    tagged = score(tmp_path, "reversed", MISSED_REFERENCE, MISSED_READ, MISSED_CIGAR,
                   flag=REVERSE_FLAG,
                   **missed_run_tags(MISSED_READ, BEFORE_GAP, 0))

    assert tagged[DELETION_RATE] < TAKEN * untagged[DELETION_RATE]


@pytest.mark.parametrize("reference, read, header", [
    (NOT_ADMITTED_REFERENCE, NOT_ADMITTED_READ, (READ_GROUP_LINE,)),
    (MISSED_REFERENCE, MISSED_READ, ()),
], ids=["adjacent-flows", "no-read-group"])
def test_a_missed_run_the_flow_order_does_not_admit_is_not_read(
        tmp_path, reference, read, header):
    untagged = score(tmp_path, "untagged", reference, read, MISSED_CIGAR)
    tagged = score(tmp_path, "tagged", reference, read, MISSED_CIGAR,
                   **missed_run_tags(read, BEFORE_GAP, 0, header=header))

    assert abs(tagged[DELETION_RATE] - untagged[DELETION_RATE]) < UNCHANGED


def test_a_missed_base_unlike_the_reference_is_read_as_a_mismatch(tmp_path):
    """Scores a read whose missing base the flow order gives as a T, against a
    reference holding a C there."""
    untagged = write_and_run(tmp_path, "untagged", OTHER_BASE_REFERENCE,
                             OTHER_BASE_READ, MISSED_CIGAR)
    tagged = write_and_run(tmp_path, "tagged", OTHER_BASE_REFERENCE, OTHER_BASE_READ,
                           MISSED_CIGAR, **missed_run_tags(OTHER_BASE_READ, BEFORE_GAP, 0))

    assert counts_of(tagged)[DELETION_RATE] < TAKEN * counts_of(untagged)[DELETION_RATE]
    assert mismatches_at(tagged, GAP) > mismatches_at(untagged, GAP) + TAKEN


def test_a_missed_base_pairs_inside_a_longer_reference_run(tmp_path):
    """Scores a read lacking a run of two reference bases, of which the
    sequencer can have missed one."""
    measured = [coverage_of(write_and_run(tmp_path, "untagged", LONGER_GAP_REFERENCE,
                                          MISSED_READ, LONGER_GAP_CIGAR), LONGER_GAP)]

    for quality in SCORES:
        measured.append(coverage_of(
            write_and_run(tmp_path, f"q{quality}", LONGER_GAP_REFERENCE, MISSED_READ,
                          LONGER_GAP_CIGAR,
                          **missed_run_tags(MISSED_READ, BEFORE_GAP, quality)),
            LONGER_GAP))

    assert_rises_with_the_chance(measured)


def test_the_flow_order_is_read_as_a_cycle(tmp_path):
    """Scores the missed run under each rotation of the flow order."""
    outputs = [
        write_and_run(tmp_path, rotation, MISSED_REFERENCE, MISSED_READ, MISSED_CIGAR,
                      **missed_run_tags(MISSED_READ, BEFORE_GAP, 0,
                                        header=(read_group_line(rotation),)))
        for rotation in ROTATIONS
    ]

    assert all(outputs_agree(outputs[0], other) for other in outputs[1:])


# ---------------------------------------------------------------------------
# The bases the read does not show
# ---------------------------------------------------------------------------

# The positions of the reference run the shorter read lacks one base of.
RUN_POSITIONS = range(IN_RUN, IN_RUN + RUN)

# A reference with a C after the run, and a read whose run is the length of
# the reference run but which lacks the C, so that an extra base of the run
# pairs the C.
UNLIKE_REFERENCE = LEFT + "A" * RUN + "C" + RIGHT
UNLIKE_READ = LEFT + "A" * RUN + RIGHT
UNLIKE_CIGAR = f"{len(LEFT) + RUN}M1D{len(RIGHT)}M"
UNLIKE = len(LEFT) + RUN


def shorter_run_case(tmp_path, name, quality=None):
    read, cigar, direction, _ = SHORTER_READ
    tags = run_length_tags(read, IN_RUN, direction, quality) if quality is not None else {}

    return write_and_run(tmp_path, name, REFERENCE, read, cigar, **tags)


def missed_run_case(tmp_path, name, quality=None):
    tags = missed_run_tags(MISSED_READ, BEFORE_GAP, quality) if quality is not None else {}

    return write_and_run(tmp_path, name, MISSED_REFERENCE, MISSED_READ, MISSED_CIGAR,
                         **tags)


@pytest.mark.parametrize("case, positions", [
    (shorter_run_case, RUN_POSITIONS),
    (missed_run_case, (GAP,)),
], ids=["shorter-run", "missed-run"])
def test_a_base_the_read_does_not_show_covers_the_base_it_pairs(
        tmp_path, case, positions):
    measured = [coverage_of(case(tmp_path, "untagged"), positions)]

    for quality in SCORES:
        measured.append(coverage_of(case(tmp_path, f"q{quality}", quality), positions))

    assert_rises_with_the_chance(measured)


def test_a_base_the_read_does_not_show_pays_for_a_mismatch(tmp_path):
    """Tags the run of the read as one base longer, where the base after the
    run in the reference is not the base of the run."""
    measured = [mismatches_at(write_and_run(tmp_path, "untagged", UNLIKE_REFERENCE,
                                            UNLIKE_READ, UNLIKE_CIGAR), UNLIKE)]

    for quality in SCORES:
        measured.append(mismatches_at(
            write_and_run(tmp_path, f"q{quality}", UNLIKE_REFERENCE, UNLIKE_READ,
                          UNLIKE_CIGAR, **run_length_tags(UNLIKE_READ, IN_RUN, LONGER,
                                                          quality)),
            UNLIKE))

    assert_rises_with_the_chance(measured)


# ---------------------------------------------------------------------------
# Tags the model cannot read
# ---------------------------------------------------------------------------


def test_a_read_without_tp_is_scored_as_untagged(tmp_path):
    """Writes t0 and a read group on a read, and no tp."""
    read, cigar, _, _ = LONGER_READ
    untagged = build_placements(tmp_path, "untagged", REFERENCE, read, [cigar])
    partial = build_placements(tmp_path, "partial", REFERENCE, read, [cigar],
                               tags=(t0_tag([0] * len(read)), READ_GROUP_TAG),
                               header=(READ_GROUP_LINE,))

    first = run_cmuts(untagged, tmp_path / "untagged.h5", band=BAND, **EVERY_POSITION)
    second = run_cmuts(partial, tmp_path / "partial.h5", band=BAND, **EVERY_POSITION)

    assert outputs_agree(first, second)


def test_tags_that_admit_no_error_score_as_no_tags(tmp_path):
    """Writes tp values of no error on every base, and no t0, and compares with
    the same read without tags and without a quality string, which is never
    misread either."""
    read, cigar, _, _ = LONGER_READ
    first = write_and_run(tmp_path, "untagged", REFERENCE, read, cigar)
    second = write_and_run(tmp_path, "tagged", REFERENCE, read, cigar,
                           tags=(tp_tag([NONE] * len(read)), READ_GROUP_TAG),
                           header=(READ_GROUP_LINE,))

    assert outputs_agree(first, second)


# ---------------------------------------------------------------------------
# Generated reads
# ---------------------------------------------------------------------------

# The dataset with a large chance of a missing base, whose reads differ from
# untagged reads in the tags only.
MISSING = "missing"


def test_stripping_the_tags_keeps_the_reads_and_changes_the_gaps(catalogue, tmp_path):
    """Scores a generated dataset with its tags, and again with the tags
    stripped."""
    data = catalogue(MISSING, NATIVE)
    stripped = strip_tags(data, tmp_path)

    with_tags = run_cmuts(data, tmp_path / "tagged.h5", **EVERY_POSITION)
    without = run_cmuts(stripped, tmp_path / "stripped.h5", **EVERY_POSITION)

    assert np.array_equal(field_of(with_tags, PRIMARY_COUNTED),
                          field_of(without, PRIMARY_COUNTED))

    gaps_with = counts_of(with_tags)
    gaps_without = counts_of(without)

    assert gaps_with[DELETION_RATE] < gaps_without[DELETION_RATE]
    assert sum(gaps_with.values()) < sum(gaps_without.values())
