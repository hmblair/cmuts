"""The run-length tags of a sequencer that measures run lengths, which an
Ultima read carries.

The quality of a base is the chance that the run it belongs to is one base
longer or shorter than the read shows, and tp says which. t0 gives the chance
that a run of one base of another letter is missing beside the base, where the
flow order of the read group admits that letter. Each test writes a read with
a tag and an alignment together, and scores it with the tag at several
chances, or with and without the tag.
"""

import gzip
import re
import struct
import zlib
from dataclasses import replace

import numpy as np
import pytest

from alignments import (NATIVE, add_alignment_type, build_placements, strip_tags,
                        tags_in_sequenced_order)
from datasets import DATASETS
from outputs import (COVERAGE, DELETION_RATE, INSERTION_RATE, PRIMARY_COUNTED, field_of,
                     fields_of, outputs_agree)
from programs import run_cmuts, samtools, try_cmuts

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
READ_GROUP_LINE = f"@RG\tID:{READ_GROUP}\tFO:{FLOW_ORDER}"
READ_GROUP_TAG = f"RG:Z:{READ_GROUP}"

REVERSE_FLAG = 16


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


def score(tmp_path, name, reference, read, cigar, **placement) -> dict:
    """Scores one read written one way, and returns its gap counts."""
    data = build_placements(tmp_path, name, reference, read, [cigar], **placement)
    output = tmp_path / f"{name}.h5"

    run_cmuts(data, output, band=BAND, **EVERY_POSITION)

    return counts_of(output)


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
# base. The flow order has G between A and C. It does not have C between A
# and G.
MISSED_REFERENCE = LEFT + "AAA" + "G" + "CCC" + RIGHT
MISSED_READ = LEFT + "AAA" + "CCC" + RIGHT
NOT_ADMITTED_REFERENCE = LEFT + "AAA" + "C" + "GGG" + RIGHT
NOT_ADMITTED_READ = LEFT + "AAA" + "GGG" + RIGHT
MISSED_CIGAR = f"{len(LEFT) + 3}M1D{3 + len(RIGHT)}M"

# The last base before the gap.
BEFORE_GAP = len(LEFT) + 2


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
], ids=["letter-off-the-flow", "no-read-group"])
def test_a_missed_run_the_flow_order_does_not_admit_is_not_read(
        tmp_path, reference, read, header):
    untagged = score(tmp_path, "untagged", reference, read, MISSED_CIGAR)
    tagged = score(tmp_path, "tagged", reference, read, MISSED_CIGAR,
                   **missed_run_tags(read, BEFORE_GAP, 0, header=header))

    assert abs(tagged[DELETION_RATE] - untagged[DELETION_RATE]) < UNCHANGED


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


# ---------------------------------------------------------------------------
# Generated reads
# ---------------------------------------------------------------------------

# The dataset with a large chance of a missing base, whose reads differ from
# untagged reads in the tags only.
MISSING = "missing"
TAGGED = "tagged"
PLAIN = "plain"


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


def test_the_substitution_score_of_tagged_reads_reaches_only_them(catalogue, tmp_path):
    """Scores a tagged dataset and an untagged one at two substitution
    scores."""
    for name, changes in ((TAGGED, True), (PLAIN, False)):
        data = catalogue(name, NATIVE)
        low = run_cmuts(data, tmp_path / f"{name}-low.h5", ultima_phred=10)
        high = run_cmuts(data, tmp_path / f"{name}-high.h5", ultima_phred=40)

        assert outputs_agree(low, high) != changes, name


# ---------------------------------------------------------------------------
# The order of the tags
# ---------------------------------------------------------------------------

# How a run reads the tags of a file whose reverse-strand reads hold them in the
# order the bases were sequenced, as a copy from the unaligned reads leaves them.
SEQUENCED = "sequenced"
STORED = "stored"

UNTAGGED_WARNING = "carry no usable tp tag"
RUN_LENGTHS = ("tp", "t0")
NO_FLOW_ORDER = "none"
ORDER_WARNING = "--tp-output-order"


def test_tags_in_sequenced_order_are_read_back_in_stored_order(data, falsifiable, tmp_path):
    """Scores a dataset, and a copy whose reverse-strand reads hold their tags in the
    order the bases were sequenced, read with and without --ultima-tag-order. The copy
    read without the option differs wherever reads are tagged, half of them being on the
    reverse strand."""
    flipped = tags_in_sequenced_order(data, tmp_path)
    tagged = DATASETS[data.name].get("tagged", 0) > 0

    stored = run_cmuts(data, tmp_path / "stored.h5", **EVERY_POSITION)
    read_back = run_cmuts(flipped, tmp_path / "read-back.h5", ultima_tag_order=SEQUENCED,
                          **EVERY_POSITION)
    misread = run_cmuts(flipped, tmp_path / "misread.h5", **EVERY_POSITION)

    falsifiable(tagged)

    assert outputs_agree(stored, read_back)
    assert outputs_agree(stored, misread) != tagged


def warnings_of(data, output, **options) -> str:
    """Counts a dataset, which must succeed, and returns what it printed to stderr."""
    result = try_cmuts(data, output, **options)

    assert result.returncode == 0, result.stderr

    return result.stderr


def test_reads_that_lost_their_run_lengths_are_reported(data, falsifiable, tmp_path):
    """Strips tp and t0 from every read and keeps the read group. The run warns exactly
    where reads that named a read group with a flow order were tagged."""
    stripped = strip_tags(data, tmp_path, tags=RUN_LENGTHS)
    parameters = DATASETS[data.name]
    lost = parameters.get("tagged", 0) > 0 and parameters.get("flow_order") != NO_FLOW_ORDER

    falsifiable(data.mapped > 0)

    assert (UNTAGGED_WARNING in warnings_of(stripped, tmp_path / "stripped.h5")) == lost


def test_reads_that_name_no_read_group_are_not_reported(catalogue, tmp_path):
    """Strips the run-length tags of a tagged dataset with and without its read group. A
    read that names no read group is not taken for an Ultima read, as an Illumina read
    counted beside Ultima reads is not."""
    tagged = catalogue(TAGGED, NATIVE)
    (tmp_path / "kept").mkdir()
    (tmp_path / "removed").mkdir()
    kept = strip_tags(tagged, tmp_path / "kept", tags=RUN_LENGTHS)
    removed = strip_tags(tagged, tmp_path / "removed")

    assert UNTAGGED_WARNING in warnings_of(kept, tmp_path / "kept.h5")
    assert UNTAGGED_WARNING not in warnings_of(removed, tmp_path / "removed.h5")


# A read group of another sequencer, which declares no flow order.
OTHER_GROUP = "other"
OTHER_GROUP_LINE = f"@RG\tID:{OTHER_GROUP}"
UNTAGGED_COUNT = re.compile(r"(\d+) of (\d+) primary reads")


def untagged_reads(data, directory, every_other=None):
    """Writes the reads of a dataset without tp and t0, under a header that also declares
    a read group with no flow order. every_other, 0 or 1, keeps every second read only, and
    with no value every read is kept and every second one is moved to the other group."""
    header, records = [], []

    for line in samtools("view", "-h", data.bam).splitlines():
        (header if line.startswith("@") else records).append(line)

    lines = []

    for at, line in enumerate(records):
        if every_other is not None and at % 2 != every_other:
            continue

        fields = [field for field in line.split("\t")
                  if not field.startswith(("tp:B:", "t0:Z:"))]

        if (every_other is None and at % 2) or every_other == 1:
            fields = [f"RG:Z:{OTHER_GROUP}" if field.startswith("RG:Z:") else field
                      for field in fields]

        lines.append("\t".join(fields))

    stem = "mixed" if every_other is None else f"part{every_other}"
    sam = directory / f"{stem}.sam"
    sam.write_text("\n".join([*header, OTHER_GROUP_LINE, *lines]) + "\n")
    output = directory / f"{stem}.bam"
    samtools("view", "-b", "-o", output, sam)

    return output


def untagged_count(data, output):
    """Returns the untagged and total reads that a run warns about, or None."""
    found = UNTAGGED_COUNT.search(warnings_of(data, output))

    return tuple(map(int, found.groups())) if found else None


def test_only_reads_of_a_read_group_with_a_flow_order_are_counted(catalogue, tmp_path):
    """Splits the reads of a tagged dataset, stripped of their run-length tags, between its
    read group and one with no flow order: as two files in either order, and within one
    file. Every run warns about the reads of the first group alone."""
    data = catalogue(TAGGED, NATIVE)
    ultima = untagged_reads(data, tmp_path, every_other=0)
    other = untagged_reads(data, tmp_path, every_other=1)
    mixed = untagged_reads(data, tmp_path)

    alone = untagged_count(replace(data, bams=(ultima,)), tmp_path / "alone.h5")
    pairs = {"forwards": (ultima, other), "backwards": (other, ultima), "mixed": (mixed,)}

    assert alone and alone[0] == alone[1] > 0
    assert untagged_count(replace(data, bams=(other,)), tmp_path / "other.h5") is None

    for name, bams in pairs.items():
        counted = untagged_count(replace(data, bams=bams), tmp_path / f"{name}.h5")

        assert counted == alone, name


def test_a_declared_tag_order_that_disagrees_is_reported(data, falsifiable, tmp_path):
    """Adds the @PG line of Ultima's aligner, which declares that the tags follow the stored
    sequence, and counts the file under each order. The same words outside the command
    line, or the option with no value, declare nothing."""
    for name in ("declared", "described", "empty"):
        (tmp_path / name).mkdir()

    declared = tags_in_sequenced_order(
        data, tmp_path / "declared",
        "@PG\tID:UA\tPN:UA\tCL:ua --tp-output-order=reference --align=true\n")
    described = tags_in_sequenced_order(
        data, tmp_path / "described",
        "@PG\tID:UA\tPN:UA\tDS:--tp-output-order=reference\tCL:ua --align=true\n")
    empty = tags_in_sequenced_order(
        data, tmp_path / "empty",
        "@PG\tID:UA\tPN:UA\tCL:ua --tp-output-order= --align=true\n")

    falsifiable(True)

    assert ORDER_WARNING in warnings_of(declared, tmp_path / "s.h5",
                                        ultima_tag_order=SEQUENCED)
    assert ORDER_WARNING not in warnings_of(declared, tmp_path / "t.h5",
                                            ultima_tag_order=STORED)
    assert ORDER_WARNING not in warnings_of(described, tmp_path / "d.h5",
                                            ultima_tag_order=SEQUENCED)

    for order in (SEQUENCED, STORED):
        assert ORDER_WARNING not in warnings_of(empty, tmp_path / f"e-{order}.h5",
                                                ultima_tag_order=order)


def test_a_tp_tag_of_another_type_before_the_run_lengths_is_passed_over(data, falsifiable,
                                                                         tmp_path):
    """Puts minimap2's tp:A alignment type ahead of the tags of every read, as minimap2 -y
    leaves a copied run-length tp:B:c, and scores the file against the original."""
    typed = add_alignment_type(data, tmp_path)

    original = run_cmuts(data, tmp_path / "original.h5", **EVERY_POSITION)
    passed_over = run_cmuts(typed, tmp_path / "typed.h5", **EVERY_POSITION)

    falsifiable(DATASETS[data.name].get("tagged", 0) > 0)

    assert outputs_agree(original, passed_over)


# A BGZF block holds at most this much data, and the file ends in an empty block.
BGZF_BLOCK = 0xff00
BGZF_END = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

# The tp field as a BAM record stores it: the tag, then the array and element types.
RUN_LENGTH_ARRAY = b"tpBc"

# A count of four-byte elements whose product with four wraps a 32-bit size to four.
OVERLONG_ARRAY = b"tpBi" + struct.pack("<I", 0x40000001)

# How many bytes of the tp field a truncated record keeps: the tag, the types, and two
# of the four bytes of the count.
TRUNCATED_AT = len(RUN_LENGTH_ARRAY) + 2


def write_bgzf(data: bytes, path):
    """Writes bytes out as BGZF, which htslib reads a BAM file as."""
    with open(path, "wb") as handle:
        for at in range(0, len(data), BGZF_BLOCK):
            chunk = data[at:at + BGZF_BLOCK]
            packer = zlib.compressobj(6, zlib.DEFLATED, -15)
            body = packer.compress(chunk) + packer.flush()
            handle.write(b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff"
                         + struct.pack("<H2sHH", 6, b"BC", 2, len(body) + 25) + body
                         + struct.pack("<II", zlib.crc32(chunk), len(chunk)))

        handle.write(BGZF_END)

    return path


def first_tagged_record(raw: bytes) -> tuple:
    """Returns where the first record holding a tp array starts and ends in an
    uncompressed BAM file, and where its tp field starts."""
    text = struct.unpack_from("<i", raw, 4)[0]
    at = 8 + text
    references = struct.unpack_from("<i", raw, at)[0]
    at += 4

    for _ in range(references):
        at += 4 + struct.unpack_from("<i", raw, at)[0] + 4

    while True:
        end = at + 4 + struct.unpack_from("<i", raw, at)[0]
        field = raw.find(RUN_LENGTH_ARRAY, at, end)

        if field >= 0:
            return at, end, field

        at = end


def without_tags_on_first_tagged_read(data, directory):
    """Copies a dataset with every tag removed from the first read that holds a tp."""
    lines, done = [], False

    for line in samtools("view", "-h", data.bam).splitlines():
        if not done and not line.startswith("@") and "\ttp:B:" in line:
            line = "\t".join(line.split("\t")[:11])
            done = True

        lines.append(line)

    sam = directory / "bare.sam"
    sam.write_text("\n".join(lines) + "\n")
    output = directory / "bare.bam"
    samtools("view", "-b", "-o", output, sam)

    return replace(data, bams=(output,))


@pytest.mark.parametrize("damage", ["overlong", "truncated"])
def test_an_array_tag_that_runs_past_its_record_is_passed_over(catalogue, tmp_path, damage):
    """Breaks the tp field of the first tagged read, by a count that runs past the end of
    the record or by cutting the record off inside the count, and counts the file. The
    read is counted as one without tags; a tag after a broken field cannot be found."""
    data = catalogue(TAGGED, NATIVE)
    raw = gzip.decompress(data.bam.read_bytes())
    start, end, field = first_tagged_record(raw)

    if damage == "overlong":
        record = raw[start:field] + OVERLONG_ARRAY + raw[field + len(OVERLONG_ARRAY):end]
    else:
        record = raw[start:field + TRUNCATED_AT]
        record = struct.pack("<i", len(record) - 4) + record[4:]

    broken = replace(data, bams=(tmp_path / "broken.bam",))
    write_bgzf(raw[:start] + record + raw[end:], broken.bams[0])
    bare = without_tags_on_first_tagged_read(data, tmp_path)

    warnings_of(broken, tmp_path / "broken.h5", **EVERY_POSITION)

    assert outputs_agree(tmp_path / "broken.h5",
                         run_cmuts(bare, tmp_path / "bare.h5", **EVERY_POSITION))
