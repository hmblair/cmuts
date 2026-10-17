"""Aligning reads with cmuts align.

The reads are taken from an alignment that cmuts gen wrote and read back as
FASTQ, so every test has a reference that its reads align to.
"""

from __future__ import annotations

import random
import shutil

import pytest

from alignments import (NATIVE, REVERSE_BIT, SECONDARY_BIT, SUPPLEMENTARY_BIT, UNMAPPED_BIT,
                        measure_dataset, reverse_run_lengths)
from outputs import outputs_agree
from oracle import FIRST_TAG_COLUMN, sequences
from programs import (CMUTS_ALIGN, attempt, run_align, run_cmuts, samtools, samtools_into,
                      try_align, try_cmuts)

# A BAM file stores the length of a read name in one byte.
NAME_LIMIT = 254

# The dataset the reads are taken from. One is enough, because no contract here
# depends on how the reads are laid out.
DATASET = "plain"

# The preset every test aligns under, other than those naming their own.
PRESET = "sr"

# How much of a read each mate holds when one read is written out as a pair.
# Anything above half leaves the two mates overlapping, which the merger
# requires.
MATE_SHARE = 0.7

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None,
                                reason="minimap2 is required to align reads")


# ---------------------------------------------------------------------------
# Reading the output
# ---------------------------------------------------------------------------


def header_of(bam) -> str:
    return samtools("view", "-H", bam)


def reference_lines(bam) -> list:
    return [line for line in header_of(bam).splitlines() if line.startswith("@SQ")]


def read_names(bam) -> list:
    return [line.split("\t")[0] for line in samtools("view", bam).splitlines()]


def checksums_in(text: str) -> set:
    """The M5 field of every @SQ line, whether of a header or of a dictionary."""
    return {field for line in text.splitlines() if line.startswith("@SQ")
            for field in line.split("\t") if field.startswith("M5:")}


# ---------------------------------------------------------------------------
# Writing the reads
# ---------------------------------------------------------------------------


def rewrite_names(reads, destination, rename):
    """Copies a FASTQ file with every read name replaced."""
    lines = reads.read_text().splitlines()

    with open(destination, "w") as handle:
        for at in range(0, len(lines), 4):
            handle.write("@" + rename(lines[at][1:]) + "\n")
            handle.write("\n".join(lines[at + 1:at + 4]) + "\n")

    return destination


def unaligned_bam(reads, destination):
    """Writes the reads out as an unaligned BAM file."""
    samtools("import", "-0", reads, "-o", destination)

    return destination


def split_into_mates(reads, directory):
    """Writes a FASTQ file out as two, each mate holding most of one read, so
    that the two mates overlap."""
    complement = str.maketrans("ACGTN", "TGCAN")
    lines = reads.read_text().splitlines()
    first = directory / "r1.fastq"
    second = directory / "r2.fastq"

    with open(first, "w") as one, open(second, "w") as two:
        for at in range(0, len(lines), 4):
            name, bases, _, qualities = lines[at:at + 4]
            take = int(len(bases) * MATE_SHARE)

            one.write(f"{name}\n{bases[:take]}\n+\n{qualities[:take]}\n")
            two.write(f"{name}\n{bases[-take:][::-1].translate(complement)}\n+\n"
                      f"{qualities[-take:][::-1]}\n")

    return first, second


@pytest.fixture(scope="session")
def alignable(catalogue, tmp_path_factory):
    """A FASTA file, and a FASTQ file holding reads that align to it."""
    data = catalogue(DATASET, NATIVE)
    reads = tmp_path_factory.mktemp("reads") / "reads.fastq"

    samtools_into(reads, "fastq", data.bam)

    return data.fasta, reads


@pytest.fixture
def aligned(alignable, tmp_path):
    """The reads aligned to the reference, as one sorted BAM file."""
    fasta, reads = alignable

    return run_align(fasta, tmp_path / "aligned.bam", PRESET, reads)


# ---------------------------------------------------------------------------
# What the output holds
# ---------------------------------------------------------------------------


def test_the_output_is_coordinate_sorted(aligned):
    assert "SO:coordinate" in header_of(aligned)


def test_the_output_holds_the_reads_that_aligned(aligned):
    assert len(read_names(aligned)) > 0


def test_every_reference_declares_the_checksum_of_its_sequence(alignable, aligned):
    fasta, _ = alignable
    declared = checksums_in(header_of(aligned))

    assert len(declared) == len(reference_lines(aligned))
    assert declared <= checksums_in(samtools("dict", fasta))


def test_a_read_name_too_long_for_a_bam_file_is_shortened(alignable, tmp_path):
    fasta, reads = alignable
    renamed = rewrite_names(reads, tmp_path / "long.fastq",
                            lambda name: name + "x" * NAME_LIMIT)

    aligned = run_align(fasta, tmp_path / "long.bam", PRESET, renamed)

    assert read_names(aligned)
    assert max(len(name) for name in read_names(aligned)) == NAME_LIMIT


def test_the_output_can_be_counted(alignable, aligned, tmp_path):
    fasta, _ = alignable
    counted = measure_dataset("aligned", (aligned,), fasta)

    assert try_cmuts(counted, tmp_path / "counts.h5").returncode == 0


# ---------------------------------------------------------------------------
# The forms the reads arrive in
# ---------------------------------------------------------------------------


def test_the_same_reads_align_alike_from_a_bam_and_from_a_fastq(
        alignable, aligned, tmp_path):
    fasta, reads = alignable
    unaligned = unaligned_bam(reads, tmp_path / "unaligned.bam")

    from_bam = run_align(fasta, tmp_path / "frombam.bam", PRESET, unaligned)

    assert samtools("view", from_bam) == samtools("view", aligned)


@pytest.mark.skipif(shutil.which("vsearch") is None,
                    reason="vsearch is required to merge a pair of mates")
def test_a_pair_of_mates_is_merged_and_aligned(alignable, tmp_path):
    fasta, reads = alignable
    first, second = split_into_mates(reads, tmp_path)

    merged = run_align(fasta, tmp_path / "merged.bam", PRESET, first, second)

    assert len(read_names(merged)) > 0


# ---------------------------------------------------------------------------
# What is refused
# ---------------------------------------------------------------------------


def test_an_already_aligned_bam_is_refused(catalogue, tmp_path):
    data = catalogue(DATASET, NATIVE)

    assert try_align(data.fasta, tmp_path / "out.bam", PRESET, data.bam).returncode != 0


def test_a_bam_given_as_one_file_of_a_pair_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    unaligned = unaligned_bam(reads, tmp_path / "unaligned.bam")

    assert try_align(fasta, tmp_path / "out.bam", PRESET,
                     reads, unaligned).returncode != 0


def test_a_pair_of_mates_is_refused_under_a_long_read_preset(alignable, tmp_path):
    fasta, reads = alignable
    first, second = split_into_mates(reads, tmp_path)

    assert try_align(fasta, tmp_path / "out.bam", "map-ont",
                     first, second).returncode != 0


def test_a_preset_that_aligns_reads_against_each_other_is_refused(alignable, tmp_path):
    fasta, reads = alignable

    assert try_align(fasta, tmp_path / "out.bam", "ava-ont", reads).returncode != 0


def test_a_truncated_bam_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    whole = unaligned_bam(reads, tmp_path / "unaligned.bam").read_bytes()
    truncated = tmp_path / "truncated.bam"
    output = tmp_path / "out.bam"

    truncated.write_bytes(whole[:len(whole) // 2])

    assert try_align(fasta, output, PRESET, truncated).returncode != 0
    assert not output.exists()


def test_a_file_holding_neither_reads_nor_alignments_is_refused(alignable, tmp_path):
    fasta, _ = alignable

    assert try_align(fasta, tmp_path / "out.bam", PRESET, fasta).returncode != 0


def test_a_minimap2_index_in_place_of_the_fasta_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    index = tmp_path / "ref.mmi"

    attempt(["minimap2", "-d", index, fasta])

    assert try_align(index, tmp_path / "out.bam", PRESET, reads).returncode != 0


def test_a_fasta_that_repeats_a_reference_name_is_refused(alignable, tmp_path):
    """The report names the file and the reference that repeats, rather than
    the failure that the repeated name causes further down the pipeline."""
    fasta, reads = alignable
    repeated = tmp_path / "repeated.fasta"
    name = next(iter(sequences(fasta)))

    repeated.write_text(fasta.read_text() * 2)
    result = try_align(repeated, tmp_path / "out.bam", PRESET, reads)
    reason = result.stderr.strip().splitlines()[-1]

    assert result.returncode != 0
    assert str(repeated) in reason
    assert f'"{name}"' in reason


def test_an_existing_output_is_kept_unless_overwriting_is_asked_for(alignable, aligned):
    fasta, reads = alignable
    before = aligned.read_bytes()

    assert try_align(fasta, aligned, PRESET, reads).returncode != 0
    assert aligned.read_bytes() == before

    run_align(fasta, aligned, PRESET, reads, overwrite=True)


def test_a_directory_in_place_of_the_output_is_refused(alignable, tmp_path):
    fasta, reads = alignable

    assert try_align(fasta, tmp_path, PRESET, reads).returncode != 0


def test_a_run_that_is_refused_leaves_no_partial_file(alignable, tmp_path):
    fasta, _ = alignable
    output = tmp_path / "out.bam"

    assert try_align(fasta, output, PRESET, fasta).returncode != 0
    assert list(tmp_path.glob("out.bam*")) == []


def test_the_program_reports_its_version():
    assert attempt([*CMUTS_ALIGN, "--version"]).returncode == 0


# ---------------------------------------------------------------------------
# The run-length tags of an Ultima read
# ---------------------------------------------------------------------------

# The tagged dataset, whose reads carry tp, t0 and a read group with a flow order.
TAGGED = "tagged"
CARRIED = ("tp:B:c,", "t0:Z:", "RG:Z:")
ALIGNMENT_TYPE = "tp:A:"
NOT_PRIMARY_FLAGS = SECONDARY_BIT | SUPPLEMENTARY_BIT

# Prefixes that make a read name and a reference name look like tags, which only the
# optional fields of an alignment are.
NAME_LIKE_A_TAG = "tp:A:"
REFERENCE_LIKE_A_TAG = "t0:Z:"

# Tagged reads added to the dataset's own: reads of random bases, which align nowhere, and
# reads joined from the ends of two references, which minimap2 splits into a primary and
# a hard-clipped supplementary alignment.
UNALIGNABLE = "unalignable"
JOINED = "joined"
ADDED_READS = 3
ADDED_LENGTH = 150
ARM = 100
SEED = 7


def carried_fields(fields) -> dict:
    """The tp, t0 and RG fields of a SAM line, by their prefix."""
    return {prefix: field for field in fields[FIRST_TAG_COLUMN:] for prefix in CARRIED
            if field.startswith(prefix)}


def sequenced_reads(data, directory):
    """Writes the reads of an alignment out as the sequencer gave them, as an unaligned
    BAM file: every reverse-strand read is turned back, its tp and t0 values with it,
    and the read group is kept. Every other read is given a name that looks like a tag.
    Returns the file, and the bases and tags of each read in the order they were
    sequenced."""
    complement = str.maketrans("ACGTN", "TGCAN")
    head = [line for line in header_of(data.bam).splitlines() if line.startswith("@RG")]
    lines = ["@HD\tVN:1.6\tSO:unknown", *head]
    sequenced = {}

    for at, line in enumerate(samtools("view", "-F", str(NOT_PRIMARY_FLAGS),
                                       data.bam).splitlines()):
        fields = line.split("\t")
        name, flag, bases, qualities = fields[0], int(fields[1]), fields[9], fields[10]
        tags = carried_fields(fields)

        if flag & REVERSE_BIT:
            bases = bases[::-1].translate(complement)
            qualities = qualities[::-1]
            tags = {prefix: reverse_run_lengths(field) for prefix, field in tags.items()}

        name = NAME_LIKE_A_TAG + name if at % 2 else name
        sequenced[name] = (bases, tags)
        lines.append("\t".join([name, str(UNMAPPED_BIT), "*", "0", "0", "*", "*", "0",
                                "0", bases, qualities, *tags.values()]))

    group = head[0].split("\t")[1][len("ID:"):]

    for name, bases in added_reads(data.fasta):
        tags = invented_tags(len(bases), group, name)
        sequenced[name] = (bases, tags)
        lines.append("\t".join([name, str(UNMAPPED_BIT), "*", "0", "0", "*", "*", "0",
                                "0", bases, "I" * len(bases), *tags.values()]))

    sam = directory / "unaligned.sam"
    sam.write_text("\n".join(lines) + "\n")
    output = directory / "unaligned.bam"
    samtools("view", "-b", "-o", output, sam)

    return output, sequenced


def added_reads(fasta):
    """Returns the reads that align nowhere and the reads joined from two references, as
    names and bases, the same on every call."""
    rng = random.Random(SEED)
    long = [bases for bases in sequences(fasta).values() if len(bases) >= 2 * ARM]

    for i in range(ADDED_READS):
        yield f"{UNALIGNABLE}{i}", "".join(rng.choice("ACGT") for _ in range(ADDED_LENGTH))
        yield f"{JOINED}{i}", long[2 * i][:ARM] + long[2 * i + 1][-ARM:]


def invented_tags(length: int, group: str, name: str) -> dict:
    """Returns tp, t0 and RG fields for a read of a length, the same for every call with
    the same name."""
    rng = random.Random(name)
    tp = ",".join(str(rng.choice((-1, 0, 1))) for _ in range(length))
    t0 = "".join(chr(ord("!") + rng.randrange(41)) for _ in range(length))

    return {"tp:B:c,": "tp:B:c," + tp, "t0:Z:": "t0:Z:" + t0, "RG:Z:": "RG:Z:" + group}


def without_run_lengths_on_supplementary(bam, directory):
    """Copies an alignment with tp and t0 removed from its supplementary alignments."""
    lines = []

    for line in samtools("view", "-h", bam).splitlines():
        fields = line.split("\t")

        if not line.startswith("@") and int(fields[1]) & SUPPLEMENTARY_BIT:
            fields = [*fields[:FIRST_TAG_COLUMN],
                      *(field for field in fields[FIRST_TAG_COLUMN:]
                        if not field.startswith(("tp:B:", "t0:Z:")))]

        lines.append("\t".join(fields))

    sam = directory / "bare-supplementary.sam"
    sam.write_text("\n".join(lines) + "\n")
    output = directory / "bare-supplementary.bam"
    samtools("view", "-b", "-o", output, sam)

    return output


def renamed_references(data, directory):
    """Writes the reference out with every name made to look like a tag."""
    lines = data.fasta.read_text().splitlines()
    renamed = directory / "renamed.fasta"
    names = (">" + REFERENCE_LIKE_A_TAG + line[1:] if line.startswith(">") else line
             for line in lines)
    renamed.write_text("".join(line + "\n" for line in names))

    return renamed


def test_the_run_length_tags_of_an_unaligned_bam_are_kept(catalogue, tmp_path):
    """Aligns the reads of the tagged dataset from an unaligned BAM, with tagged reads that
    align nowhere and tagged reads joined from two references. Every read comes out once,
    as a primary alignment or unmapped. A primary alignment carries the tp and t0 values
    the read was given, in the order of its stored sequence on either strand, and an
    unmapped record in the order they were sequenced; both carry the read group, and the
    header keeps its flow order. A hard-clipped supplementary alignment carries the
    values of the whole read in the order of its strand, which do not cover its bases, so
    cmuts hmm counts it as it counts one without them. Read and reference names that look
    like tags pass through unchanged."""
    data = catalogue(TAGGED, NATIVE)
    reads, sequenced = sequenced_reads(data, tmp_path)
    fasta = renamed_references(data, tmp_path)
    aligned = run_align(fasta, tmp_path / "aligned.bam", PRESET, reads)

    groups = [line for line in header_of(aligned).splitlines() if line.startswith("@RG")]
    assert groups and all("FO:" in line for line in groups)

    reported = []
    strands = set()
    unmapped = supplementary = 0

    for line in samtools("view", aligned).splitlines():
        fields = line.split("\t")
        flag = int(fields[1])
        reverse = bool(flag & REVERSE_BIT)
        bases, expected = sequenced[fields[0]]
        tags = carried_fields(fields)
        oriented = {prefix: reverse_run_lengths(field) if reverse else field
                    for prefix, field in expected.items()}

        assert not any(field.startswith(ALIGNMENT_TYPE)
                       for field in fields[FIRST_TAG_COLUMN:])

        if flag & SUPPLEMENTARY_BIT:
            supplementary += "H" in fields[5]
            assert len(tags["tp:B:c,"].split(",")) - 1 == len(bases) > len(fields[9])
            assert tags == oriented, fields[0]
            continue

        reported.append(fields[0])

        if flag & UNMAPPED_BIT:
            unmapped += bool(tags)
            assert tags == expected, fields[0]
            continue

        strands.add(reverse)

        assert fields[2].startswith(REFERENCE_LIKE_A_TAG)
        assert len(tags["tp:B:c,"].split(",")) - 1 == len(fields[9])
        assert tags == oriented, fields[0]

    assert sorted(reported) == sorted(sequenced)
    assert strands == {False, True}
    assert unmapped > 0 and supplementary > 0

    both = dict(alignment_type="primary,supplementary", min_depth=0)
    counted = measure_dataset(TAGGED, (aligned,), fasta)
    stripped = without_run_lengths_on_supplementary(aligned, tmp_path)
    bare = measure_dataset(TAGGED, (stripped,), fasta)

    assert outputs_agree(run_cmuts(counted, tmp_path / "carried.h5", **both),
                         run_cmuts(bare, tmp_path / "bare.h5", **both))


def test_a_fastq_still_aligns_without_its_comments(alignable, tmp_path):
    """Aligns reads whose FASTQ comments are not SAM tags, which are not copied."""
    fasta, reads = alignable
    lines = reads.read_text().splitlines()
    commented = tmp_path / "commented.fastq"
    commented.write_text("\n".join(f"{line} zzcomment" if at % 4 == 0 else line
                                   for at, line in enumerate(lines)) + "\n")

    output = run_align(fasta, tmp_path / "aligned.bam", PRESET, commented)

    assert len(read_names(output)) > 0
    assert "zzcomment" not in samtools("view", output)
