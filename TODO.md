# `cmuts run` — unified pipeline command

## Interface

```
# Pre-demultiplexed
cmuts run \
  --fasta ref.fasta \
  --mod 2A3_rep1.fq 2A3_rep2.fq \
  --nomod DMSO_rep1.fq DMSO_rep2.fq \
  -o profiles.h5

# Barcoded (names in CSV must be "mod" and "nomod")
cmuts run \
  --fasta ref.fasta \
  --mod pool1.fq pool2.fq \
  --nomod pool1.fq pool2.fq \
  --barcodes barcodes.csv \
  -o profiles.h5
```

## Arguments

```
Required:
  --fasta FILE              Reference FASTA
  --mod FILE [FILE ...]     Modified sample FASTQ(s)
  -o/--output FILE          Output HDF5 file

Optional:
  --nomod FILE [FILE ...]   Control FASTQ(s)
  --group NAME              Probe name in output HDF5 (default: "")
  --barcodes FILE           Barcode CSV (names must be "mod" and "nomod")
  --threads N               Parallelism (default: 1)
  --overwrite               Overwrite existing output

Alignment:
  --trim-5 SEQ              5' adapter to trim
  --trim-3 SEQ              3' adapter to trim
  --local                   Bowtie2 local alignment mode
  --pairs FILE [FILE ...]   Paired-end mates

Counting:
  --no-insertions           Exclude insertions
  --no-deletions            Exclude deletions
  --no-mismatches           Exclude mismatches
  --min-mapq N              Min mapping quality (default: 10)
  --min-phred N             Min base quality (default: 10)
  --min-length N            Min read length (default: 2)
  --max-length N            Max read length (default: 1024)
  --no-reverse              Ignore reverse-complement reads
  --only-reverse            Only use reverse-complement reads

Normalization:
  --norm METHOD             ubr|outlier|raw (default: ubr)
  --clip-low                Clamp negatives to 0
  --clip-high               Cap at 1.0
  --blank-5p N              NaN out N 5' terminal bases
  --blank-3p N              NaN out N 3' terminal bases
  --blank-cutoff N          Min reads per position (default: 10)
  --norm-cutoff N           Min reads for normalization (default: 500)
```

## Design decisions

### Pure function, no side effects
The current pipeline (`cmuts align` → `cmuts core` → `cmuts normalize`) leaves intermediate artifacts (index directory, BAM files, raw counts HDF5) on disk. For `cmuts run`, all intermediates go into a temp dir and are cleaned up on success. The command is a pure function: FASTA + FASTQs in, HDF5 + figures out. This avoids confusing wet-lab users with files they don't need and prevents stale intermediates from biting them.

### No config files or passthrough args
Early designs considered config files (TOML/CSV) for specifying experimental layouts and `--core-args`/`--align-args` passthrough flags for advanced options. Both were rejected. Config files add complexity for the primary audience (wet-lab users who just want reactivities). Passthrough args are awkward and hard to document. Instead, the most commonly needed knobs from each stage are promoted to top-level flags on `cmuts run`. The full subcommands remain available for power users who need the long tail of options.

### One group per invocation
Early designs tried to handle all probe types (2A3, DMS, etc.) in a single invocation, which required either a config file or a complex multi-group CLI syntax. Neither was simpler than just running the command twice. One group per invocation keeps the interface flat and predictable. The second invocation appends to the same output HDF5 via `--group`.

### Consistent interface for demultiplexing
The demultiplexed and pre-demultiplexed cases use the same interface. `--mod` and `--nomod` always take FASTQ files. When `--barcodes` is present, the FASTQs are demultiplexed first, and the barcode CSV names (`mod`/`nomod`) determine which demuxed samples are which. If `--mod` and `--nomod` point to the same pooled FASTQs, the backend deduplicates so demuxing only happens once. This avoids having two different modes with different flag semantics — the user doesn't need to think about whether their data is pre-demuxed or not beyond adding `--barcodes`.

### Flag subset matches cmuts-space
The exposed flags are the same subset that `cmuts-space` (the Hugging Face Spaces web UI) exposes. This is a battle-tested set of options that covers typical MaP-seq experiments without overwhelming users. Options like `--pairwise`, `--tokenize`, ambiguous deletion handling, and per-event quality filter overrides are intentionally omitted — users who need those know enough to use the individual subcommands.

## Behavior

1. Parse args, resolve all FASTQ paths.
2. If `--barcodes`: deduplicate `--mod`/`--nomod` FASTQs, demux via ultraplex. Barcode CSV names `mod`/`nomod` determine which demuxed samples are which.
3. Align all FASTQs to reference (bowtie2) in a temp dir.
4. Run `cmuts core` over all BAMs → raw counts HDF5 in temp dir.
5. Run `cmuts normalize` with the mod/nomod BAM groups → write to output HDF5 under `--group`.
6. Generate figures alongside the output HDF5.
7. Clean up temp dir.

All intermediates (index, BAMs, raw counts) live in a temp dir and are cleaned up on success. Output is just the HDF5 and figures.

## Multi-probe usage

Run once per probe, appending to the same output HDF5:

```bash
cmuts run --fasta ref.fasta --mod 2A3/*.fq --nomod DMSO/*.fq -o profiles.h5 --group 2A3
cmuts run --fasta ref.fasta --mod DMS/*.fq --nomod ETH/*.fq  -o profiles.h5 --group DMS
```

## Implementation

- New Python script at `src/python/cmuts-run`.
- Shells out to `cmuts align`, `cmuts core`, `cmuts normalize` (stays in sync with existing tools).
- Add `run` subcommand to `src/scripts/cmuts` dispatcher.
