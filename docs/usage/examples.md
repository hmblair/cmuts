## Overview

### Environment Setup

These variables will be used throughout the examples; feel free to set them to a value appropriate for your computer. The latter three specify the location of various intermediate files and outputs from the pipeline.

```bash
THREADS=8
ALIGNMENTS="examples/alignments"
COUNTS="examples/counts.h5"
PROFILES="examples/profiles.h5"
```

All files used in these examples can be found in the repo, under `examples`.

### Pipeline Steps

The `cmuts` pipeline consists of three parts:

  1. The `cmuts align` wrapper,
  2. The `cmuts core` program, and
  3. Either `cmuts normalize` or `cmuts cov`, depending on whether 1D or 2D MaP-seq is being performed.

# MaP-seq

```bash
PARENT="examples/map-seq"
```

A directory of FASTQ files and a reference FASTA file is required.

```bash
FASTQ="$PARENT/fastq"
FASTA="$PARENT/ref.fasta"
```

In addition, `cmuts normalize` will need to know which datasets correspond to which condition. How this is specified depends on whether demultiplexing is performed or not.

## Without Demultiplexing

The following arrays contain the basenames of the FASTQ files corresponding to each condition. A parallel array containing the name of each condition is necessary as well.

```bash
MODS=("bicine-2A3", "bicine-DMS")
NOMODS=("bicine-DMSO", "bicine-ETH")
NAMES=("2A3", "DMS")
```

`NOMODS` may be empty, in the case where background subtraction is not to be performed.


One can run `cmuts align` directly on the FASTQ files and then pass the output to `cmuts core` to locate modifications.

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*

cmuts core \
  -f "$FASTA" \
  -o "$COUNTS" \
  --filter-coverage \
  --no-insertions \
  "$ALIGNMENTS"/*
```

This will output a single HDF5 file, `$COUNTS`, which will contain one dataset for each input. One can then loop over all input files, and pass the respective dataset to `cmuts normalize`.

```bash
for ((IX=0; IX<${#MODS[@]}; IX++)); do
  MOD_DS="$ALIGNMENTS"/${MODS[IX]}
  NOMOD_DS="$ALIGNMENTS"/${NOMODS[IX]}
  NAME=${NAMES[IX]}
  cmuts normalize \
    -o "$PROFILES" \
    --mod "$MOD_DS" \
    --nomod "$NOMOD_DS" \
    --group "$NAME" \
    "$COUNTS"
done
```

## With Demultiplexing

```bash
PARENT="examples/map-seq-dmux"
```

The main difference when working with FASTQ files which must be demultiplexed first is that `MODS` and `NOMODS` should refer to the **barcodes** for each condition, rather than the files themselves.

```bash
MODS=("GCCTGGGTGGCT", "TGACCATGTATA")
NOMODS=("AAGGACCACTGG", "CTTATTACACAC")
NAMES=("2A3", "DMS")
```

We use a different directory of FASTQ files for this example. Naturally, in addition to the previous inputs, a comma-separated file containing the barcodes to use for demultiplexing must be provided to `cmuts align`.

```bash
FASTQ="$PARENT/fastq"
FASTA="$PARENT/ref.fasta"
BARCODES="$PARENT/barcodes.fasta"
```

Otherwise, the pipeline remains the same:

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --barcodes "$BARCODES" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*

cmuts core \
  -f "$FASTA" \
  -o "$COUNTS" \
  --filter-coverage \
  --no-insertions \
  "$ALIGNMENTS"/*

for ((IX=0; IX<${#MODS[@]}; IX++)); do
  MOD_DS="$ALIGNMENTS"/${MODS[IX]}
  NOMOD_DS="$ALIGNMENTS"/${NOMODS[IX]}
  NAME=${NAMES[IX]}
  cmuts normalize \
    -o "$PROFILES" \
    --mod "$MOD_DS" \
    --nomod "$NOMOD_DS" \
    --group "$NAME" \
    "$COUNTS"
done
```

# Two-Dimensional MaP-seq

## Residue-Residue Correlations

```bash
PARENT="examples/mohca"
```

This example uses **MOHCA-seq** data.

```bash
FASTQ="$PARENT/fastq"
FASTA="$PARENT/ref.fasta"
```

The first two steps are the same, save for the `--joint` flag passed to `cmuts core`.

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*

cmuts core \
  --joint \
  -f "$FASTA" \
  -o "$COUNTS" \
  --filter-coverage \
  --no-insertions \
  "$ALIGNMENTS"/*
```

The third step uses `cmuts cov` to postprocess the 2D data.

```bash
for ((IX=0; IX<${#MODS[@]}; IX++)); do
  MOD_DS="$ALIGNMENTS"/${MODS[IX]}
  cmuts cov \
    -o "$PROFILES" \
    --dataset "$MOD_DS" \
    --mutual-information \
    "$COUNTS"
done
```

## Mutate-And-Map

```bash
PARENT="examples/m2"
```

```bash
FASTQ="$PARENT/fastq"
FASTA="$PARENT/ref.fasta"
```

