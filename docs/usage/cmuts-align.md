
## Purpose

`cmuts align` is a simple wrapper script for performing trimming, demultiplexing, and alignment of one or more FASTQ files against a reference FASTA file, and sorting of the resulting SAM files. This is required if running `cmuts` directly from raw sequencing data.

## Usage

### Alignment

By defaut, `cmuts align` will perform alignment only, which also includes building the `bowtie2` index files and sorting the output SAM files. The syntax is

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*
```

The sorted alignments will be in BAM format and can be found in the directory `$ALIGNMENTS`, with the same basename as the original FASTQ files.

### Paired-End Alignment

To specify mates for paired-end sequencing, pass the additional FASTQ files to the `--pairs` argument *after* the forward reads.

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq* \
  --pairs "$PAIRS"/*.fastq*
```

!!! warning
    The `--pairs` inputs **must** come after the forward reads, otherwise `cmuts align` will not be able to distinguish what you intend to be the mate.

### Demultiplexing

Pre-alignment demultiplexing (via `ultraplex`) can be achieved by passing a CSV of barcodes to the `--barcodes` argument.

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --barcodes "$BARCODES" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*
```

The sorted alignments will again be in BAM format under the directory `$ALIGNMENTS`. However, their names will be given by the respective barcode which they correspond to, not the name of the original FASTQ file.

!!! warning
    If the barcodes are named in the CSV file, then the outputs will have those as names rather than the respective barcodes (see `ultraplex` documentation).

### Trimming

Pre-alignment trimming (via `cutadapt`) can be achieved by passing a sequence to either the `--trim-5` or `--trim-3` argument, depending on which side the adapter belongs to. Separate adapters can be passed for each.

```bash
cmuts align \
  --fasta "$FASTA" \
  --threads "$THREADS" \
  --trim-5 "$TRIM5" \
  --trim-3 "$TRIM3" \
  --output "$ALIGNMENTS" \
  "$FASTQ"/*.fastq*
```

Trimming does not change the format of the output.

!!! warning
    If at least one of `--trim-5` and `--trim-3` are passed with `--barcodes` specified, the trimming will occur first.
