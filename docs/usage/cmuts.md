`cmuts core` performs the main job of the `cmuts` pipeline, which is determining the location and type of mutations, insertions, and deletions in a collection of aligned reads.

## Inputs and Usage

All modes of `cmuts core` require two inputs to run:

1. **Reference sequences**, stored in a FASTA file,
2. **Aligned reads**, stored in one or more SAM/BAM/CRAM files

A generic call will look like
```bash
cmuts core -o OUTPUT -f FASTA [OPTIONAL ARGUMENTS] FILE
```

!!! warning
    The alignments must be sorted before passing to `cmuts`. If you are generating them via `cmuts align`, then they will automatically be sorted for you.

## Outputs

The output of `cmuts core` depends on which mode it is run in.

### Standard

Most uses of `cmuts core` (without any of the flags specified below) will output an HDF5 file with a dataset of shape \(n \times l \times 4 \times 7\). The former two dimensions specify the reference sequence and residue respectively, and the latter two specify the type of modification seen in accordance with the following array:

![cmuts core heatmap](figures/heatmap.png)


### Joint

The `--joint` flag instructs `cmuts core` to count *pairs* of mutations

### Low-Mem

### Required Inputs

1. **Reference sequences**: A FASTA file specified with the `-f` option
2. **Aligned reads**: One or more SAM/BAM/CRAM files

The program will automatically create necessary index files (.bai or .crai) and a custom binary .binfa file from the FASTA file. If input files are not sorted, `samtools sort` will be called automatically.

### Output Format

The output HDF5 file contains one dataset per input file, named by the file path (without extension). Each dataset has shape `n × l × 4 × 7` where:

- `n` = number of sequences
- `l` = maximum sequence length  
- Dimension 2 = original base (A, C, G, T)
- Dimension 3 = mutated base or inserted base

## Command Line Options

### Core Options

**`-o, --output`** : Output HDF5 filename (required)

**`-f, --fasta`** : Reference FASTA file (required)

### Output Control

**`--overwrite`** : Overwrite existing HDF5 file

**`--compression`** : HDF5 compression level (0-9, default: 3)

### Analysis Modes

**`--tokenize`** : Tokenize reference sequences and store in `sequence` dataset

**`--joint`** : Compute joint distribution of mutations over all position pairs

**`--low-mem`** : Record only modification locations and coverage (not type), reducing memory usage by 2x+

### Quality Filtering

**`--min-mapq`** : Mapping quality threshold (default: 10)

**`--min-phred`** : PHRED score threshold (default: 10)

**`--quality-window`** : Quality check window size around each base (default: 2)

**`--min-length`** : Minimum alignment length (default: 2)

**`--max-length`** : Maximum alignment length (default: 10,000)

### Modification Detection

**`--max-indel-length`** : Maximum indel length to consider (default: 10)

**`--collapse`** : Collapse modifications within this distance in a read (default: 2)

**`--no-mismatches`** : Exclude mismatches from modification counts

**`--no-insertions`** : Exclude insertions from modification counts (recommended for MaP-seq)

**`--no-deletions`** : Exclude deletions from modification counts

### Deletion Handling

**`--uniform-spread`** : Uniformly spread ambiguous deletions

**`--mutation-spread`** : Spread ambiguous deletions according to existing mutation profile

**`--disable-ambiguous`** : Use alignment-provided deletions only

**`--contiguous-ambiguous`** : Allow only contiguous regions as ambiguous deletions

### Sampling and Coverage

**`--subsample`** : Random read acceptance probability (default: 1.0)

**`--filter-coverage`** : Apply modification filters to matches (recommended for MaP-seq)

### Performance

**`--chunk-size`** : Internal buffer size per thread in references (default: 128)

**`--print-every`** : Progress update frequency in reads processed (default: 1,000)


## Example Use Cases

### Tokenization Only

Generate tokenized sequences without processing alignments:

```bash
cmuts --tokenize -f seq.fasta -o counts.h5
```

### Joint Modification Analysis

Compute joint modification distributions:

```bash
cmuts --joint -o joint.h5 -f seq.fasta sorted.bam
```

Output shape: `n × l × l × 2 × 2` where the final dimensions indicate modification presence (1) or absence (0) at each position pair.

### Memory-Efficient Analysis

For large datasets, use low-memory mode:

```bash
cmuts --low-mem -o out.h5 -f seq.fasta sorted.bam
```

## Recommended MaP-seq Settings

For standard MaP-seq analysis, use:

```bash
cmuts --no-insertions --filter-coverage -o out.h5 -f seq.fasta sorted.bam
```

These flags:
- Exclude insertions which are not relevant for chemical modification detection
- Apply consistent quality filters to both modifications and coverage calculations
