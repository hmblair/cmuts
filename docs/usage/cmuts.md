# cmuts

`cmuts` is a program for counting mutations and computing reactivity profiles in MaP-seq experiments. It features fast, compiled C++ code with native multithreading support, streamed IO, and direct output to compressed HDF5 files.

## Basic Usage

### Modification Counting

To count modifications in a single aligned HTS file:

```bash
cmuts -o out.h5 -f seq.fasta sorted.bam
```

For the MPI-enabled version, use multiple threads:

```bash
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted.bam
```

Multiple SAM/BAM/CRAM files can be processed simultaneously:

```bash
mpirun -np 8 cmuts -o out.h5 -f seq.fasta sorted1.bam sorted2.bam ...
```

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

## Special Use Cases

### Tokenization Only

Generate tokenized sequences without processing alignments:

```bash
cmuts --tokenize -f seq.fasta -o out.h5
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
