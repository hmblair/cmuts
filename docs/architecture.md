# Architecture Overview

This document describes the high-level architecture of cmuts for developers and maintainers.

## Data Flow Pipeline

```
FASTQ → Alignment → BAM/CRAM → cmuts-core → HDF5 → cmuts-normalize → Reactivity
```

### Stage 1: Mutation Counting (C++)

The C++ core (`cmuts-core`) processes aligned reads to count mutations:

1. **Input Parsing**: Reads BAM/CRAM files using htslib with streaming I/O
2. **Reference Loading**: Loads FASTA reference sequences into memory
3. **Mutation Detection**: For each aligned read:
   - Parses CIGAR string to identify mismatches, insertions, deletions
   - Records mutation type and position in count matrices
   - Handles ambiguous deletions using mutation-informed heuristics
4. **Output**: Writes HDF5 files containing count tensors

### Stage 2: Normalization (Python)

The Python package (`cmuts-normalize`) transforms counts to reactivity:

1. **Count Aggregation**: Combines counts across replicates/conditions
2. **Coverage Filtering**: Masks positions below coverage threshold
3. **Reactivity Calculation**: `reactivity = modifications / coverage`
4. **Normalization**: Applies scheme (raw, percentile, or outlier-based)
5. **Pairwise Analysis**: Computes correlation matrices and mutual information

## Module Responsibilities

### C++ Core (`src/cpp/`)

| File | Responsibility |
|------|----------------|
| `main.cpp` | Entry point, argument parsing |
| `bam.cpp` | BAM file reading and iteration |
| `cram.cpp` | CRAM file reading with reference handling |
| `cmuts.cpp` | Core mutation counting logic |
| `hdf5.cpp` | HDF5 output with chunking and compression |
| `mpi.cpp` | MPI coordination for distributed processing |
| `fasta.cpp` | FASTA reference sequence parsing |
| `utils.cpp` | Utility functions |

### Python Package (`src/python/cmuts/`)

| Module | Responsibility |
|--------|----------------|
| `internal.py` | Core data structures (`ProbingData`, `Opts`) and normalization |
| `normalize/` | CLI for reactivity normalization |
| `visualize/` | Plotting and visualization tools |

## Data Structures

### HDF5 Count File Structure

```
/
├── sequences          # Tokenized reference sequences (optional)
├── mod/               # Modified condition
│   ├── counts-1d      # Shape: (refs, seqlen, 4, 7) - 1D mutation counts
│   └── counts-2d      # Shape: (refs, seqlen, seqlen, 2, 2) - pairwise counts
└── nomod/             # Unmodified control (optional)
    ├── counts-1d
    └── counts-2d
```

### 1D Count Tensor Axes

- Axis 0: Reference sequences
- Axis 1: Sequence positions
- Axis 2: Reference nucleotide (A, C, G, T)
- Axis 3: Observation type (A, C, G, T, del, ins, term)

### 2D Count Tensor Axes

- Axes 0-2: Same as 1D (ref, pos_i, pos_j)
- Axes 3-4: Binary mutation state at each position (0=match, 1=mutation)

## MPI Coordination Model

When built with MPI support, cmuts distributes work across ranks:

1. **Rank 0 (Coordinator)**:
   - Opens input BAM/CRAM file
   - Distributes reference sequence assignments to workers
   - Aggregates final results and writes output

2. **Worker Ranks**:
   - Receive assigned reference sequences
   - Process reads mapping to assigned references
   - Send local counts back to coordinator

3. **Load Balancing**:
   - References assigned by estimated read count
   - Larger references may be split across multiple ranks

## Extension Points

### Adding New Mutation Types

1. Extend the observation type enum in `src/cpp/cmuts.hpp`
2. Update CIGAR parsing in `cmuts.cpp` to detect new type
3. Modify HDF5 schema and Python loader

### Adding New Normalization Schemes

1. Add enum value to `NormScheme` in `internal.py`
2. Implement `_get_norm_<scheme>()` function
3. Update `_get_norm()` dispatch logic

### Adding New Output Formats

1. Create new writer class in `src/cpp/`
2. Add CLI flag in argument parser
3. Wire up in `main.cpp`

## Performance Considerations

- **Memory**: Count tensors can be large; use chunked HDF5 I/O
- **I/O Bound**: BAM parsing is typically the bottleneck
- **Thread Safety**: OpenMP parallelization requires thread-local accumulators
- **MPI Overhead**: Small references may not benefit from distribution

## Dependencies

### C++ (System)

- htslib: BAM/CRAM parsing
- HDF5: Output format
- zlib: Compression
- OpenMP: Threading
- MPI (optional): Distributed computing

### C++ (Fetched)

- xtensor: N-dimensional arrays
- xtl: Template utilities
- argparse: CLI parsing

### Python

- numpy, scipy: Numerical computing
- dask: Lazy array operations
- h5py: HDF5 file access
- matplotlib: Visualization
- biopython: Sequence handling
