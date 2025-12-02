
The `cmuts` pipeline consists of three main programs:

  1. **`cmuts align`** - Alignment and preprocessing wrapper (uses bowtie2, samtools, fastp)
  2. **`cmuts core`** - Mutation counting engine (use `--pairwise` flag for 2D MaP-seq data)
  3. **`cmuts normalize`** - Normalization and reactivity profile generation (handles both 1D and 2D data automatically)

!!! warning
    Some sequencing vendors (e.g. Ultima) will provide pre-aligned BAM/CRAM files. In such a case, simply ignore the `cmuts align` step.

## Optional Tools

- **`cmuts plot`** - Plot reactivity profiles from HDF5 files
- **`cmuts visualize`** - Overlay reactivities on 3D structures (requires ChimeraX)
- **`cmuts generate`** - Generate synthetic test data for validation
- **`cmuts test`** - Run the integration test suite
