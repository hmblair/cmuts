
The functionality of `cmuts` is distributed among several subprograms. All accept the `--log` argument to redirect their outputs to a log file.

<!-- BEGIN AUTO-GENERATED USAGE -->
```text
Usage: cmuts [--log LOGFILE] <command> [options]

Commands:
  align       Trim, demultiplex, and align FASTQ/BAM reads to a reference FASTA
  core        Count mutations, insertions, and deletions in aligned reads
  generate    Create synthetic test data for validating the pipeline
  normalize   Compute normalized reactivity profiles from core counts
  plot        Plot reactivity profiles from HDF5 files
  visualize   Overlay a reactivity profile onto a 3D RNA structure
  test        Run the integration test suite

Run 'cmuts <command> -h' for command-specific options.
```
<!-- END AUTO-GENERATED USAGE -->

The main `cmuts` pipeline comprises the three subprograms

  1. [`cmuts align`](/cmuts/usage/cmuts-align),
  3. [`cmuts core`](/cmuts/usage/cmuts-core), and
  2. [`cmuts normalize`](/cmuts/usage/cmuts-normalize).

!!! warning
    Some sequencing vendors (e.g. Ultima) will provide pre-aligned BAM/CRAM files. In such a case, simply ignore the `cmuts align` step.

See [Examples](/cmuts/usage/examples) for detailed usage instructions.

[`cmuts visualize`](/cmuts/usage/cmuts-visualize) may be used for overlaying generated reactivities on tertiary structures.

[`cmuts plot`](/cmuts/usage/cmuts-plot) may be used for plotting reactivity profiles.

[`cmuts test`](/cmuts/usage/cmuts-test) runs the integration test suite to verify installation.
