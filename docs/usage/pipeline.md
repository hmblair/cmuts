
The `cmuts` pipeline consists of three independent programs:

  1. The `cmuts align` wrapper,
  2. The `cmuts core` program, and
  3. Either `cmuts normalize` or `cmuts pairwise`, depending on whether 1D or 2D MaP-seq is being performed.

!!! warning
    Some sequencing vendors (e.g. Ultima) will provide pre-aligned BAM/CRAM files. In such a case, simply ignore the `cmuts align` step.
