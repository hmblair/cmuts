
The functionality of `cmuts` is distributed among several subprograms. All accept the `--log` argument to redirect their outputs to a log file.

```bash
cmuts [--log LOGFILE] {align|core|normalize|visualize|test}
```

The main `cmuts` pipeline comprises the three subprograms

  1. [`cmuts align`](/cmuts/usage/cmuts-align),
  3. [`cmuts core`](/cmuts/usage/cmuts-core), and
  2. [`cmuts normalize`](/cmuts/usage/cmuts-normalize).

!!! warning
    Some sequencing vendors (e.g. Ultima) will provide pre-aligned BAM/CRAM files. In such a case, simply ignore the `cmuts align` step.

See [Examples](/cmuts/usage/examples) for detailed usage instructures.

[`cmuts visualize`](/cmuts/usage/cmuts-visualize) may be used for overlaying generated reactivities on tertiary structures.
