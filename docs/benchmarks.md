## Purpose

The `benchmarks/` directory holds the scripts that compare `cmuts` against the
two other MaP-seq mutation counters, **rf-count** (RNAFramework) and
**shapemapper2**. There are three benchmarks:

- **profile** — wall-clock time and peak memory of each tool across sweeps of
  query count, reference count, and reference length.
- **correctness** — per-nucleotide agreement: with deletion spreading off and
  every tool configured to count the same events, how close is each tool's raw
  mutation rate to cmuts'. This is the Fig 2a control that isolates deletion
  spreading.
- **accuracy** — how well each tool's reactivity predicts base pairing in
  deposited 3D structures (PDB130), with each tool tuned for accuracy.

These are development tools, not part of the installed `cmuts` package, so they
are run as plain scripts rather than `cmuts` subcommands.

## Shared tool layer

rf-count and shapemapper2 read different formats and take different flags from
cmuts. [`benchmarks/external.py`](https://github.com/hmblair/cmuts/blob/master/benchmarks/external.py)
is the single place that knows how to *invoke* and *parse* each of them; all
three benchmarks import it, so the flag spelling and the `.rc` / counts.txt
parsers live in exactly one place.

Quality and processing knobs are carried in `external.Params` and mapped to each
tool's own flags. Every benchmark passes the `Params` matching its intent:

| Benchmark   | Intent                                                   |
|-------------|----------------------------------------------------------|
| profile     | mirror the synthetic-data generation (counts insertions) |
| correctness | match cmuts defaults so the rate is the same quantity    |
| accuracy    | tune each tool for accuracy (its standard mutation set)  |

## Environment

The external tools are resolved from the environment, so a machine without them
still runs the cmuts-only paths (an absent tool is skipped):

| Variable     | Meaning                          | Default          |
|--------------|----------------------------------|------------------|
| `CMUTS`      | cmuts dispatcher                 | `cmuts` on PATH  |
| `RF_COUNT`   | RNAFramework `rf-count`          | `rf-count`       |
| `RF_RCTOOLS` | RNAFramework `rf-rctools`        | `rf-rctools`     |
| `RF_NORM`    | RNAFramework `rf-norm` (accuracy)| `rf-norm`        |
| `SAMTOOLS`   | `samtools`                       | `samtools`       |
| `SM2_DIR`    | shapemapper2 install directory   | unset → skipped  |

`accuracy` additionally requires **ciffy**, **biopython**, and **scipy** (cmuts
already depends on the latter two). BAMs passed to `correctness` and `accuracy`
for the per-reference shapemapper2 loop must be **sorted and indexed**.

## Running

### profile

```bash
python benchmarks/profile.py --references --runs 3 --threads 1 8 \
  --shapemapper-dir /path/to/shapemapper2 --output profile-results.csv
```

Sweeps are selected with `--queries`, `--references`, `--lengths` (at least one
required). Synthetic data is generated with `cmuts generate` and cached under
`./.cases`. rf-count is skipped if not on PATH; shapemapper2 runs only when
`--shapemapper-dir` is given. Results are appended to the CSV; the figure
sources read it.

### correctness

```bash
python benchmarks/correctness.py alignments.bam -f references.fasta \
  -o correctness.csv
```

Writes one row per (reference, position) with the cmuts / rf-count /
shapemapper2 raw mutation rates and the absolute difference of each external
tool from cmuts. As a sanity check on the parsers, rf-count and shapemapper2
should agree closely with each other on deep coverage (Pearson ≳ 0.99).

### accuracy

```bash
python benchmarks/accuracy.py \
  --fasta references.fasta \
  --structures /path/to/cif/dir \
  --dms-mod dms_mod.bam --dms-nomod dms_nomod.bam \
  --2a3-mod 2a3_mod.bam --2a3-nomod 2a3_nomod.bam \
  --per-ref per_reference.tsv --summary summary.tsv
```

Each condition takes a modified and an unmodified (nomod) sample, and **each tool
uses its own native normalization**:

- **cmuts** — `cmuts core` + `cmuts normalize --mod --nomod` for all three spread
  modes (`nospread`, `default`, `uniform`).
- **rf-count** — `rf-norm` (treated-only Zubradt scoring, 2-8% normalization) on
  the modified sample, counted without deletions (rf-count cannot spread them, so
  excluding them is its accuracy-optimal config; matches the reference pipeline).
- **shapemapper2** — its `make_reactivity_profiles` + `normalize_profiles` step on
  the modified + untreated samples (run with shapemapper's bundled Python; the
  normalization factor is computed library-wide across references).

Pass `--profiles profiles.h5` in place of the BAMs to load precomputed cmuts
profiles (e.g. from `generate_profiles.sh`) and skip recomputing them. The
summary reports mean/median AUC, Pearson, and Spearman per (tool/mode ×
condition); `--blank-5p` / `--blank-3p` / `--downsample` tune the cmuts
normalization.

## Sherlock notes

Two issues found while validating these on the cluster:

1. **Do not `source .venv/bin/activate`.** The venv activation drops samtools
   from `PATH` and rf-count then fails with "samtools is not in PATH". Instead
   source the cluster deps and call the interpreter directly:

   ```bash
   source "$CMUTS_HOME/__sherlock_deps"
   export PATH="$CMUTS_HOME/bin:/path/to/RNAFramework:$PATH"
   export SM2_DIR=/path/to/shapemapper2
   "$CMUTS_HOME/.venv/bin/python" benchmarks/correctness.py ...
   ```

2. `correctness`/`accuracy` call `cmuts normalize`, which on older builds
   allocated the full SNR-scaling array and could run out of memory on inputs
   with thousands of references. This is fixed in current builds (normalization
   no longer plots); deploy a current `cmuts` before large runs.
