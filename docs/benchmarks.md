## Purpose

The `benchmarks/` directory holds the scripts that compare `cmuts` against the
two other MaP-seq mutation counters, **rf-count** (RNAFramework) and
**shapemapper2**:

- **profile** — wall-clock time and peak memory of each tool across sweeps of
  query count, reference count, and reference length. Standalone.
- **profiles** — runs each tool's full pipeline (counting, background
  subtraction, and normalization) once on the modified/unmodified BAMs and writes
  one HDF5 of reactivity profiles: cmuts at each spread mode and tuned to each
  tool, plus rf-count and shapemapper2. This is the shared input for the two
  scoring benchmarks, so the tools are run once rather than per benchmark.
- **correctness** — reads the profiles HDF5; per-nucleotide agreement of cmuts
  (with spreading off, tuned to each tool) against that tool's own profile. This
  is the Fig 2a control that isolates deletion spreading.
- **accuracy** — reads the profiles HDF5; how well each tool's reactivity
  predicts base pairing in deposited 3D structures (PDB130), scored as AUC /
  Pearson / Spearman against a ground truth read from the deposited mmCIF
  hydrogen-bond annotations.

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

| Benchmark | Intent                                                          |
|-----------|-----------------------------------------------------------------|
| profile   | mirror the synthetic-data generation (counts insertions)        |
| profiles  | one shared map-seq counting config across every tool (`MAP`)    |

`profiles.py` drives both scoring benchmarks, so all tools share a single
map-seq config (`MAP` in `profiles.py`): insertions off, deletions on and
right-aligned, collapse 2, mapq/phred 10, eval-surrounding (quality window 1),
duplicates kept. The only thing that varies across the cmuts datasets is the
deletion spread mode.

## Environment

The external tools are resolved from the environment, so a machine without them
still runs the cmuts-only paths (an absent tool is skipped):

| Variable     | Meaning                          | Default          |
|--------------|----------------------------------|------------------|
| `CMUTS`      | cmuts dispatcher                 | `cmuts` on PATH  |
| `RF_COUNT`   | RNAFramework `rf-count`          | `rf-count`       |
| `RF_RCTOOLS` | RNAFramework `rf-rctools`        | `rf-rctools`     |
| `RF_NORM`    | RNAFramework `rf-norm` (profiles)| `rf-norm`        |
| `SAMTOOLS`   | `samtools`                       | `samtools`       |
| `SM2_DIR`    | shapemapper2 install directory   | unset → skipped  |

`accuracy` additionally requires **ciffy**, **biopython**, and **scipy** (cmuts
already depends on the latter two). BAMs passed to `profiles` for the
per-reference shapemapper2 loop must be **sorted and indexed**.

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

### profiles

```bash
python benchmarks/profiles.py \
  --fasta references.fasta \
  --dms-mod dms_mod.bam --dms-nomod dms_nomod.bam \
  --2a3-mod 2a3_mod.bam --2a3-nomod 2a3_nomod.bam \
  -o profiles.h5 [--downsample N]
```

Each condition (DMS, 2A3) takes a modified and an unmodified (nomod) sample, and
**each tool uses its own native pipeline** (counting + background subtraction +
normalization), so the comparison is between fully processed reactivities:

- **rf-count** + **rf-norm** (Siegfried `-sm 3`, 2-8% `-nm 1`, treated vs
  untreated).
- **shapemapper2** — its `make_reactivity_profiles` + `normalize_profiles` step on
  the modified + untreated samples (run with shapemapper's bundled Python; the
  normalization factor is computed library-wide across references).
- **cmuts** — `cmuts core` (once per spread mode) + `cmuts normalize --mod
  --nomod`. The `cmuts-match-rf` dataset uses `--no-spread` with outlier
  normalization and `--independent-norm`; `cmuts-{nospread,uniform,default}` use
  UBR normalization across the three spread modes.

The output HDF5 has one `(n_references, length)` reactivity array per dataset and
condition (`/<dataset>/<DMS|2A3>`), in FASTA order. `--downsample N` caps reads
per reference in `cmuts core`.

### correctness

```bash
python benchmarks/correctness.py --profiles profiles.h5 -f references.fasta \
  -o correctness.csv
```

Reads the profiles HDF5 and writes one row per (comparison, condition, reference,
position) with the cmuts-tuned-to-the-tool reactivity, that tool's reactivity,
and their absolute difference. With spreading off and matched settings, cmuts
should reproduce each tool's profile nucleotide-by-nucleotide; the residual
isolates deletion spreading.

### accuracy

```bash
python benchmarks/accuracy.py \
  --profiles profiles.h5 \
  --fasta references.fasta \
  --structures /path/to/cif/dir \
  --per-ref per_reference.tsv --summary summary.tsv
```

Reads each dataset's reactivity from the profiles HDF5 and scores it against a
base-pairing ground truth read from the deposited mmCIF hydrogen-bond
annotations (references whose structure carries no such annotations are skipped).
The summary reports mean/median AUC, Pearson, and Spearman per (dataset ×
condition), and each dataset's delta against cmuts' default spread mode.

## Sherlock notes

Two issues found while validating these on the cluster:

1. **Do not `source .venv/bin/activate`.** The venv activation drops samtools
   from `PATH` and rf-count then fails with "samtools is not in PATH". Instead
   source the cluster deps and call the interpreter directly:

   ```bash
   source "$CMUTS_HOME/__sherlock_deps"
   export PATH="$CMUTS_HOME/bin:/path/to/RNAFramework:$PATH"
   export SM2_DIR=/path/to/shapemapper2
   "$CMUTS_HOME/.venv/bin/python" benchmarks/profiles.py ...
   ```

   (Only `profiles.py` invokes the external tools; `correctness.py` and
   `accuracy.py` read the HDF5 it writes and need no samtools/rf-count/sm2.)

2. `profiles` calls `cmuts normalize`, which on older builds
   allocated the full SNR-scaling array and could run out of memory on inputs
   with thousands of references. This is fixed in current builds (normalization
   no longer plots); deploy a current `cmuts` before large runs.
