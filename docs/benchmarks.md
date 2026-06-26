## Purpose

The `benchmarks/` directory holds the scripts that compare `cmuts` against the
two other MaP-seq mutation counters, **rf-count** (RNAFramework) and
**shapemapper2**:

- **profile** — wall-clock time and peak memory of each tool's *full reactivity
  pipeline* (count + normalize) across sweeps of query count, reference count,
  and reference length. Generates its own synthetic data.
- **tool-vs-tool** — per-nucleotide agreement of cmuts (with spreading off, tuned
  to each tool) against that tool's own reactivity. This is the Fig 2a control
  that isolates deletion spreading.
- **tool-vs-structure** — how well each tool's reactivity predicts base pairing
  in deposited 3D structures (PDB130), scored as AUC / Pearson / Spearman against
  a ground truth read from the deposited mmCIF hydrogen-bond annotations.

These are development tools, not part of the installed `cmuts` package, so they
are run as plain scripts rather than `cmuts` subcommands.

## Shared tool layer

cmuts, rf-count, and shapemapper2 read different formats and take different flags.
[`benchmarks/external.py`](https://github.com/hmblair/cmuts/blob/master/benchmarks/external.py)
is the single place that knows how to *run* each of them, so no benchmark
re-implements any tool invocation, format conversion, or output parsing.

The entry point is `reactivity(tool, count, norm, inputs, condition, out_h5)`:
given an `Inputs` struct (modified + untreated alignments in any of bam/cram/sam)
and the two stage-parameter structs, it adapts each input to the format that tool
needs (cmuts reads all natively; rf-count converts to BAM; shapemapper2
sorts/indexes to an indexed BAM), counts, and normalizes, and writes a unified
HDF5 (`/names`, `/reactivity` `(n_ref, length)`, and `/reads` `(n_ref,)` —
per-reference depth, the max post-filter per-position coverage). It applies no
coverage floor of its own; filtering on `/reads` is the consumer's choice. Load
it with `read_profiles`. Two stage structs map to each tool's flags:

- **`CountParams`** — count-stage knobs (quality thresholds, collapse, cmuts
  spread mode, and a few tool-specific knobs).
- **`NormParams`** — normalize-stage knobs (rf-norm scoring/method, cmuts norm
  method / per-reference).

`external.py` holds **no** named-dataset registry — each benchmark owns its own
configurations. The two scoring scripts share one map-seq counting config (`MAP`,
defined in each): insertions off, deletions on and right-aligned, collapse 2,
mapq/phred 10, eval-surrounding (cmuts quality window 1), duplicates kept; only
the cmuts spread mode and the normalization vary across the cmuts datasets.
`external.build_profiles(datasets, inputs_by_condition, outdir, ...)` runs a dict
of `(tool, CountParams, NormParams)` recipes for each condition, writes one
unified HDF5 per (dataset, condition) under `outdir`, and returns their paths.

## Environment

The external tools are resolved from the environment, so a machine without them
still runs the cmuts-only paths (an absent tool is skipped):

| Variable     | Meaning                          | Default          |
|--------------|----------------------------------|------------------|
| `CMUTS`      | cmuts dispatcher                 | `cmuts` on PATH  |
| `RF_COUNT`   | RNAFramework `rf-count`          | `rf-count`       |
| `RF_RCTOOLS` | RNAFramework `rf-rctools`        | `rf-rctools`     |
| `RF_NORM`    | RNAFramework `rf-norm`           | `rf-norm`        |
| `SAMTOOLS`   | `samtools`                       | `samtools`       |
| `SM2_DIR`    | shapemapper2 install directory   | unset → skipped  |

All three scripts invoke the external tools (each builds the reactivity it
needs), so each needs samtools, optionally rf-count, and `SM2_DIR` for the
shapemapper2 datasets. `tool-vs-structure` additionally requires **ciffy**,
**biopython**, and **scipy** (cmuts already depends on the latter two). Input
alignments should be coordinate-sorted; `external.py` indexes/converts them per
tool as needed.

## Running

### profile

```bash
python benchmarks/profile.py --references --runs 3 --threads 1 8 \
  --format bam,cram --shapemapper-dir /path/to/shapemapper2 \
  --output profile-results.csv
```

Sweeps are selected with `--queries`, `--references`, `--lengths` (at least one
required). Synthetic modified + untreated alignments are generated with `cmuts
generate` and cached under `./.cases`; each tool is timed on each `--format`, and
any per-tool format conversion is part of the measured pipeline. rf-count is
skipped if not on PATH; shapemapper2 runs only when `--shapemapper-dir` is given.
Results are appended to the CSV; the figure sources read it.

### tool-vs-tool

```bash
python benchmarks/tool-vs-tool.py \
  -f references.fasta \
  --dms-mod dms_mod.bam --dms-nomod dms_nomod.bam \
  --2a3-mod 2a3_mod.bam --2a3-nomod 2a3_nomod.bam \
  --threads 8 -o tool-vs-tool.csv
```

Builds each tool's reactivity from the alignments via `external.reactivity` and
writes one row per (comparison, condition, reference, position) with the
cmuts-tuned-to-the-tool reactivity, that tool's reactivity, and their absolute
difference. With spreading off and matched settings, cmuts should reproduce each
tool's profile nucleotide-by-nucleotide; the residual isolates deletion
spreading. The `cmuts-match-rf` vs `rnaframework` comparison needs rf-count; the
`cmuts-match-sm` vs `shapemapper2` comparison needs `SM2_DIR` (skipped if unset).

### tool-vs-structure

```bash
python benchmarks/tool-vs-structure.py \
  --fasta references.fasta \
  --structures /path/to/cif/dir \
  --dms-mod dms_mod.bam --dms-nomod dms_nomod.bam \
  --2a3-mod 2a3_mod.bam --2a3-nomod 2a3_nomod.bam \
  --threads 8 --per-ref per_reference.tsv --summary summary.tsv
```

Builds each dataset's reactivity (cmuts at each spread mode, plus rf-count and
shapemapper2) from the alignments and scores it against a base-pairing ground
truth read from the deposited mmCIF hydrogen-bond annotations (references whose
structure carries no such annotations are skipped). The summary reports
mean/median AUC, Pearson, and Spearman per (dataset × condition), and each
dataset's delta against cmuts' default spread mode.

## Sherlock notes

Two issues found while validating these on the cluster:

1. **Do not `source .venv/bin/activate`.** The venv activation drops samtools
   from `PATH` and rf-count then fails with "samtools is not in PATH". Instead
   source the cluster deps and call the interpreter directly:

   ```bash
   source "$CMUTS_HOME/__sherlock_deps"
   export PATH="$CMUTS_HOME/bin:/path/to/RNAFramework:$PATH"
   export SM2_DIR=/path/to/shapemapper2
   "$CMUTS_HOME/.venv/bin/python" benchmarks/tool-vs-tool.py ...
   ```

   All three scripts invoke the external tools, so each needs samtools (and
   rf-count / `SM2_DIR` for those datasets) on `PATH`.

2. The scoring scripts call `cmuts normalize`, which on older builds allocated
   the full SNR-scaling array and could run out of memory on inputs with
   thousands of references. This is fixed in current builds (normalization no
   longer plots); deploy a current `cmuts` before large runs.
