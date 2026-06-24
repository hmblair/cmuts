#!/usr/bin/env python3
"""Generate the benchmark reactivity profiles once, for both downstream benchmarks.

Runs each tool's full pipeline (counting + background subtraction + normalization)
on the per-condition modified/unmodified BAMs and writes one HDF5, so the
correctness and accuracy benchmarks read precomputed profiles instead of
re-running the tools (which lets the structure ground truth be re-scored cheaply).

Datasets, each a (n_references, length) reactivity array per probe condition
(DMS, 2A3), in FASTA order:

    rnaframework     rf-count + rf-norm (Siegfried -sm 3, 2-8% -nm 1, treated vs untreated)
    shapemapper2     shapemapper2 native reactivity (modified vs untreated, 2-8%)
    cmuts-match-rf   cmuts --no-spread, normalize --norm outlier --per-reference-norm (mod/nomod)
    cmuts-match-sm   cmuts --no-spread, normalize --norm sm-dms (DMS) / sm-shape (2A3) (mod/nomod)
    cmuts-nospread   cmuts --no-spread, normalize --norm ubr (mod/nomod)
    cmuts-uniform    cmuts --uniform-spread, --norm ubr (mod/nomod)
    cmuts-default    cmuts (mutation-informed spread), --norm ubr (mod/nomod)

correctness consumes {rnaframework, shapemapper2, cmuts-match-rf, cmuts-match-sm};
accuracy consumes {rnaframework, shapemapper2, cmuts-nospread, cmuts-uniform,
cmuts-default}.

All tools share one map-seq counting config (`MAP`): insertions off, deletions on
and right-aligned, collapse 2, mapq/phred 10, eval-surrounding (window 1),
duplicates kept. The BAMs must be sorted and indexed. See `external` for the tool
environment variables (RF_COUNT, RF_RCTOOLS, RF_NORM, SM2_DIR, SAMTOOLS).

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import tempfile
from pathlib import Path

import external
import h5py
import numpy as np

CONDITIONS = ["DMS", "2A3"]
BLANK_CUTOFF = 10  # min reads per position (shared with the other tools' MIN_DEPTH)

# One map-seq counting config for every tool (cmuts core flags are set explicitly
# below; this drives rf-count and the shapemapper parser).
MAP = external.Params(
    insertions=False,
    deletions=True,
    right_align_deletions=True,
    collapse=2,
    eval_surrounding=True,
    cov_low_qual=True,  # do not pass rf-count --no-cov-low-qual (matches COMMANDS)
    discard_duplicates=False,  # keep duplicates: map-seq reads pile up at one coordinate
    min_mapq=10,
    min_phred=10,
    min_length=2,
    median_quality=0,
    max_edit_distance=1.0,
    max_internal_match=2,
)

# cmuts core flag per spread mode ("" = the mutation-informed default).
SPREAD = {"nospread": "--no-spread", "uniform": "--uniform-spread", "default": ""}


def _strip_ext(bam: Path) -> str:
    s = str(bam)
    for ext in (".bam", ".cram", ".sam"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


def _to_array(per_ref: dict[str, np.ndarray], names: list[str], length: int) -> np.ndarray:
    """Stack a {name: 1-D reactivity} dict into a (n_references, length) array in
    FASTA order; missing references / positions are NaN."""
    out = np.full((len(names), length), np.nan)
    for i, name in enumerate(names):
        v = per_ref.get(name)
        if v is not None:
            out[i, : len(v)] = v[:length]
    return out


# ---------------------------------------------------------------------------
# Per-tool full pipelines
# ---------------------------------------------------------------------------


def rnaframework(mod, nomod, fasta, workdir) -> dict[str, np.ndarray]:
    """rf-count both samples, then rf-norm (Siegfried) with background subtraction."""
    rc_mod = external.run_rfcount(mod, fasta, workdir / "rf-mod", MAP)
    rc_nomod = external.run_rfcount(nomod, fasta, workdir / "rf-nomod", MAP)
    norm_dir = external.run_rfnorm(
        rc_mod, workdir / "rf-norm", scoring_method=3, norm_method=1, untreated=rc_nomod
    )
    return external.parse_rfnorm(norm_dir)


def shapemapper(mod, nomod, fasta_list, sm_dir, dms, workdir) -> dict[str, np.ndarray]:
    """shapemapper2's native reactivity (modified vs untreated, normalized)."""
    return external.run_shapemapper_reactivity(
        mod, nomod, fasta_list, sm_dir, MAP, dms=dms, workdir=workdir / "sm"
    )


def cmuts_core(spread_flag: str, all_bams, fasta, counts, downsample, threads) -> None:
    """cmuts core with the shared map-seq config and a spread flag.

    `--threads N` makes the cmuts dispatcher run `mpirun -np N cmuts-core`, which
    parallelizes across references; without it cmuts-core runs as a single MPI
    rank (one core), which is ~N times slower on a many-reference dataset.
    """
    cmuts = os.environ.get("CMUTS", "cmuts")
    cmd = [
        cmuts,
        "core",
        "--threads",
        str(threads),
        "-f",
        str(fasta),
        "-o",
        str(counts),
        "--overwrite",
        "--no-insertions",
        "--min-mapq",
        "10",
        "--min-phred",
        "10",
        "--min-length",
        "2",
        "--collapse",
        "2",
        "--quality-window",
        "1",
    ]
    if spread_flag:
        cmd.append(spread_flag)
    if downsample:
        cmd += ["--downsample", str(downsample)]
    cmd += [str(b) for b in all_bams]
    external.run_checked(cmd)


def cmuts_normalize(counts, mod, nomod, fasta, cond, norm, independent, out_h5) -> np.ndarray:
    """cmuts normalize (mod vs nomod) -> the reactivity array for one condition.

    The condition name is the experiment (and output HDF5 group). `independent`
    selects per-reference normalization (rf-norm-like), matching cmuts-match-rf.
    """
    cmuts = os.environ.get("CMUTS", "cmuts")
    experiment = [cond, f"mod={_strip_ext(mod)}", f"nomod={_strip_ext(nomod)}"]
    cmd = [
        cmuts,
        "normalize",
        str(counts),
        "-o",
        str(out_h5),
        "--fasta",
        str(fasta),
        "--norm",
        norm,
        "--blank-cutoff",
        str(BLANK_CUTOFF),
        "--experiment",
        *experiment,
        "--overwrite",
    ]
    if independent:
        cmd.append("--per-reference-norm")
    external.run_checked(cmd)
    with h5py.File(out_h5, "r") as f:
        return np.asarray(f[f"{cond}/reactivity"])


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--dms-mod", dest="mod_dms", required=True, type=Path)
    p.add_argument("--dms-nomod", dest="nomod_dms", required=True, type=Path)
    p.add_argument("--2a3-mod", dest="mod_2a3", required=True, type=Path)
    p.add_argument("--2a3-nomod", dest="nomod_2a3", required=True, type=Path)
    p.add_argument("-o", "--out", required=True, type=Path, help="Output profiles HDF5")
    p.add_argument(
        "--downsample", type=int, default=0, help="Per-reference read cap for cmuts core"
    )
    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="MPI ranks for cmuts core and worker threads for rf-count",
    )
    args = p.parse_args()

    # rf-count runs single-threaded by default in MAP; give it the same budget.
    global MAP
    MAP = dataclasses.replace(MAP, threads=args.threads)

    fasta = _read_fasta(args.fasta)
    names = [n for n, _ in fasta]
    length = max(len(s) for _, s in fasta)
    cond_samples = {
        "DMS": (args.mod_dms, args.nomod_dms),
        "2A3": (args.mod_2a3, args.nomod_2a3),
    }
    sm_dir = os.environ.get("SM2_DIR")
    all_bams = [b for pair in cond_samples.values() for b in pair]

    # cmuts core once per spread mode (shared by the cmuts datasets below).
    cmuts_dir = Path(tempfile.mkdtemp(prefix="profiles-cmuts-"))
    counts = {}
    for mode, flag in SPREAD.items():
        print(f"cmuts core [{mode}] ...")
        counts[mode] = cmuts_dir / f"counts-{mode}.h5"
        cmuts_core(flag, all_bams, args.fasta, counts[mode], args.downsample, args.threads)

    # (dataset, builder(cond, mod, nomod)) -> {name: reactivity}
    datasets: dict[str, dict[str, np.ndarray]] = {}
    for cond, (mod, nomod) in cond_samples.items():
        wd = Path(tempfile.mkdtemp(prefix=f"profiles-{cond.lower()}-"))
        print(f"=== {cond} ===")

        print("  rnaframework ...")
        datasets.setdefault("rnaframework", {})[cond] = _to_array(
            rnaframework(mod, nomod, args.fasta, wd), names, length
        )
        if sm_dir:
            print("  shapemapper2 ...")
            datasets.setdefault("shapemapper2", {})[cond] = _to_array(
                shapemapper(mod, nomod, fasta, sm_dir, cond == "DMS", wd), names, length
            )

        prof = wd / "cmuts.h5"
        print("  cmuts-match-rf ...")
        datasets.setdefault("cmuts-match-rf", {})[cond] = _pad(
            cmuts_normalize(
                counts["nospread"], mod, nomod, args.fasta, cond, "outlier", True, prof
            ),
            length,
        )
        # cmuts tuned to shapemapper2: global per-base 75th-pct (DMS) or boxplot
        # (2A3) normalization, both library-wide (not per-reference).
        print("  cmuts-match-sm ...")
        sm_norm = "sm-dms" if cond == "DMS" else "sm-shape"
        datasets.setdefault("cmuts-match-sm", {})[cond] = _pad(
            cmuts_normalize(counts["nospread"], mod, nomod, args.fasta, cond, sm_norm, False, prof),
            length,
        )
        for mode, dataset in [
            ("nospread", "cmuts-nospread"),
            ("uniform", "cmuts-uniform"),
            ("default", "cmuts-default"),
        ]:
            print(f"  {dataset} ...")
            datasets.setdefault(dataset, {})[cond] = _pad(
                cmuts_normalize(counts[mode], mod, nomod, args.fasta, cond, "ubr", False, prof),
                length,
            )

    with h5py.File(args.out, "w") as f:
        f.create_dataset("names", data=np.array(names, dtype="S"))
        for dataset, by_cond in datasets.items():
            for cond, arr in by_cond.items():
                f.create_dataset(f"{dataset}/{cond}", data=arr)
    print(f"\nwrote {len(datasets)} datasets x {len(CONDITIONS)} conditions to {args.out}")


def _pad(arr: np.ndarray, length: int) -> np.ndarray:
    """Pad a (n_refs, L) cmuts array to (n_refs, length) with NaN."""
    if arr.shape[1] >= length:
        return arr[:, :length]
    out = np.full((arr.shape[0], length), np.nan)
    out[:, : arr.shape[1]] = arr
    return out


def _read_fasta(path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


if __name__ == "__main__":
    main()
