"""How the benchmarks drive cmuts: `cmuts core` (count) then `cmuts normalize`.

cmuts reads BAM/CRAM/SAM natively, so no input conversion is needed. `reactivity`
at the bottom is the full per-tool pipeline.
"""

from __future__ import annotations

import os
from pathlib import Path

from .common import FLOOR, CountParams, run_checked, strip_ext


def cmuts_bin() -> str:
    return os.environ.get("CMUTS", "cmuts")


# spread mode -> `cmuts core` flag (None = the mutation-informed default).
_SPREAD = {"default": None, "nospread": "--no-spread", "uniform": "--uniform-spread"}


def core_command(bams, fasta, out_h5, params: CountParams, *, downsample: int = 0) -> list[str]:
    """Build the `cmuts core` argv from `params` (counts one or more BAM/CRAM/SAM,
    which cmuts reads natively, into one HDF5)."""
    cmd = [cmuts_bin(), "core", "--threads", str(params.threads), "-f", str(fasta),
           "-o", str(out_h5), "--overwrite"]  # fmt: skip
    spread = _SPREAD[params.cmuts_spread]
    if spread:
        cmd.append(spread)
    if not params.insertions:
        cmd.append("--no-insertions")
    if params.min_mapq is not None:
        cmd += ["--min-mapq", str(params.min_mapq)]
    if params.min_phred is not None:
        cmd += ["--min-phred", str(params.min_phred)]
    if params.min_length is not None:
        cmd += ["--min-length", str(params.min_length)]
    if params.max_length is not None:
        cmd += ["--max-length", str(params.max_length)]
    if params.max_indel is not None:
        cmd += ["--max-indel-length", str(params.max_indel)]
    if params.collapse and params.collapse > 1:
        cmd += ["--collapse", str(params.collapse)]
    cmd += ["--quality-window", "1" if params.eval_surrounding else "0"]
    if downsample:
        cmd += ["--downsample", str(downsample)]
    cmd += [str(b) for b in bams]
    return cmd


def normalize(counts_h5, mod, nomod, fasta, condition, method, per_reference, out_h5):
    """`cmuts normalize` (mod vs nomod) -> (reactivity (n_ref, length), reads
    (n_ref,)). Samples are identified by BAM basename, as `cmuts core` writes them.
    `reads` is cmuts' per-reference depth (max per-position coverage)."""
    import h5py
    import numpy as np

    experiment = [condition, f"mod={strip_ext(mod)}"]
    if nomod is not None:
        experiment.append(f"nomod={strip_ext(nomod)}")
    cmd = [cmuts_bin(), "normalize", str(counts_h5), "-o", str(out_h5), "--fasta", str(fasta),
           "--norm", method, "--blank-cutoff", str(FLOOR), "--experiment", *experiment,
           "--overwrite"]  # fmt: skip
    if per_reference:
        cmd.append("--per-reference-norm")
    run_checked(cmd)
    with h5py.File(out_h5, "r") as f:
        return np.asarray(f[f"{condition}/reactivity"]), np.asarray(f[f"{condition}/reads"])


def reactivity(count, norm, mod, nomod, condition, fasta, names, wd):
    """cmuts full pipeline: core (count) -> normalize -> ({ref: reactivity},
    {ref: depth}). cmuts reads any format, so mod/nomod pass straight through."""
    wd = Path(wd)
    bams = [mod] + ([nomod] if nomod is not None else [])
    counts_h5 = wd / "cmuts-counts.h5"
    run_checked(core_command(bams, fasta, counts_h5, count))
    method = norm.cmuts_norm
    if method == "sm":
        method = "sm-dms" if condition == "DMS" else "sm-shape"
    react, reads = normalize(
        counts_h5,
        mod,
        nomod,
        fasta,
        condition,
        method,
        norm.cmuts_per_reference,
        wd / "cmuts-norm.h5",
    )
    n_react = min(len(names), react.shape[0])
    n_reads = min(len(names), len(reads))
    return (
        {names[i]: react[i] for i in range(n_react)},
        {names[i]: float(reads[i]) for i in range(n_reads)},
    )
