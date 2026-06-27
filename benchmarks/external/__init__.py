"""The unified tool layer for the benchmarks.

cmuts, rf-count (RNAFramework), and shapemapper2 read different formats and take
different flags. Each tool's *complete* usage lives in its own module
(`cmuts`, `rnaframework`, `shapemapper`), exposing one `reactivity(...) ->
({ref: reactivity}, {ref: depth})` that adapts the input format, counts, and
normalizes. This package ties them together:

    reactivity(tool, count, norm, inputs, condition, out_h5)
        run one tool's pipeline and write a unified HDF5
        (/names, /reactivity (n_ref, length), /reads (n_ref,) per-reference depth)
    build_profiles(datasets, inputs_by_condition, outdir)
        run a dict of {label: (tool, CountParams, NormParams)} recipes per condition
    read_profiles(h5) -> (names, reactivity, reads)

`__init__` only dispatches and does the shared HDF5 I/O; per-tool logic stays in
the tool modules. Count- and normalize-stage knobs live in `CountParams` /
`NormParams`; callers own the configurations (which datasets, which params).

Hamish M. Blair, 2026
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from . import cmuts, rnaframework, shapemapper
from .common import TOOLS, CountParams, Inputs, NormParams, read_fasta, to_array
from .rnaframework import available as rfcount_available

if TYPE_CHECKING:
    import numpy as np

__all__ = [
    "CountParams",
    "NormParams",
    "Inputs",
    "TOOLS",
    "reactivity",
    "build_profiles",
    "read_profiles",
    "read_fasta",
    "rfcount_available",
    "to_array",
]


def reactivity(
    tool: str,
    count: CountParams,
    norm: NormParams,
    inputs: Inputs,
    condition: str,
    out_h5,
    *,
    sm_dir: str | None = None,
    workdir=None,
) -> Path:
    """Run `tool`'s full reactivity pipeline on `inputs` and write `out_h5`.

    Reads the FASTA, dispatches to the tool module (which adapts the input format,
    counts, and normalizes), and writes a unified HDF5: `/names`, `/reactivity`
    (n_ref, length), `/reads` (n_ref,) per-reference depth. No coverage floor of
    its own -- consumers filter on `/reads`. `condition` (DMS / 2A3) names the
    experiment and selects condition-dependent knobs.
    """
    import h5py
    import numpy as np

    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="react-"))
    wd.mkdir(parents=True, exist_ok=True)
    fasta_list = read_fasta(inputs.fasta)
    names = [n for n, _ in fasta_list]
    length = max((len(s) for _, s in fasta_list), default=0)

    if tool == "cmuts":
        react_d, reads_d = cmuts.reactivity(
            count, norm, inputs.mod, inputs.nomod, condition, inputs.fasta, names, wd
        )
    elif tool == "rnaframework":
        react_d, reads_d = rnaframework.reactivity(
            count, norm, inputs.mod, inputs.nomod, inputs.fasta, wd
        )
    elif tool == "shapemapper2":
        if sm_dir is None:
            raise ValueError("shapemapper2 requires sm_dir")
        react_d, reads_d = shapemapper.reactivity(
            count, inputs.mod, inputs.nomod, condition, inputs.fasta, sm_dir, wd
        )
    else:
        raise ValueError(f"unknown tool: {tool}")

    react = to_array(react_d, names, length).astype(np.float32)
    reads = np.array([reads_d.get(n, 0.0) for n in names], dtype=np.float32)
    with h5py.File(out_h5, "w") as f:
        f.create_dataset("names", data=np.array(names, dtype="S"))
        f.create_dataset("reactivity", data=react)
        f.create_dataset("reads", data=reads)
    return Path(out_h5)


def read_profiles(h5) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Load a unified profiles HDF5 -> (names, reactivity (n_ref, length), reads
    (n_ref,)), written by `reactivity`."""
    import h5py
    import numpy as np

    with h5py.File(h5, "r") as f:
        names = [n.decode() for n in f["names"][:]]
        return names, np.asarray(f["reactivity"]), np.asarray(f["reads"])


def build_profiles(
    datasets: dict[str, tuple[str, CountParams, NormParams]],
    inputs_by_condition: dict[str, Inputs],
    outdir,
    *,
    sm_dir: str | None = None,
) -> dict[str, dict[str, Path]]:
    """Run each `(tool, count, norm)` recipe for each condition, writing one
    unified HDF5 per (label, condition) under `outdir`. Returns
    {label: {condition: h5_path}} (load with `read_profiles`).

    The caller owns `datasets` (the recipes and their labels); this is just the
    loop. Recipes whose tool is unavailable (rf-count missing, no `sm_dir`) are
    skipped, so the result holds only datasets that could be built.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, dict[str, Path]] = {}
    for label, (tool, count, norm) in datasets.items():
        if tool == "rnaframework" and not rfcount_available():
            continue
        if tool == "shapemapper2" and sm_dir is None:
            continue
        out[label] = {}
        for cond, inputs in inputs_by_condition.items():
            h5 = outdir / f"{label}__{cond}.h5"
            reactivity(
                tool,
                count,
                norm,
                inputs,
                cond,
                h5,
                sm_dir=sm_dir,
                workdir=outdir / f"_wd-{label}-{cond}",
            )
            out[label][cond] = h5
    return out
