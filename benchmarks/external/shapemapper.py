"""How the benchmarks drive shapemapper2.

shapemapper's counter is single-reference (one `-n length`, no RNAME), so each
reference is extracted from an indexed BAM/CRAM and parsed+counted on its own
(`count_references`); the per-reference profiles are then background-subtracted
and normalized library-wide with shapemapper's bundled Python (`_normalized_
profiles`). `reactivity` at the bottom is the full per-tool pipeline.

shapemapper has no PATH convention, so its install dir (`sm_dir`) is passed in.
Its prebuilt binaries link the vendored boost, scoped onto LD_LIBRARY_PATH by
`env()`; its reactivity scripts need the bundled Python 3.9 (`python()`).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from .common import (
    FLOOR,
    CountParams,
    covered_references,
    ensure_index,
    max_combined,
    read_fasta,
    run_checked,
    samtools_bin,
)

if TYPE_CHECKING:
    import numpy as np


def paths(sm_dir: str) -> tuple[str, str]:
    """(parser, counter) binary paths for a shapemapper2 install dir."""
    base = f"{sm_dir}/internals/bin"
    return f"{base}/shapemapper_mutation_parser", f"{base}/shapemapper_mutation_counter"


def _libdir(sm_dir: str) -> str:
    return f"{sm_dir}/internals/thirdparty/miniconda/envs/shapemapper_make/lib"


def python(sm_dir: str) -> str:
    """shapemapper's bundled Python 3.9 (its reactivity scripts use `open(..., "rU")`)."""
    return f"{sm_dir}/internals/thirdparty/miniconda/envs/shapemapper_make/bin/python"


def env(sm_dir: str) -> dict[str, str]:
    """Environment with the vendored boost on LD_LIBRARY_PATH, scoped to these
    binaries (exporting it process-wide shadows system libs and breaks samtools)."""
    e = dict(os.environ)
    libdir = _libdir(sm_dir)
    if os.path.isdir(libdir):
        e["LD_LIBRARY_PATH"] = libdir + ":" + e.get("LD_LIBRARY_PATH", "")
    return e


def parser_flags(params: CountParams) -> list[str]:
    flags = ["--input_is_unpaired"]
    if params.min_mapq is not None:
        flags += ["-m", str(params.min_mapq)]
    if params.min_phred is not None:
        flags += ["--min_qual", str(params.min_phred)]
    if params.max_internal_match is not None:
        flags += ["--max_internal_match", str(params.max_internal_match)]
    flags.append("--right_align_ambig_dels")  # always right-align ambiguous deletions
    return flags


def parse_profile(path, length: int) -> np.ndarray:
    """Read the reactivity column from a shapemapper profile / normalized file
    (`Norm_profile`, else `Reactivity_profile`); non-finite entries become NaN."""
    import numpy as np

    out = np.full(length, np.nan)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = next(
            (header.index(c) for c in ("Norm_profile", "Reactivity_profile") if c in header), None
        )
        if col is None:
            return out
        for pos, line in enumerate(f):
            if pos >= length:
                break
            cells = line.rstrip("\n").split("\t")
            if col < len(cells):
                try:
                    v = float(cells[col])
                except ValueError:
                    v = float("nan")
                out[pos] = v if np.isfinite(v) else np.nan
    return out


def _counts_depth(path) -> np.ndarray:
    """Per-position `effective_depth` from a shapemapper counts file."""
    import numpy as np

    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "effective_depth" not in header:
            return np.zeros(0)
        idx = header.index("effective_depth")
        vals: list[float] = []
        for line in f:
            cells = line.rstrip("\n").split("\t")
            if idx < len(cells):
                try:
                    vals.append(float(cells[idx]))
                except ValueError:
                    vals.append(0.0)
    return np.asarray(vals)


def _count_reference(
    bam,
    name,
    length,
    *,
    parser,
    counter,
    params,
    env,
    samtools,
    workdir,
    counts,
    sam_prefix,
    reference=None,
) -> Path | None:
    """Extract one reference's reads and run the parser then counter.

    Counts go to `counts` (kept per reference so they coexist); the SAM/`.mut`
    intermediates (named from `sam_prefix`) are transient. `reference` is the FASTA
    for `samtools view -T` so CRAM decodes (ignored for BAM). None on failure.
    """
    sam = workdir / f"{sam_prefix}.sam"
    parsed = workdir / f"{sam_prefix}.mut"
    counts = Path(counts)
    view_cmd = [samtools, "view", "-h", "-F", "4"]
    if reference is not None:
        view_cmd += ["-T", str(reference)]
    with open(sam, "w") as fh:
        view = subprocess.run([*view_cmd, str(bam), name], stdout=fh)
    if view.returncode != 0:
        return None
    subprocess.run(
        [parser, "-i", str(sam), "-o", str(parsed), *parser_flags(params)],
        env=env,
        capture_output=True,
    )
    cp = subprocess.run(
        [counter, "-i", str(parsed), "-c", str(counts), "-n", str(length), "--warn_on_no_mapped"],
        env=env,
        capture_output=True,
    )
    return counts if (cp.returncode == 0 and counts.exists()) else None


def count_references(
    bam,
    fasta_list,
    sm_dir,
    params: CountParams,
    *,
    references=None,
    reference=None,
    samtools=None,
    min_depth=FLOOR,
    workdir,
    prefix="sm",
) -> dict[str, Path]:
    """Per-reference parser+counter over a multi-reference indexed BAM/CRAM, one
    `-n length` invocation per covered reference. Returns {ref_name: counts_file}
    (distinct files, so all coexist). `references` overrides the idxstats coverage
    pre-filter (used to count a control over its modified mate's references)."""
    parser, counter = paths(sm_dir)
    st = samtools or samtools_bin()
    e = env(sm_dir)
    covered = (
        references
        if references is not None
        else covered_references(bam, min_reads=min_depth, samtools=st)
    )
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for idx, (name, seq) in enumerate(fasta_list):
        if name not in covered:
            continue
        counts = _count_reference(
            bam,
            name,
            len(seq),
            parser=parser,
            counter=counter,
            params=params,
            env=e,
            samtools=st,
            workdir=workdir,
            counts=workdir / f"{prefix}-{idx}-counts.txt",
            sam_prefix=prefix,
            reference=reference,
        )
        if counts is not None:
            out[name] = counts
    return out


def _normalized_profiles(
    mod_counts, nomod_counts, fasta_list, sm_dir, dms, wd
) -> dict[str, np.ndarray]:
    """Background-subtract each reference (`make_reactivity_profiles.py`), normalize
    all together with one library-wide factor (`normalize_profiles.py`), and read
    each reference's profile (normalized if produced, else the raw subtraction)."""
    py = python(sm_dir)
    e = env(sm_dir)
    make_prof = f"{sm_dir}/internals/bin/make_reactivity_profiles.py"
    norm_prof = f"{sm_dir}/internals/bin/normalize_profiles.py"
    (wd / "prof").mkdir(parents=True, exist_ok=True)
    (wd / "norm").mkdir(parents=True, exist_ok=True)
    lengths = {name: len(seq) for name, seq in fasta_list}

    order: list[str] = []
    for name, seq in fasta_list:
        if name not in mod_counts:
            continue
        count_files = [str(mod_counts[name])]
        if name in nomod_counts:
            count_files.append(str(nomod_counts[name]))
        ref_fa = wd / "ref.fa"
        ref_fa.write_text(f">{name}\n{seq}\n")
        prof_path = wd / "prof" / f"{len(order)}.txt"
        cmd = [py, make_prof, "--fa", str(ref_fa), "--counts", *count_files,
               "--out", str(prof_path), "--mindepth", str(FLOOR)]  # fmt: skip
        if dms:
            cmd.append("--dms")
        subprocess.run(cmd, env=e, capture_output=True)
        if prof_path.exists():
            order.append(name)

    tonorm = [str(wd / "prof" / f"{i}.txt") for i in range(len(order))]
    normout = [str(wd / "norm" / f"{i}.txt") for i in range(len(order))]
    if tonorm:
        ncmd = [py, norm_prof, "--tonorm", *tonorm, "--normout", *normout, "--warn-on-error"]
        if dms:
            ncmd.append("--dms")
        subprocess.run(ncmd, env=e, capture_output=True)

    out: dict[str, np.ndarray] = {}
    for i, name in enumerate(order):
        normed = Path(normout[i])
        src = normed if normed.exists() else Path(tonorm[i])
        out[name] = parse_profile(src, lengths[name])
    return out


def _to_indexed(path, fasta, wd, samtools) -> Path:
    """shapemapper needs an indexed BAM/CRAM (a SAM is sorted+indexed to BAM)."""
    path = Path(path)
    if path.suffix.lower() in (".bam", ".cram"):
        ensure_index(path, samtools)
        return path
    out = Path(wd) / f"{path.stem}.bam"
    run_checked([samtools, "sort", "-o", str(out), str(path)])
    run_checked([samtools, "index", str(out)])
    return out


def reactivity(count, mod, nomod, condition, fasta, sm_dir, wd):
    """shapemapper2 full pipeline: ensure indexed inputs -> per-reference
    parser+counter (mod and untreated over the same references) -> background-
    subtract + normalize -> ({ref: reactivity}, {ref: depth}). Depth is the
    combined (mod + untreated) max per-position effective_depth. Returns ({}, {})
    if the install lacks the reactivity scripts."""
    wd = Path(wd)
    if not os.path.exists(f"{sm_dir}/internals/bin/make_reactivity_profiles.py"):
        return {}, {}
    st = samtools_bin()
    fasta_list = read_fasta(fasta)
    mod = _to_indexed(mod, fasta, wd, st)
    nomod = _to_indexed(nomod, fasta, wd, st) if nomod is not None else None

    mod_counts = count_references(
        mod, fasta_list, sm_dir, count, reference=str(fasta), samtools=st, workdir=wd, prefix="mod"
    )
    nomod_counts = (
        count_references(
            nomod,
            fasta_list,
            sm_dir,
            count,
            references=set(mod_counts),
            reference=str(fasta),
            samtools=st,
            workdir=wd,
            prefix="nomod",
        )
        if nomod is not None
        else {}
    )
    react = _normalized_profiles(
        mod_counts, nomod_counts, fasta_list, sm_dir, condition == "DMS", wd
    )
    reads = {
        name: max_combined(
            [_counts_depth(p)]
            + ([_counts_depth(nomod_counts[name])] if name in nomod_counts else [])
        )
        for name, p in mod_counts.items()
    }
    return react, reads
