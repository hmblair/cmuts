"""The unified tool layer for the benchmarks.

cmuts, rf-count (RNAFramework), and shapemapper2 read different formats and take
different flags. This module is the single place that knows how to *run* each
tool's full reactivity pipeline and return a uniform per-reference array, so the
timing and scoring benchmarks re-implement no tool invocation or format
conversion.

The entry point is `reactivity(tool, count, norm, inputs, condition)`: it adapts
the inputs to whatever format `tool` needs, counts, and normalizes. Count- and
normalize-stage knobs live in `CountParams` / `NormParams` and are mapped to each
tool's flags here; callers own the actual configurations (which datasets, which
params).

Tools are resolved from the environment so a machine without them still runs the
cmuts-only paths:

    CMUTS       cmuts          (default: "cmuts" on PATH)
    RF_COUNT    rf-count       (default: "rf-count" on PATH)
    RF_RCTOOLS  rf-rctools     (default: "rf-rctools" on PATH)
    RF_NORM     rf-norm        (default: "rf-norm" on PATH)
    SAMTOOLS    samtools       (default: "samtools" on PATH)

shapemapper2 has no PATH convention upstream, so its install dir is passed in.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

# numpy is imported lazily in the functions that build arrays, so paths that only
# build/run commands do not pay its import cost. The type annotations still refer
# to np.ndarray, resolved here for type-checkers/linters without a runtime import.
if TYPE_CHECKING:
    import numpy as np

# A position needs at least this many reads for a rate to be meaningful (used by
# the standalone rate parsers).
MIN_DEPTH = 10

# `reactivity` imposes no coverage floor of its own -- it only blanks zero-coverage
# positions, the most permissive each tool allows -- and exposes per-reference
# depth as `/reads` so consumers filter as they choose.
_FLOOR = 1


@dataclass
class CountParams:
    """Count-stage knobs, mapped to each tool's count flags by the builders below.

    Fields defaulting to None omit their flag (the tool's own default applies);
    fields with concrete defaults are always emitted, so `CountParams()` is a
    sensible baseline. A few fields apply to only one tool (noted inline) and are
    ignored by the others.
    """

    insertions: bool = True
    deletions: bool = True
    right_align_deletions: bool = True
    collapse: int = 2  # >1 -> collapse mods within (collapse - 1) of each other
    eval_surrounding: bool = True
    cov_low_qual: bool = False  # count low-quality positions toward coverage
    discard_duplicates: bool = True
    fast: bool = False  # rf-count --fast
    threads: int = 1
    min_mapq: int | None = None
    min_phred: int | None = None
    min_length: int | None = None
    max_length: int | None = None  # cmuts --max-length (max read length)
    max_indel: int | None = None
    quality_window: int | None = None  # cmuts --quality-window
    cmuts_spread: str = "default"  # cmuts core: default | nospread | uniform
    max_edit_distance: float | None = None  # rf-count only
    median_quality: int | None = None  # rf-count only
    max_internal_match: int | None = None  # shapemapper parser only


@dataclass
class NormParams:
    """Normalize-stage knobs, mapped to each tool's normalize step.

    cmuts uses `cmuts_norm`/`cmuts_per_reference`; rf-count uses
    `rf_scoring`/`rf_norm_method`; shapemapper2 normalizes internally (its only
    knob, DMS vs SHAPE, is derived from the condition). `cmuts_norm == "sm"`
    resolves to `sm-dms`/`sm-shape` by condition.
    """

    cmuts_norm: str = "ubr"  # ubr | outlier | sm-dms | sm-shape | sm
    cmuts_per_reference: bool = False
    rf_scoring: int = 3  # rf-norm -sm
    rf_norm_method: int = 1  # rf-norm -nm


@dataclass
class Inputs:
    """Alignment inputs for a reactivity run (any of bam/cram/sam)."""

    fasta: Path
    mod: Path
    nomod: Path | None = None


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def run_checked(cmd: list[str], **kwargs) -> None:
    """Run a command, raising with captured output on a non-zero exit."""
    proc = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}\n"
            f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}"
        )


def samtools_bin() -> str:
    return os.environ.get("SAMTOOLS", "samtools")


def _strip_ext(path) -> str:
    """Drop a trailing .bam/.cram/.sam extension (cmuts sample identifiers)."""
    s = str(path)
    for ext in (".bam", ".cram", ".sam"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


def covered_references(bam, *, min_reads: int, samtools: str | None = None) -> set[str]:
    """References in `bam` with at least `min_reads` mapped reads (from idxstats).

    A position's depth cannot exceed a reference's total reads, so references
    below this floor cannot reach `min_reads` depth anywhere and can be skipped.
    Requires an indexed BAM.
    """
    idx = subprocess.run(
        [samtools or samtools_bin(), "idxstats", str(bam)], capture_output=True, text=True
    )
    return {
        cols[0]
        for cols in (ln.split("\t") for ln in idx.stdout.splitlines())
        if len(cols) >= 3 and cols[2].isdigit() and int(cols[2]) >= min_reads
    }


def _max_combined(arrays) -> float:
    """Max over positions of the element-wise sum of per-position depth arrays --
    the combined (mod + untreated) reference depth, matching cmuts' `reads`."""
    import numpy as np

    arrays = [np.asarray(a, dtype=float) for a in arrays if len(a)]
    if not arrays:
        return 0.0
    total = np.zeros(max(len(a) for a in arrays))
    for a in arrays:
        total[: len(a)] += a
    return float(total.max())


# ---------------------------------------------------------------------------
# cmuts
# ---------------------------------------------------------------------------


def cmuts_bin() -> str:
    return os.environ.get("CMUTS", "cmuts")


# cmuts core spread mode -> flag ("" / None = the mutation-informed default).
_CMUTS_SPREAD = {"default": None, "nospread": "--no-spread", "uniform": "--uniform-spread"}


def cmuts_core_command(
    bams, fasta, out_h5, params: CountParams, *, downsample: int = 0
) -> list[str]:
    """Build the `cmuts core` argv from `params` (counts one or more BAM/CRAM/SAM,
    which cmuts reads natively, into one HDF5)."""
    cmd = [
        cmuts_bin(),
        "core",
        "--threads",
        str(params.threads),
        "-f",
        str(fasta),
        "-o",
        str(out_h5),
        "--overwrite",
    ]
    spread = _CMUTS_SPREAD[params.cmuts_spread]
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
    if params.quality_window is not None:
        cmd += ["--quality-window", str(params.quality_window)]
    if downsample:
        cmd += ["--downsample", str(downsample)]
    cmd += [str(b) for b in bams]
    return cmd


def run_cmuts_normalize(
    counts_h5,
    mod,
    nomod,
    fasta,
    condition: str,
    norm: str,
    per_reference: bool,
    out_h5,
):
    """`cmuts normalize` (mod vs nomod) -> (reactivity, reads) for the experiment.

    `reactivity` is the (n_ref, length) array; `reads` is the (n_ref,) per-
    reference depth cmuts reports (the max per-position coverage). `condition`
    names the experiment and output HDF5 group; `per_reference` selects rf-norm-
    like per-reference normalization. Samples are identified by their BAM
    basename (without extension), as written by `cmuts core`.
    """
    import h5py
    import numpy as np

    experiment = [condition, f"mod={_strip_ext(mod)}"]
    if nomod is not None:
        experiment.append(f"nomod={_strip_ext(nomod)}")
    cmd = [
        cmuts_bin(),
        "normalize",
        str(counts_h5),
        "-o",
        str(out_h5),
        "--fasta",
        str(fasta),
        "--norm",
        norm,
        "--blank-cutoff",
        str(_FLOOR),
        "--experiment",
        *experiment,
        "--overwrite",
    ]
    if per_reference:
        cmd.append("--per-reference-norm")
    run_checked(cmd)
    with h5py.File(out_h5, "r") as f:
        return np.asarray(f[f"{condition}/reactivity"]), np.asarray(f[f"{condition}/reads"])


# ---------------------------------------------------------------------------
# rf-count (RNAFramework)
# ---------------------------------------------------------------------------


def rfcount_available() -> bool:
    return (
        shutil.which(os.environ.get("RF_COUNT", "rf-count")) is not None
        and shutil.which(os.environ.get("RF_RCTOOLS", "rf-rctools")) is not None
    )


def rfcount_command(
    bam, fasta, outdir, params: CountParams, *, overwrite: bool = True
) -> list[str]:
    """Build the rf-count argv from `params`. Long flags throughout; this is the
    single home for the flag spelling the profile and correctness benchmarks
    both relied on.
    """
    rf = os.environ.get("RF_COUNT", "rf-count")
    cmd = [rf, "-m", "-f", str(fasta), "-o", str(outdir), "-wt", str(params.threads)]
    if overwrite:
        cmd.append("--overwrite")
    if params.fast:
        cmd.append("--fast")
    if params.eval_surrounding:
        cmd.append("--eval-surrounding")
    if not params.insertions:
        cmd.append("--no-insertions")
    if not params.deletions:
        cmd.append("--no-deletions")
    if params.right_align_deletions:
        cmd.append("--right-deletion")
    if not params.cov_low_qual:
        cmd.append("--no-cov-low-qual")
    if not params.discard_duplicates:
        cmd.append("--no-discard-duplicates")
    if params.collapse and params.collapse > 1:
        cmd += ["--collapse-consecutive", "--max-collapse-distance", str(params.collapse - 1)]
    if params.min_mapq is not None:
        cmd += ["--map-quality", str(params.min_mapq)]
    if params.min_phred is not None:
        cmd += ["--min-quality", str(params.min_phred)]
    if params.min_length is not None:
        cmd += ["--discard-shorter", str(params.min_length)]
    if params.max_indel is not None:
        cmd += ["--max-deletion-len", str(params.max_indel)]
    if params.max_edit_distance is not None:
        cmd += ["--max-edit-distance", str(params.max_edit_distance)]
    if params.median_quality is not None:
        cmd += ["--median-quality", str(params.median_quality)]
    cmd.append(str(bam))
    return cmd


def run_rfcount(bam, fasta, outdir, params: CountParams) -> Path:
    """Run rf-count and return the path to the resulting .rc file."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_checked(rfcount_command(bam, fasta, outdir, params))
    base = Path(bam).name
    base = base[:-4] if base.endswith(".bam") else base
    return outdir / f"{base}.rc"


def parse_rfcount(rc_path, *, min_depth: int = MIN_DEPTH) -> dict[str, np.ndarray]:
    """Per-position mutation rate per reference from an rf-count .rc file.

    Dumps the binary .rc with `rf-rctools view` and parses its four-line blocks.
    """
    rc = os.environ.get("RF_RCTOOLS", "rf-rctools")
    proc = subprocess.run([rc, "view", str(rc_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"rf-rctools view failed: {proc.stderr[-2000:]}")
    return _parse_rctools_view(proc.stdout, min_depth=min_depth)


def _parse_rctools_view(text: str, *, min_depth: int = MIN_DEPTH) -> dict[str, np.ndarray]:
    """Parse `rf-rctools view` output into {name: rate}.

    Each transcript is four lines: id, sequence, comma-separated raw counts,
    comma-separated coverage. rate = count / coverage (NaN below min_depth).
    """
    import numpy as np

    rates: dict[str, np.ndarray] = {}
    lines = [ln for ln in text.splitlines() if ln.strip() != ""]
    i = 0
    while i + 4 <= len(lines):
        name = lines[i].strip()
        counts = np.array([float(x) for x in lines[i + 2].split(",") if x != ""])
        coverage = np.array([float(x) for x in lines[i + 3].split(",") if x != ""])
        with np.errstate(divide="ignore", invalid="ignore"):
            rates[name] = np.where(coverage >= min_depth, counts / coverage, np.nan)
        i += 4
    return rates


def rfcount_coverage(rc_path) -> dict[str, np.ndarray]:
    """Per-position coverage per reference from an rf-count .rc (rf-rctools view)."""
    import numpy as np

    rc = os.environ.get("RF_RCTOOLS", "rf-rctools")
    proc = subprocess.run([rc, "view", str(rc_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"rf-rctools view failed: {proc.stderr[-2000:]}")
    out: dict[str, np.ndarray] = {}
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip() != ""]
    i = 0
    while i + 4 <= len(lines):
        out[lines[i].strip()] = np.array([float(x) for x in lines[i + 3].split(",") if x != ""])
        i += 4
    return out


def rfnorm_available() -> bool:
    return shutil.which(os.environ.get("RF_NORM", "rf-norm")) is not None


def run_rfnorm(
    rc_path,
    outdir,
    *,
    scoring_method: int = 4,
    norm_method: int = 1,
    untreated=None,
) -> Path:
    """Normalize an rf-count .rc with rf-norm and return the output directory.

    Defaults reproduce the reference pipeline: mutational-profiling scoring
    (-sm 4, "Zubradt") with 2-8% normalization (-nm 1) and coverage filters off.
    That scoring is treated-only; `untreated` is for the Ding/Siegfried methods.
    """
    rfn = os.environ.get("RF_NORM", "rf-norm")
    outdir = Path(outdir)
    cmd = [
        rfn,
        "--treated",
        str(rc_path),
        "-sm",
        str(scoring_method),
        "-nm",
        str(norm_method),
        "--mean-coverage",
        "0",
        "--median-coverage",
        "0",
        "--output-dir",
        str(outdir),
        "--overwrite",
    ]
    if untreated is not None:
        cmd += ["--untreated", str(untreated)]
    run_checked(cmd)
    return outdir


def parse_rfnorm(norm_dir) -> dict[str, np.ndarray]:
    """Parse an rf-norm output directory into {name: reactivity}.

    rf-norm writes one <id>.xml per transcript (id = FASTA header); the
    <reactivity> element is comma-separated with NaN for unscored positions.
    """
    import numpy as np

    rates: dict[str, np.ndarray] = {}
    for xml in Path(norm_dir).glob("*.xml"):
        root = ET.parse(xml).getroot()
        for transcript in root.findall("transcript"):
            react = transcript.find("reactivity")
            if react is None or not react.text:
                continue
            vals = [v.strip() for v in react.text.split(",") if v.strip() != ""]
            rates[xml.stem] = np.array(vals, dtype=float)
    return rates


# ---------------------------------------------------------------------------
# shapemapper2
# ---------------------------------------------------------------------------

# Counter columns: cols 1-30 are mutation categories, then read_depth,
# effective_depth, off_target_mapped_depth, low_mapq_mapped_depth, mapped_depth.
# The rate is (selected mutations) / effective_depth. The depth columns are never
# mutations; the insertion columns are excluded when params.insertions is False.
_SM_DEPTH_COLS = {
    "read_depth",
    "effective_depth",
    "off_target_mapped_depth",
    "low_mapq_mapped_depth",
    "mapped_depth",
}
_SM_INSERTION_COLS = {"-A", "-T", "-G", "-C", "-N", "multinuc_insertion", "complex_insertion"}


def shapemapper_paths(sm_dir: str) -> tuple[str, str]:
    """(parser, counter) binary paths for a shapemapper2 install dir."""
    base = f"{sm_dir}/internals/bin"
    return f"{base}/shapemapper_mutation_parser", f"{base}/shapemapper_mutation_counter"


def _shapemapper_libdir(sm_dir: str) -> str:
    return f"{sm_dir}/internals/thirdparty/miniconda/envs/shapemapper_make/lib"


def shapemapper_python(sm_dir: str) -> str:
    """shapemapper's bundled Python 3.9. Its reactivity/normalize scripts use
    `open(..., "rU")`, which a modern system Python no longer accepts."""
    return f"{sm_dir}/internals/thirdparty/miniconda/envs/shapemapper_make/bin/python"


def shapemapper_env(sm_dir: str) -> dict[str, str]:
    """Environment with the vendored boost on LD_LIBRARY_PATH.

    shapemapper's prebuilt binaries link against the boost shipped in its
    vendored conda env. Scope it to these binaries: exporting it process-wide
    shadows system libs (libattr/libacl) and breaks samtools/coreutils.
    """
    env = dict(os.environ)
    libdir = _shapemapper_libdir(sm_dir)
    if os.path.isdir(libdir):
        env["LD_LIBRARY_PATH"] = libdir + ":" + env.get("LD_LIBRARY_PATH", "")
    return env


def _sm_parser_flags(params: CountParams) -> list[str]:
    flags = ["--input_is_unpaired"]
    if params.min_mapq is not None:
        flags += ["-m", str(params.min_mapq)]
    if params.min_phred is not None:
        flags += ["--min_qual", str(params.min_phred)]
    if params.max_internal_match is not None:
        flags += ["--max_internal_match", str(params.max_internal_match)]
    if params.right_align_deletions:
        flags.append("--right_align_ambig_dels")
    return flags


def read_fasta(path) -> list[tuple[str, str]]:
    """Parse a FASTA into [(name, sequence)], the name taken up to first whitespace."""
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


def parse_shapemapper_counts(
    path, length: int, params: CountParams, *, min_depth: int = MIN_DEPTH
) -> np.ndarray:
    """Per-position mutation rate from a shapemapper2 counts file.

    rate = (selected mutation categories) / effective_depth. Insertion columns
    are excluded when params.insertions is False (to match cmuts --no-insertions).
    """
    import numpy as np

    rate = np.full(length, np.nan)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if "effective_depth" not in header:
            return rate
        depth_idx = header.index("effective_depth")
        skip = set(_SM_DEPTH_COLS)
        if not params.insertions:
            skip |= _SM_INSERTION_COLS
        mut_idx = [i for i, name in enumerate(header) if name not in skip]
        for pos, line in enumerate(f):
            if pos >= length:
                break
            cells = line.rstrip("\n").split("\t")
            if len(cells) <= depth_idx:
                continue
            depth = float(cells[depth_idx])
            if depth >= min_depth:
                rate[pos] = sum(float(cells[i]) for i in mut_idx) / depth
    return rate


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
    name: str,
    length: int,
    *,
    parser: str,
    counter: str,
    params: CountParams,
    env: dict[str, str],
    samtools: str,
    workdir: Path,
    counts: Path,
    sam_prefix: str = "sm",
    reference: str | None = None,
) -> Path | None:
    """Extract one reference's reads and run the parser then counter.

    Counts go to `counts` (kept per reference so they coexist); the SAM/`.mut`
    intermediates (named from `sam_prefix`) are transient. `reference` is the
    FASTA passed to `samtools view -T` so CRAM input decodes (ignored for BAM).
    Returns the counts file, or None on failure.
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
        [parser, "-i", str(sam), "-o", str(parsed), *_sm_parser_flags(params)],
        env=env,
        capture_output=True,
    )
    cp = subprocess.run(
        [counter, "-i", str(parsed), "-c", str(counts), "-n", str(length), "--warn_on_no_mapped"],
        env=env,
        capture_output=True,
    )
    return counts if (cp.returncode == 0 and counts.exists()) else None


def run_shapemapper_counts(
    bam,
    fasta: list[tuple[str, str]],
    sm_dir: str,
    params: CountParams,
    *,
    references: set[str] | None = None,
    reference: str | None = None,
    samtools: str | None = None,
    min_depth: int = MIN_DEPTH,
    workdir=None,
    prefix: str = "sm",
) -> dict[str, Path]:
    """Per-reference shapemapper parser+counter over a multi-reference BAM/CRAM.

    The counter is single-reference (one `-n length`, no RNAME), so each covered
    reference is extracted to its own SAM and counted separately. Returns
    {ref_name: counts_file}, each distinct so all coexist. References below
    `min_depth` reads are skipped via an idxstats pre-filter; `references`
    overrides it (to count a control over its modified mate's references).
    `reference` is the FASTA for `samtools view -T` (CRAM decode; ignored for
    BAM). Requires an indexed BAM/CRAM.
    """
    parser, counter = shapemapper_paths(sm_dir)
    st = samtools or samtools_bin()
    env = shapemapper_env(sm_dir)
    covered = (
        references
        if references is not None
        else covered_references(bam, min_reads=min_depth, samtools=st)
    )

    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="sm-counts-"))
    wd.mkdir(parents=True, exist_ok=True)

    out: dict[str, Path] = {}
    for idx, (name, seq) in enumerate(fasta):
        if name not in covered:
            continue
        counts = _count_reference(
            bam,
            name,
            len(seq),
            parser=parser,
            counter=counter,
            params=params,
            env=env,
            samtools=st,
            workdir=wd,
            counts=wd / f"{prefix}-{idx}-counts.txt",
            sam_prefix=prefix,
            reference=reference,
        )
        if counts is not None:
            out[name] = counts
    return out


def run_shapemapper_rates(
    bam,
    fasta: list[tuple[str, str]],
    sm_dir: str,
    params: CountParams,
    *,
    reference: str | None = None,
    samtools: str | None = None,
    min_depth: int = MIN_DEPTH,
    workdir=None,
) -> dict[str, np.ndarray]:
    """Per-position raw mutation rate per reference from shapemapper2.

    Thin wrapper over `run_shapemapper_counts`: count each covered reference,
    then read its per-position rate. (For the normalized reactivity profile, see
    `run_shapemapper_reactivity`.)
    """
    lengths = {name: len(seq) for name, seq in fasta}
    counts = run_shapemapper_counts(
        bam,
        fasta,
        sm_dir,
        params,
        reference=reference,
        samtools=samtools,
        min_depth=min_depth,
        workdir=workdir,
    )
    return {
        name: parse_shapemapper_counts(path, lengths[name], params, min_depth=min_depth)
        for name, path in counts.items()
    }


def parse_shapemapper_profile(path, length: int) -> np.ndarray:
    """Read the reactivity column from a shapemapper profile / normalized file.

    Prefers `Norm_profile` (normalized) and falls back to `Reactivity_profile`
    (background-subtracted but un-normalized). Non-finite entries become NaN.
    """
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


def run_shapemapper_reactivity(
    mod_bam,
    nomod_bam,
    fasta: list[tuple[str, str]],
    sm_dir: str,
    params: CountParams,
    *,
    dms: bool = False,
    reference: str | None = None,
    samtools: str | None = None,
    min_depth: int = MIN_DEPTH,
    workdir=None,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """shapemapper2 native reactivity per reference (mod + untreated control).

    Per reference: parser+counter on the modified and (if given) untreated
    samples, then `make_reactivity_profiles.py` to background-subtract. All
    per-reference profiles are then normalized together in a single
    `normalize_profiles.py` call -- its factor is computed across the whole set,
    since one short reference is too sparse to normalize alone. Returns
    ({ref: normalized reactivity}, {ref: depth}), the reactivity falling back per
    reference to the un-normalized Reactivity_profile where normalization could
    not be applied, and depth being the combined (mod + untreated) max
    per-position effective_depth.
    """
    env = shapemapper_env(sm_dir)
    py = shapemapper_python(sm_dir)
    make_prof = f"{sm_dir}/internals/bin/make_reactivity_profiles.py"
    norm_prof = f"{sm_dir}/internals/bin/normalize_profiles.py"
    st = samtools or samtools_bin()
    if not os.path.exists(make_prof):
        return {}, {}

    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="sm-react-"))
    (wd / "prof").mkdir(parents=True, exist_ok=True)
    (wd / "norm").mkdir(parents=True, exist_ok=True)
    lengths = {name: len(seq) for name, seq in fasta}

    # 1. Count mod and the untreated control over the same references.
    mod_counts = run_shapemapper_counts(
        mod_bam,
        fasta,
        sm_dir,
        params,
        reference=reference,
        samtools=st,
        min_depth=min_depth,
        workdir=wd,
        prefix="mod",
    )
    nomod_counts = (
        run_shapemapper_counts(
            nomod_bam,
            fasta,
            sm_dir,
            params,
            references=set(mod_counts),
            reference=reference,
            samtools=st,
            min_depth=min_depth,
            workdir=wd,
            prefix="nomod",
        )
        if nomod_bam is not None
        else {}
    )

    # 2. Per reference: build a reactivity profile (background-subtract).
    order: list[str] = []
    for name, seq in fasta:
        if name not in mod_counts:
            continue
        count_files = [str(mod_counts[name])]
        if name in nomod_counts:
            count_files.append(str(nomod_counts[name]))
        ref_fa = wd / "ref.fa"
        ref_fa.write_text(f">{name}\n{seq}\n")
        idx = len(order)
        prof_path = wd / "prof" / f"{idx}.txt"
        cmd = [
            py,
            make_prof,
            "--fa",
            str(ref_fa),
            "--counts",
            *count_files,
            "--out",
            str(prof_path),
            "--mindepth",
            str(min_depth),
        ]
        if dms:
            cmd.append("--dms")
        subprocess.run(cmd, env=env, capture_output=True)
        if prof_path.exists():
            order.append(name)

    # 3. Normalize all profiles together (one library-wide factor).
    tonorm = [str(wd / "prof" / f"{i}.txt") for i in range(len(order))]
    normout = [str(wd / "norm" / f"{i}.txt") for i in range(len(order))]
    if tonorm:
        ncmd = [py, norm_prof, "--tonorm", *tonorm, "--normout", *normout, "--warn-on-error"]
        if dms:
            ncmd.append("--dms")
        subprocess.run(ncmd, env=env, capture_output=True)

    # 4. Read each reference's profile (normalized if produced, else raw), and
    #    its per-reference depth (max effective_depth of the modified counts).
    out: dict[str, np.ndarray] = {}
    for i, name in enumerate(order):
        normed = Path(normout[i])
        src = normed if normed.exists() else Path(tonorm[i])
        out[name] = parse_shapemapper_profile(src, lengths[name])
    reads = {
        name: _max_combined(
            [_counts_depth(path)]
            + ([_counts_depth(nomod_counts[name])] if name in nomod_counts else [])
        )
        for name, path in mod_counts.items()
    }
    return out, reads


# ---------------------------------------------------------------------------
# Unified interface
# ---------------------------------------------------------------------------

TOOLS = ("cmuts", "rnaframework", "shapemapper2")


def _ensure_index(path, samtools: str) -> None:
    """Index a coordinate-sorted BAM/CRAM if its companion index is missing."""
    path = Path(path)
    suffix = ".crai" if path.suffix.lower() == ".cram" else ".bai"
    if not Path(str(path) + suffix).exists():
        run_checked([samtools, "index", str(path)])


def _prepare_input(tool: str, path, fasta, workdir, *, samtools: str | None = None) -> Path:
    """Adapt one alignment to the format `tool` consumes, converting in `workdir`.

    cmuts reads bam/cram/sam natively (passthrough); rf-count needs a BAM;
    shapemapper2 needs an indexed bam/cram (a SAM is sorted+indexed to BAM). This
    is the single home for the per-tool format handling, so callers just pass
    whatever format they have.
    """
    st = samtools or samtools_bin()
    path = Path(path)
    ext = path.suffix.lower()
    if tool == "cmuts":
        return path
    if tool == "rnaframework":
        if ext == ".bam":
            return path
        out = Path(workdir) / f"{path.stem}.bam"
        run_checked([st, "view", "-b", "-T", str(fasta), "-o", str(out), str(path)])
        return out
    if tool == "shapemapper2":
        if ext in (".bam", ".cram"):
            _ensure_index(path, st)
            return path
        out = Path(workdir) / f"{path.stem}.bam"
        run_checked([st, "sort", "-o", str(out), str(path)])
        run_checked([st, "index", str(out)])
        return out
    raise ValueError(f"unknown tool: {tool}")


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

    The single entry point: adapt each input to the format the tool needs
    (conversion is part of the work), then count + normalize. `condition` (DMS /
    2A3) names the experiment and selects condition-dependent knobs (shapemapper2
    DMS mode; cmuts `cmuts_norm == "sm"` -> sm-dms/sm-shape). Writes a unified
    HDF5 (`/names`, `/reactivity` (n_ref, length), `/reads` (n_ref,) per-reference
    depth) and returns its path. It applies no coverage floor of its own --
    downstream filtering on `/reads` is the caller's choice.
    """
    import h5py
    import numpy as np

    st = samtools_bin()
    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="react-"))
    wd.mkdir(parents=True, exist_ok=True)
    fasta_list = read_fasta(inputs.fasta)
    names = [n for n, _ in fasta_list]
    length = max((len(s) for _, s in fasta_list), default=0)
    mod = _prepare_input(tool, inputs.mod, inputs.fasta, wd, samtools=st)
    nomod = (
        _prepare_input(tool, inputs.nomod, inputs.fasta, wd, samtools=st)
        if inputs.nomod is not None
        else None
    )
    if tool == "cmuts":
        react_d, reads_d = _cmuts_reactivity(
            count, norm, inputs.fasta, mod, nomod, condition, names, wd
        )
    elif tool == "rnaframework":
        react_d, reads_d = _rf_reactivity(count, norm, inputs.fasta, mod, nomod, wd)
    elif tool == "shapemapper2":
        if sm_dir is None:
            raise ValueError("shapemapper2 requires sm_dir")
        react_d, reads_d = run_shapemapper_reactivity(
            mod,
            nomod,
            fasta_list,
            sm_dir,
            count,
            dms=(condition == "DMS"),
            reference=str(inputs.fasta),
            samtools=st,
            min_depth=_FLOOR,
            workdir=wd,
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


def _cmuts_reactivity(count, norm, fasta, mod, nomod, condition, names, wd):
    bams = [mod] + ([nomod] if nomod is not None else [])
    counts_h5 = Path(wd) / "cmuts-counts.h5"
    run_checked(cmuts_core_command(bams, fasta, counts_h5, count))
    method = norm.cmuts_norm
    if method == "sm":
        method = "sm-dms" if condition == "DMS" else "sm-shape"
    react, reads = run_cmuts_normalize(
        counts_h5,
        mod,
        nomod,
        fasta,
        condition,
        method,
        norm.cmuts_per_reference,
        Path(wd) / "cmuts-norm.h5",
    )
    n = min(len(names), react.shape[0])
    return (
        {names[i]: react[i] for i in range(n)},
        {names[i]: float(reads[i]) for i in range(min(len(names), len(reads)))},
    )


def _rf_reactivity(count, norm, fasta, mod, nomod, wd):
    rc_mod = run_rfcount(mod, fasta, Path(wd) / "rf-mod", count)
    rc_nomod = (
        run_rfcount(nomod, fasta, Path(wd) / "rf-nomod", count) if nomod is not None else None
    )
    norm_dir = run_rfnorm(
        rc_mod,
        Path(wd) / "rf-norm",
        scoring_method=norm.rf_scoring,
        norm_method=norm.rf_norm_method,
        untreated=rc_nomod,
    )
    cov_mod = rfcount_coverage(rc_mod)
    cov_nomod = rfcount_coverage(rc_nomod) if rc_nomod is not None else {}
    reads = {
        name: _max_combined([c] + ([cov_nomod[name]] if name in cov_nomod else []))
        for name, c in cov_mod.items()
    }
    return parse_rfnorm(norm_dir), reads


def to_array(per_ref: dict[str, np.ndarray], names: list[str], length: int) -> np.ndarray:
    """Stack {ref_name: 1-D reactivity} into a (len(names), length) array in
    `names` order; missing references/positions are NaN."""
    import numpy as np

    out = np.full((len(names), length), np.nan)
    for i, name in enumerate(names):
        v = per_ref.get(name)
        if v is not None:
            out[i, : len(v)] = v[:length]
    return out


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


# ---------------------------------------------------------------------------
# Command-line entry point (lets profile.py time `reactivity` under GNU time)
# ---------------------------------------------------------------------------


def _main(argv: list[str]) -> int:
    """Run `reactivity` from the command line so the profile benchmark can time a
    tool's whole pipeline in a subprocess. Output is discarded (timing only)."""
    import argparse
    import json

    p = argparse.ArgumentParser(prog="external.py")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("reactivity", help="run a tool's full reactivity pipeline")
    r.add_argument("--tool", required=True, choices=TOOLS)
    r.add_argument("--count", required=True, help="JSON file of CountParams fields")
    r.add_argument("--norm", required=True, help="JSON file of NormParams fields")
    r.add_argument("--fasta", required=True)
    r.add_argument("--mod", required=True)
    r.add_argument("--nomod", default=None)
    r.add_argument("--condition", default="DMS")
    r.add_argument("--sm-dir", default=None)
    r.add_argument("--workdir", required=True)
    r.add_argument("--out", required=True, help="output profiles HDF5")
    args = p.parse_args(argv)

    if args.cmd == "reactivity":
        with open(args.count) as fh:
            count = CountParams(**json.load(fh))
        with open(args.norm) as fh:
            norm = NormParams(**json.load(fh))
        reactivity(
            args.tool,
            count,
            norm,
            Inputs(Path(args.fasta), Path(args.mod), Path(args.nomod) if args.nomod else None),
            args.condition,
            args.out,
            sm_dir=args.sm_dir,
            workdir=args.workdir,
        )
    return 0


if __name__ == "__main__":
    import sys

    raise SystemExit(_main(sys.argv[1:]))
