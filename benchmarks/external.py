"""Shared invocation and output parsing for the external comparison tools.

cmuts is benchmarked against rf-count (RNAFramework) and shapemapper2. Those two
read different formats and take different flags; this module is the single place
that knows how to *call* each tool and *parse* its output, so the profile,
correctness, and accuracy benchmarks do not each re-implement it.

Quality/processing knobs live in `Params` and are mapped to each tool's own
flags here. Every benchmark passes the `Params` matching its intent:

    profile      mirrors the synthetic-data generation (counts insertions)
    correctness  matches cmuts defaults so the rate is the same quantity
    accuracy     tunes each tool for accuracy (its standard mutation set)

Tools are resolved from the environment so a machine without them still runs the
cmuts-only paths:

    RF_COUNT    rf-count       (default: "rf-count" on PATH)
    RF_RCTOOLS  rf-rctools     (default: "rf-rctools" on PATH)
    RF_NORM     rf-norm        (default: "rf-norm" on PATH)
    SAMTOOLS    samtools       (default: "samtools" on PATH)

shapemapper2 has no PATH convention upstream, so its install dir is passed in.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# A position needs at least this many reads for a rate to be meaningful.
MIN_DEPTH = 10


@dataclass
class Params:
    """Quality/processing knobs, mapped to each tool's flags by the builders below.

    The quality fields default to None: the corresponding flag is then omitted
    and the tool's own default applies. Fields with concrete defaults are always
    emitted, so a `Params()` is a sensible baseline (insertions counted, deletions
    right-aligned, consecutive mods collapsed within 1, low-quality coverage not
    counted).
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
    max_indel: int | None = None
    max_edit_distance: float | None = None  # rf-count only
    median_quality: int | None = None  # rf-count only
    max_internal_match: int | None = None  # shapemapper parser only


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


def samtools_view_mapped(
    src, out_sam, *, region: str | None = None, samtools: str | None = None
) -> None:
    """Write the mapped reads of `src` (BAM/CRAM/SAM) to `out_sam` as SAM, with
    header. Optionally restrict to one reference `region`. Untimed input prep
    that mirrors the unmapped-read filtering rf-count/cmuts apply internally.
    """
    cmd = [samtools or samtools_bin(), "view", "-h", "-F", "4", str(src)]
    if region:
        cmd.append(region)
    with open(out_sam, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


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


# ---------------------------------------------------------------------------
# rf-count (RNAFramework)
# ---------------------------------------------------------------------------


def rfcount_available() -> bool:
    return (
        shutil.which(os.environ.get("RF_COUNT", "rf-count")) is not None
        and shutil.which(os.environ.get("RF_RCTOOLS", "rf-rctools")) is not None
    )


def rfcount_command(bam, fasta, outdir, params: Params, *, overwrite: bool = True) -> list[str]:
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


def run_rfcount(bam, fasta, outdir, params: Params) -> Path:
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


def _sm_parser_flags(params: Params) -> list[str]:
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


def shapemapper_command(sam, parsed, counts, length: int, sm_dir: str, params: Params) -> str:
    """A bash one-liner (parser -> counter), with LD scoping, for timing in the
    profile benchmark. The same parser flags as the per-reference runner.
    """
    parser, counter = shapemapper_paths(sm_dir)
    libdir = _shapemapper_libdir(sm_dir)
    ld = (
        f'export LD_LIBRARY_PATH="{libdir}:${{LD_LIBRARY_PATH:-}}"; '
        if Path(libdir).is_dir()
        else ""
    )
    pflags = " ".join(_sm_parser_flags(params))
    return (
        f"{ld}{parser} -i {sam} -o {parsed} {pflags} && "
        f"{counter} -i {parsed} -c {counts} -n {length} --warn_on_no_mapped"
    )


def parse_shapemapper_counts(
    path, length: int, params: Params, *, min_depth: int = MIN_DEPTH
) -> np.ndarray:
    """Per-position mutation rate from a shapemapper2 counts file.

    rate = (selected mutation categories) / effective_depth. Insertion columns
    are excluded when params.insertions is False (to match cmuts --no-insertions).
    """
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


def _count_reference(
    bam,
    name: str,
    length: int,
    *,
    parser: str,
    counter: str,
    params: Params,
    env: dict[str, str],
    samtools: str,
    workdir: Path,
    prefix: str = "sm",
) -> Path | None:
    """Extract one reference's reads and run the shapemapper parser then counter.

    Returns the per-position counts file, or None if extraction/counting failed.
    """
    sam = workdir / f"{prefix}.sam"
    parsed = workdir / f"{prefix}.mut"
    counts = workdir / f"{prefix}-counts.txt"
    with open(sam, "w") as fh:
        view = subprocess.run([samtools, "view", "-h", "-F", "4", str(bam), name], stdout=fh)
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


def run_shapemapper_rates(
    bam,
    fasta: list[tuple[str, str]],
    sm_dir: str,
    params: Params,
    *,
    samtools: str | None = None,
    min_depth: int = MIN_DEPTH,
    workdir=None,
) -> dict[str, np.ndarray]:
    """Per-position raw mutation rate per reference from shapemapper2.

    The shapemapper2 counter is single-reference (length via -n), so each
    reference with coverage is processed independently: extract its reads, run
    the parser then the counter, and read the per-position counts. References
    below `min_depth` total reads are skipped (idxstats pre-filter), so the loop
    is bounded by covered references rather than the whole FASTA. Requires an
    indexed BAM. (For the normalized reactivity profile, see
    `run_shapemapper_reactivity`.)
    """
    parser, counter = shapemapper_paths(sm_dir)
    st = samtools or samtools_bin()
    covered = covered_references(bam, min_reads=min_depth, samtools=st)

    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="sm-"))
    wd.mkdir(parents=True, exist_ok=True)
    count = functools.partial(
        _count_reference,
        parser=parser,
        counter=counter,
        params=params,
        env=shapemapper_env(sm_dir),
        samtools=st,
        workdir=wd,
    )

    rates: dict[str, np.ndarray] = {}
    for name, seq in fasta:
        if name not in covered:
            continue
        counts = count(bam, name, len(seq))
        if counts is not None:
            rates[name] = parse_shapemapper_counts(counts, len(seq), params, min_depth=min_depth)
    return rates


def parse_shapemapper_profile(path, length: int) -> np.ndarray:
    """Read the reactivity column from a shapemapper profile / normalized file.

    Prefers `Norm_profile` (normalized) and falls back to `Reactivity_profile`
    (background-subtracted but un-normalized). Non-finite entries become NaN.
    """
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
    params: Params,
    *,
    dms: bool = False,
    samtools: str | None = None,
    min_depth: int = MIN_DEPTH,
    workdir=None,
) -> dict[str, np.ndarray]:
    """shapemapper2 native reactivity per reference (mod + untreated control).

    Per reference: parser+counter on the modified and (if given) untreated
    samples, then `make_reactivity_profiles.py` to background-subtract. All
    per-reference profiles are then normalized together in a single
    `normalize_profiles.py` call -- its factor is computed across the whole set,
    since one short reference is too sparse to normalize alone. Returns
    {ref_name: normalized reactivity}, falling back per reference to the
    un-normalized Reactivity_profile where normalization could not be applied.
    """
    parser, counter = shapemapper_paths(sm_dir)
    env = shapemapper_env(sm_dir)
    py = shapemapper_python(sm_dir)
    make_prof = f"{sm_dir}/internals/bin/make_reactivity_profiles.py"
    norm_prof = f"{sm_dir}/internals/bin/normalize_profiles.py"
    st = samtools or samtools_bin()
    if not os.path.exists(make_prof):
        return {}
    covered = covered_references(mod_bam, min_reads=min_depth, samtools=st)

    wd = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="sm-react-"))
    (wd / "prof").mkdir(parents=True, exist_ok=True)
    (wd / "norm").mkdir(parents=True, exist_ok=True)
    lengths = {name: len(seq) for name, seq in fasta}
    count = functools.partial(
        _count_reference,
        parser=parser,
        counter=counter,
        params=params,
        env=env,
        samtools=st,
        workdir=wd,
    )

    # 1. Per reference: count mod (+ untreated) and build a reactivity profile.
    order: list[str] = []
    for name, seq in fasta:
        if name not in covered:
            continue
        mod_counts = count(mod_bam, name, len(seq), prefix="mod")
        if mod_counts is None:
            continue
        count_files = [str(mod_counts)]
        if nomod_bam is not None:
            nomod_counts = count(nomod_bam, name, len(seq), prefix="nomod")
            if nomod_counts is not None:
                count_files.append(str(nomod_counts))
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

    # 2. Normalize all profiles together (one library-wide factor).
    tonorm = [str(wd / "prof" / f"{i}.txt") for i in range(len(order))]
    normout = [str(wd / "norm" / f"{i}.txt") for i in range(len(order))]
    if tonorm:
        ncmd = [py, norm_prof, "--tonorm", *tonorm, "--normout", *normout, "--warn-on-error"]
        if dms:
            ncmd.append("--dms")
        subprocess.run(ncmd, env=env, capture_output=True)

    # 3. Read each reference's profile (normalized if produced, else raw).
    out: dict[str, np.ndarray] = {}
    for i, name in enumerate(order):
        normed = Path(normout[i])
        src = normed if normed.exists() else Path(tonorm[i])
        out[name] = parse_shapemapper_profile(src, lengths[name])
    return out
