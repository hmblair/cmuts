#!/usr/bin/env python3
"""Benchmark the full reactivity pipeline of cmuts, rf-count, and shapemapper2.

Measures wall-clock time and peak memory for each tool's whole pipeline (count +
normalize, via `external.reactivity`), across three sweeps: number of queries,
number of references, and reference length. Timing the full pipeline -- not just
counting -- reflects what users actually run.

Each tool is timed on each requested input format (--format, default bam).
external.py adapts each format to what the tool needs (cmuts reads all natively;
rf-count and shapemapper2 convert as required), so any conversion cost is part of
the measured pipeline. Synthetic modified + untreated alignments are generated
with `cmuts generate` and cached under ./.cases; the untreated sample is a copy
of the modified one, so the reactivity values are degenerate but the work is
representative for timing.

Results are printed and appended to a CSV (--output). Plotting is left to the
figure sources, which read the CSV.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import external

# Quality/processing parameters, set to cmuts' defaults so cmuts runs out of the
# box; rf-count and shapemapper2 are configured to match (see _count_params), and
# the synthetic data is generated with the same thresholds.
MIN_MAPQ = 10
MIN_PHRED = 10
WINDOW = 2
MIN_LENGTH = 2
MAX_INDEL = 10
COLLAPSE = 2

FORMATS = ("bam", "cram", "sam")  # canonical order for parsing/dedup

EXP_H5 = "expected.h5"  # cmuts generate ground-truth counts (unused downstream)

SEP = " ─────────────────────────────"


def _gnu_time_bin() -> str:
    """Locate GNU time (gtime on macOS, /usr/bin/time on Linux)."""
    for cand in ("gtime", "/usr/bin/time"):
        path = cand if cand.startswith("/") else shutil.which(cand)
        if not path or not Path(path).exists():
            continue
        try:
            probe = subprocess.run(
                [path, "--format", "%e", "true"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except OSError:
            continue
        if probe.returncode == 0:
            return path
    sys.exit(
        "Error: GNU time not found. Install with 'brew install gnu-time' (macOS) "
        "or ensure GNU 'time' is installed (Linux)."
    )


TIME_BIN = ""  # set in main once samtools/time are confirmed available


def measure(cmd, runs, *, cwd=None, shell_script=None, pre=None):
    """Run a command `runs` times under GNU time; return (min_time_s, min_mem_mb).

    Time and memory minima are taken independently (as in the original harness).
    `pre` is an untimed callback run before each iteration (e.g. to clear output).
    `shell_script` runs a string via `bash -c` instead of an argv list.
    """
    best_t = None
    best_m = None
    for _ in range(runs):
        if pre is not None:
            pre()
        fd, tmp = tempfile.mkstemp()
        os.close(fd)
        try:
            inner = ["bash", "-c", shell_script] if shell_script is not None else list(cmd)
            proc = subprocess.run(
                [TIME_BIN, "--format", "%e %M", "-o", tmp, *inner],
                cwd=cwd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            raw = Path(tmp).read_text()
        finally:
            os.unlink(tmp)
        # GNU time exits with the wrapped command's status; on failure it writes a
        # diagnostic line instead of the "%e %M" numbers. Record NaN so one failed
        # point is visible in the CSV without aborting the whole sweep.
        try:
            if proc.returncode != 0:
                raise ValueError(f"exit {proc.returncode}")
            fields = raw.split()
            t = float(fields[0])
            mem_mb = float(fields[1]) / 1024
        except (ValueError, IndexError) as e:
            sys.stderr.write(f"  [warn] measurement failed ({e}): {' '.join(inner)[:120]}\n")
            return float("nan"), float("nan")
        if best_t is None or t < best_t:
            best_t = t
        if best_m is None or mem_mb < best_m:
            best_m = mem_mb
    return best_t, best_m


@dataclass
class Case:
    mode: str
    references: int
    queries: int
    length: int
    threads: int

    @property
    def base(self) -> str:
        return f"aln_{self.references}_{self.queries}_{self.length}"

    @property
    def fasta(self) -> str:
        return f"seq_{self.references}_{self.queries}_{self.length}.fasta"

    @property
    def max_length(self) -> int:
        return 2 * self.length


@dataclass
class Result:
    case: Case
    tool: str
    fmt: str
    time_s: float
    mem_mb: float


def _print_tool(label: str, time_s: float, mem_mb: float) -> None:
    print(f"   {label}:")
    print(f"    Execution time:    {time_s:.2f} seconds")
    print(f"    Peak memory usage: {mem_mb:.2f} MB")


def generate(case: Case, formats: list[str]) -> None:
    """Generate synthetic modified + untreated alignments (+ FASTA) for each
    format, unless cached. The untreated sample is a copy of the modified one
    (a distinct name so cmuts sees two samples); reactivity is degenerate but the
    pipeline work is representative for timing."""
    have = all(
        Path(f"{case.base}.{fmt}").exists() and Path(f"{case.base}.nomod.{fmt}").exists()
        for fmt in formats
    )
    if have:
        return
    fmt_flags = [f"--{fmt}" for fmt in formats]
    # Removing EXP_H5 first avoids cmuts generate's "dataset already exists" when
    # the cached ground-truth HDF5 is reused across cases; it is unused downstream.
    Path(EXP_H5).unlink(missing_ok=True)
    subprocess.run(
        [
            "cmuts",
            "generate",
            "--length",
            str(case.length),
            "--queries",
            str(case.queries),
            "--references",
            str(case.references),
            "--out-fasta",
            case.fasta,
            "-o",
            case.base,
            "--out-h5",
            EXP_H5,
            "--min-mapq",
            str(MIN_MAPQ),
            "--min-phred",
            str(MIN_PHRED),
            "--min-length",
            str(MIN_LENGTH),
            "--max-length",
            str(case.max_length),
            "--max-indel-length",
            str(MAX_INDEL),
            "--quality-window",
            str(WINDOW),
            "--collapse",
            str(COLLAPSE),
            *fmt_flags,
        ],
        check=False,
    )
    missing = [f"{case.base}.{fmt}" for fmt in formats if not Path(f"{case.base}.{fmt}").exists()]
    if missing:
        raise RuntimeError(f"cmuts generate did not produce {missing}")
    for fmt in formats:
        src, dst = f"{case.base}.{fmt}", f"{case.base}.nomod.{fmt}"
        shutil.copyfile(src, dst)
        for ext in (".bai", ".crai"):  # copy the index companion if cmuts generate wrote one
            if Path(src + ext).exists():
                shutil.copyfile(src + ext, dst + ext)


def _count_params(case: Case) -> external.CountParams:
    """Count-stage settings for every tool, mirroring the cmuts core flags the
    synthetic data was generated with so all three measure the same work."""
    return external.CountParams(
        insertions=True,
        collapse=COLLAPSE,
        cov_low_qual=False,
        threads=case.threads,
        min_mapq=MIN_MAPQ,
        min_phred=MIN_PHRED,
        min_length=MIN_LENGTH,
        max_length=case.max_length,
        max_indel=MAX_INDEL,
        quality_window=WINDOW,
        cmuts_spread="default",
    )


def _datasets(case: Case) -> list[tuple[str, external.CountParams, external.NormParams]]:
    """One timing dataset per tool: (tool, count params, norm params)."""
    count = _count_params(case)
    return [
        ("cmuts", count, external.NormParams(cmuts_norm="ubr")),
        ("rnaframework", count, external.NormParams(rf_scoring=3, rf_norm_method=1)),
        ("shapemapper2", count, external.NormParams()),
    ]


def _bench(
    dataset: tuple[str, external.CountParams, external.NormParams],
    fmt: str,
    case: Case,
    sm_dir: str,
    runs: int,
) -> Result | None:
    """Time one tool's full reactivity pipeline on one input format, via the
    `external.py reactivity` subprocess under GNU time. Returns None if the tool
    is unavailable. `--min-depth 1` counts every reference with reads."""
    tool, count, norm = dataset
    if tool == "rnaframework" and not external.rfcount_available():
        return None
    if tool == "shapemapper2" and not sm_dir:
        return None

    work, count_json, norm_json, out_h5 = "_react_work", "_count.json", "_norm.json", "_react.h5"
    with open(count_json, "w") as fh:
        json.dump(dataclasses.asdict(count), fh)
    with open(norm_json, "w") as fh:
        json.dump(dataclasses.asdict(norm), fh)

    cmd = [
        sys.executable,
        external.__file__,
        "reactivity",
        "--tool",
        tool,
        "--count",
        count_json,
        "--norm",
        norm_json,
        "--fasta",
        case.fasta,
        "--mod",
        f"{case.base}.{fmt}",
        "--nomod",
        f"{case.base}.nomod.{fmt}",
        "--condition",
        "DMS",
        "--workdir",
        work,
        "--out",
        out_h5,
    ]
    if sm_dir:
        cmd += ["--sm-dir", sm_dir]

    def reset() -> None:
        shutil.rmtree(work, ignore_errors=True)
        Path(out_h5).unlink(missing_ok=True)

    t, m = measure(cmd, runs, pre=reset)
    shutil.rmtree(work, ignore_errors=True)
    for f in (count_json, norm_json, out_h5):
        Path(f).unlink(missing_ok=True)
    _print_tool(f"{tool} [{fmt}]", t, m)
    return Result(case, tool, fmt, t, m)


def profile(case: Case, formats: list[str], runs: int, sm_dir: str) -> list[Result]:
    """Time every (tool, format) for one case."""
    results: list[Result] = []
    for dataset in _datasets(case):
        for fmt in formats:
            r = _bench(dataset, fmt, case, sm_dir, runs)
            if r is not None:
                results.append(r)
    return results


# Sweep definitions: (mode, swept label, values, fixed-header lines, case factory).
# The case factory takes (swept value, thread count).
def _sweeps():
    return [
        (
            "queries",
            "QUERIES",
            [10**2, 10**4, 10**6],
            [("REFERENCES", 1024), ("REFERENCE LENGTH", 200)],
            lambda v, t: Case("queries", 1024, v, 200, t),
        ),
        (
            "references",
            "REFERENCES",
            [10**1, 10**2, 10**3, 10**4, 10**5],
            [("QUERIES", 1024), ("REFERENCE LENGTH", 200)],
            lambda v, t: Case("references", v, 1024, 200, t),
        ),
        (
            "lengths",
            "REFERENCE LENGTH",
            [10**1, 10**2, 10**3, 10**4],
            [("QUERIES", 1024), ("REFERENCES", 1024)],
            lambda v, t: Case("lengths", 1024, 1024, v, t),
        ),
    ]


def _run_sweep(spec, args, output: Path, all_rows: list[Result]) -> None:
    mode, var_label, values, header, make_case = spec
    cap = {
        "queries": args.max_queries,
        "references": args.max_references,
        "lengths": args.max_length,
    }[mode]
    if cap is not None:
        values = [v for v in values if v <= cap]
        if not values:
            sys.exit(f"No {mode} values <= {cap}; lower the cap or widen the sweep.")
    print()
    print(SEP)
    for label, val in header:
        print(f"   {label + ':':<18}{val}")
    print(SEP)
    for v in values:
        print(f"   {var_label + ':':<18}{v}")
        print(SEP)
        # Data is thread-independent, so generate once and reuse for each count.
        generate(make_case(v, args.threads[0]), args.formats)
        for t in args.threads:
            print(f"   {'THREADS:':<18}{t}")
            case = make_case(v, t)
            rows = profile(case, args.formats, args.runs, args.shapemapper_dir)
            _append_csv(output, rows, args.runs)
            all_rows.extend(rows)
            print(SEP)


def _append_csv(path: Path, rows: list[Result], runs: int) -> None:
    """Append rows to the CSV, writing the header if the file is new/empty.

    Called per case as results come in so a timeout or crash mid-sweep keeps
    everything measured so far. The handle is flushed and fsync'd before return.
    """
    new = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(
                [
                    "mode",
                    "tool",
                    "format",
                    "threads",
                    "references",
                    "queries",
                    "length",
                    "runs",
                    "time_s",
                    "mem_mb",
                ]
            )
        for r in rows:
            w.writerow(
                [
                    r.case.mode,
                    r.tool,
                    r.fmt,
                    r.case.threads,
                    r.case.references,
                    r.case.queries,
                    r.case.length,
                    runs,
                    f"{r.time_s:.4f}",
                    f"{r.mem_mb:.2f}",
                ]
            )
        fh.flush()
        os.fsync(fh.fileno())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs", type=int, default=1, help="Iterations per measurement (min is kept)")
    p.add_argument(
        "--threads",
        type=int,
        nargs="+",
        default=[1],
        metavar="N",
        help="Thread counts for cmuts/rf-count; multiple allowed (e.g. --threads 1 8)",
    )
    p.add_argument("--queries", action="store_true", help="Sweep number of queries")
    p.add_argument("--references", action="store_true", help="Sweep number of references")
    p.add_argument("--lengths", action="store_true", help="Sweep reference length")
    p.add_argument(
        "--max-queries",
        type=int,
        default=None,
        metavar="N",
        help="Cap the queries sweep to values <= N",
    )
    p.add_argument(
        "--max-references",
        type=int,
        default=None,
        metavar="N",
        help="Cap the references sweep to values <= N",
    )
    p.add_argument(
        "--max-length",
        type=int,
        default=None,
        metavar="N",
        help="Cap the reference-length sweep to values <= N",
    )
    p.add_argument(
        "--formats",
        default="bam",
        metavar="LIST",
        help="Comma-separated input formats to benchmark (bam,cram,sam); default: bam",
    )
    p.add_argument("--shapemapper-dir", default="", metavar="DIR", help="ShapeMapper2 install dir")
    p.add_argument(
        "--output",
        default="outputs/profile-results.csv",
        help="CSV results path (appended); parent dir is created. Default is the gitignored outputs/.",
    )
    args = p.parse_args()

    if not shutil.which("samtools"):
        sys.exit("Samtools is required but it's not installed. Exiting.")
    if args.shapemapper_dir:
        for tool in ("shapemapper_mutation_parser", "shapemapper_mutation_counter"):
            if not os.access(f"{args.shapemapper_dir}/internals/bin/{tool}", os.X_OK):
                sys.exit(f"{tool} not found in {args.shapemapper_dir}/internals/bin/. Exiting.")

    requested = [f.strip().lower() for f in args.formats.split(",") if f.strip()]
    bad = [f for f in requested if f not in FORMATS]
    if bad:
        sys.exit(f"Unknown format(s): {', '.join(bad)}. Choose from {', '.join(FORMATS)}.")
    args.formats = [f for f in FORMATS if f in requested]  # canonical order, deduped
    if not args.formats:
        sys.exit("Nothing to do: --formats is empty.")

    global TIME_BIN
    TIME_BIN = _gnu_time_bin()

    enabled = {"queries": args.queries, "references": args.references, "lengths": args.lengths}
    selected = [s for s in _sweeps() if enabled[s[0]]]
    if not selected:
        sys.exit("Nothing to do: pass at least one of --queries, --references, --lengths.")

    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    cases_dir = Path.cwd() / ".cases"
    cases_dir.mkdir(exist_ok=True)

    all_rows: list[Result] = []
    prev = Path.cwd()
    os.chdir(cases_dir)
    try:
        for spec in selected:
            _run_sweep(spec, args, output, all_rows)
    finally:
        os.chdir(prev)

    # Rows are appended per case in _run_sweep; nothing left to flush here.
    print()
    print(f"Wrote {len(all_rows)} rows to {output}")


if __name__ == "__main__":
    main()
