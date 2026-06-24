#!/usr/bin/env python3
"""Benchmark cmuts core against rf-count and shapemapper2.

Measures wall-clock time and peak memory for cmuts core, rf-count
(RNAFramework), and the shapemapper2 mutation parser/counter, across three
sweeps: number of queries, number of references, and reference length.

cmuts reads BAM, CRAM, and SAM natively; the formats to benchmark are chosen
with --formats (comma-separated, default bam). rf-count reads only BAM and
shapemapper2 only SAM, so a requested non-native format adds a one-off samtools
conversion to their native runtime (timed separately and added); the tools
themselves are never re-run.

With --convert, each non-BAM format is also reported as a "cmuts-via-bam" row:
the samtools cost of converting it to BAM plus cmuts' BAM runtime, to compare
pre-converting against reading the format directly.

Synthetic data is generated with `cmuts generate` and cached under ./.cases.
Results are printed and appended to a CSV (--output). Plotting is left to the
figure sources, which can read the CSV.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import external

# Quality/processing parameters, set to cmuts' defaults so cmuts runs out of the
# box; rf-count and shapemapper2 are configured to match (see _profile_params),
# and the synthetic data is generated with the same thresholds.
MIN_MAPQ = 10
MIN_PHRED = 10
WINDOW = 2
MIN_LENGTH = 2
MAX_INDEL = 10
COLLAPSE = 2

FORMATS = ("bam", "cram", "sam")  # canonical order for parsing/dedup

EXP_H5 = "expected.h5"  # cmuts generate ground-truth counts (unused downstream)
CMUTS_H5 = "cmuts.h5"  # cmuts core output, removed after each run
RF_DIR = "_rf_count"  # rf-count writes here; recreated each run

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
    mem_mb: float | None  # None for convert+run sums (memory is not additive)


def _print_tool(label: str, time_s: float, mem_mb: float | None = None) -> None:
    print(f"   {label}:")
    print(f"    Execution time:    {time_s:.2f} seconds")
    if mem_mb is not None:
        print(f"    Peak memory usage: {mem_mb:.2f} MB")


def generate(case: Case, formats: list[str]) -> None:
    """Generate the requested synthetic alignments + reference FASTA, unless
    already cached. All formats are produced in one invocation so the reads are
    identical across them, so regenerate only when some requested file is missing."""
    if all(Path(f"{case.base}.{fmt}").exists() for fmt in formats):
        return
    fmt_flags = [f"--{fmt}" for fmt in formats]
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
        check=True,
    )


def _bench_cmuts(case: Case, fmt: str, runs: int) -> Result:
    cmd = [
        "cmuts",
        "core",
        "--fasta",
        case.fasta,
        "--output",
        CMUTS_H5,
        "--min-mapq",
        str(MIN_MAPQ),
        "--min-phred",
        str(MIN_PHRED),
        "--quality-window",
        str(WINDOW),
        "--min-length",
        str(MIN_LENGTH),
        "--max-length",
        str(case.max_length),
        "--max-indel-length",
        str(MAX_INDEL),
        "--threads",
        str(case.threads),
        "--overwrite",
        f"{case.base}.{fmt}",
    ]
    t, m = measure(cmd, runs, pre=lambda: shutil.rmtree(CMUTS_H5, ignore_errors=True))
    shutil.rmtree(CMUTS_H5, ignore_errors=True)
    _print_tool(f"cmuts [{fmt}]", t, m)
    return Result(case, "cmuts", fmt, t, m)


def _profile_params(case: Case) -> external.Params:
    """Quality/processing settings applied to every external tool, matching the
    cmuts core flags above so all three measure the same work."""
    return external.Params(
        insertions=True,
        right_align_deletions=True,
        collapse=COLLAPSE,
        cov_low_qual=False,
        discard_duplicates=False,
        fast=True,
        threads=case.threads,
        min_mapq=MIN_MAPQ,
        min_phred=MIN_PHRED,
        min_length=MIN_LENGTH,
        max_indel=MAX_INDEL,
        max_edit_distance=1.0,
        median_quality=0,
    )


def _bench_rf_count(case: Case, runs: int) -> tuple[float, float]:
    """Time rf-count on its native BAM input; returns (time_s, mem_mb)."""

    def reset() -> None:
        shutil.rmtree(RF_DIR, ignore_errors=True)
        os.mkdir(RF_DIR)

    cmd = external.rfcount_command(f"{case.base}.bam", case.fasta, RF_DIR, _profile_params(case))
    t, m = measure(cmd, runs, pre=reset)
    _print_tool("rf-count [bam] (native run)", t, m)
    return t, m


def _bench_shapemapper(case: Case, sm_dir: str, runs: int) -> tuple[float, float]:
    """Time shapemapper2 on its native SAM input; returns (time_s, mem_mb)."""
    parsed, counts, sam_in = "_sm_parsed.mut", "_sm_counts.txt", "_sm_input.sam"

    # Drop unmapped reads up front (untimed), matching the input rf-count/cmuts see.
    external.samtools_view_mapped(f"{case.base}.sam", sam_in)

    script = external.shapemapper_command(
        sam_in, parsed, counts, case.length, sm_dir, _profile_params(case)
    )

    def reset() -> None:
        for f in (parsed, counts):
            Path(f).unlink(missing_ok=True)

    t, m = measure([], runs, shell_script=script, pre=reset)
    for f in (parsed, counts, sam_in):
        Path(f).unlink(missing_ok=True)
    _print_tool("shapemapper2 [sam] (native run)", t, m)
    return t, m


def _time_conversion(case: Case, src_fmt: str, dst_fmt: str, runs: int) -> tuple[float, float]:
    """Time a one-off samtools conversion of the case's src_fmt file to dst_fmt
    ('bam' or 'sam'); returns (time_s, mem_mb)."""
    src = f"{case.base}.{src_fmt}"
    if dst_fmt == "bam":
        out = "_conv.bam"
        cmd = ["samtools", "view", "-@", str(case.threads), "-T", case.fasta, "-b", "-o", out, src]
        t, m = measure(cmd, runs, pre=lambda: Path(out).unlink(missing_ok=True))
    else:  # sam: keep mapped reads only, matching shapemapper2's native input
        out = "_conv.sam"
        script = f"samtools view -h -F 4 -T {case.fasta} {src} > {out}"
        t, m = measure([], runs, shell_script=script, pre=lambda: Path(out).unlink(missing_ok=True))
    Path(out).unlink(missing_ok=True)
    return t, m


def _bench_external_formats(
    case: Case,
    formats: list[str],
    runs: int,
    rf: tuple[float, float] | None,
    sm: tuple[float, float] | None,
) -> list[Result]:
    """For BAM-only rf-count and SAM-only shapemapper2, report each requested
    format: the native run alone when the format matches the tool's input, else
    the native runtime plus a one-off samtools conversion. The tools are not re-run."""
    results: list[Result] = []
    if rf is None and sm is None:
        return results

    # A FASTA index is needed for CRAM decode and samtools -T.
    subprocess.run(
        ["samtools", "faidx", case.fasta],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    for tool, native, base in (("rf-count", "bam", rf), ("shapemapper2", "sam", sm)):
        if base is None:
            continue
        base_t, base_m = base
        for fmt in formats:
            if fmt == native:
                results.append(Result(case, tool, native, base_t, base_m))
                continue
            ct, cm = _time_conversion(case, fmt, native, runs)
            _print_tool(f"samtools {fmt}->{native}", ct, cm)
            _print_tool(f"{tool} [{fmt}] (= {fmt}->{native} + {native} run)", ct + base_t)
            results.append(Result(case, "samtools", f"{fmt}->{native}", ct, cm))
            results.append(Result(case, tool, fmt, ct + base_t, None))

    return results


def _bench_cmuts_via_bam(
    case: Case, formats: list[str], runs: int, cmuts_times: dict[str, float]
) -> list[Result]:
    """For each non-BAM format, measure the alternative of converting it to BAM
    with samtools and then running cmuts, to compare against reading the format
    directly. cmuts reads every format natively, so the question is whether the
    one-off conversion ever beats the (small) direct-read overhead. The cmuts BAM
    runtime is reused, not re-run; the summed convert+process time is recorded
    under the "cmuts-via-bam" tool so the CSV holds both sides of the comparison."""
    targets = [f for f in formats if f != "bam"]
    bam_time = cmuts_times.get("bam")
    if not targets or bam_time is None:
        return []

    results: list[Result] = []
    # A FASTA index is needed for CRAM decode and samtools -T.
    subprocess.run(
        ["samtools", "faidx", case.fasta],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    for fmt in targets:
        ct, cm = _time_conversion(case, fmt, "bam", runs)
        _print_tool(f"samtools {fmt}->bam", ct, cm)
        _print_tool(f"cmuts via bam [{fmt}] (= {fmt}->bam + bam run)", ct + bam_time)
        results.append(Result(case, "samtools", f"{fmt}->bam", ct, cm))
        results.append(Result(case, "cmuts-via-bam", fmt, ct + bam_time, None))

    return results


def profile(case: Case, formats: list[str], runs: int, sm_dir: str, convert: bool) -> list[Result]:
    results: list[Result] = []

    # cmuts reads BAM, CRAM, and SAM natively, so benchmark every requested format.
    cmuts_times: dict[str, float] = {}
    for fmt in formats:
        r = _bench_cmuts(case, fmt, runs)
        cmuts_times[fmt] = r.time_s
        results.append(r)

    # rf-count (native BAM) and shapemapper2 (native SAM) run once on their input.
    rf = _bench_rf_count(case, runs) if shutil.which("rf-count") else None
    sm = _bench_shapemapper(case, sm_dir, runs) if sm_dir else None
    results += _bench_external_formats(case, formats, runs, rf, sm)

    # Optionally compare direct reads against converting to BAM first.
    if convert:
        results += _bench_cmuts_via_bam(case, formats, runs, cmuts_times)

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
        generate(make_case(v, args.threads[0]), args.gen_formats)
        for t in args.threads:
            print(f"   {'THREADS:':<18}{t}")
            case = make_case(v, t)
            rows = profile(case, args.formats, args.runs, args.shapemapper_dir, args.convert)
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
                    "" if r.mem_mb is None else f"{r.mem_mb:.2f}",
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
    p.add_argument(
        "--convert",
        action="store_true",
        help="Also measure converting each non-BAM format to BAM (samtools) then running "
        "cmuts, to compare against reading it directly. Forces BAM into the benchmark.",
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

    # The --convert comparison needs cmuts' direct BAM runtime as the baseline, so
    # ensure BAM is benchmarked (and generated) even if the user didn't request it.
    if args.convert and "bam" not in args.formats:
        args.formats = [f for f in FORMATS if f in (*args.formats, "bam")]

    # rf-count and shapemapper2 run on their native BAM/SAM, so those must exist
    # even when not benchmarked directly (they back the conversion-cost rows).
    gen = set(args.formats)
    if shutil.which("rf-count"):
        gen.add("bam")
    if args.shapemapper_dir:
        gen.add("sam")
    args.gen_formats = [f for f in FORMATS if f in gen]

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
