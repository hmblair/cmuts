#!/usr/bin/env python3
"""Per-nucleotide correctness comparison: cmuts vs rf-count and shapemapper2.

Runs all three frameworks on one BAM against one FASTA, extracts the per-position
raw mutation rate for every reference, and writes the absolute per-position
difference between cmuts and each of the other two tools. cmuts is run with
deletion spreading OFF (--no-spread) and insertions off (--no-insertions); the
other two are configured to match (see `MATCH`), so all three compute the same
quantity. This is the Fig 2a control that isolates deletion spreading, not an
accuracy result.

The external tools (rf-count, shapemapper2) and their output formats are handled
by the shared `external` module; this script only wires cmuts to them and writes
the comparison. A tool that is absent is skipped, so this develops locally with
cmuts alone and runs fully where all three are installed. See `external` for the
environment variables (RF_COUNT, RF_RCTOOLS, SM2_DIR, SAMTOOLS).

Quantity compared: raw mutation rate = mutations / read depth at each position,
NaN where depth is below MIN_DEPTH. Reference identity comes from the FASTA
(cmuts indexes by FASTA order; the others key by reference name).

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
import tempfile
from pathlib import Path

import external
import h5py
import numpy as np
from external import MIN_DEPTH

# rf-count and shapemapper2 settings that match cmuts' defaults (min-mapq 10,
# min-phred 10, min-length 2, max-indel 10, collapse 2) with insertions off and
# deletions right-aligned, so the three tools count the same events.
MATCH = external.Params(
    insertions=False,
    right_align_deletions=True,
    collapse=2,
    min_mapq=10,
    min_phred=10,
    min_length=2,
    max_indel=10,
    max_internal_match=2,
)


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Ordered (name, sequence) list. Names are the header up to the first space."""
    out: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def cmuts_rates(bam: str, fasta: str, names: list[str], workdir: Path) -> dict[str, np.ndarray]:
    """Raw per-position mutation rate per reference, deletion spreading off.

    cmuts core counts with --no-spread --no-insertions (its default min-mapq/phred
    of 10 match `MATCH`); cmuts normalize --norm raw makes the reactivity the raw
    mutation rate.
    """
    cmuts = os.environ.get("CMUTS", "cmuts")
    counts = workdir / "cmuts-counts.h5"
    rate_h5 = workdir / "cmuts-rate.h5"
    external.run_checked(
        [
            cmuts,
            "core",
            "-f",
            fasta,
            "-o",
            str(counts),
            "--overwrite",
            "--no-spread",
            "--no-insertions",
            bam,
        ]
    )
    mod = bam[:-4] if bam.endswith(".bam") else bam
    external.run_checked(
        [
            cmuts,
            "normalize",
            "-o",
            str(rate_h5),
            "--fasta",
            fasta,
            "--mod",
            mod,
            "--group",
            "rate",
            "--norm",
            "raw",
            "--overwrite",
            str(counts),
        ]
    )
    with h5py.File(rate_h5, "r") as f:
        reactivity = np.asarray(f["rate/reactivity"])  # (n_refs, length), FASTA order
    return {name: reactivity[i] for i, name in enumerate(names)}


def rf_count_rates(bam: str, fasta: str, workdir: Path) -> dict[str, np.ndarray]:
    """Raw per-position mutation rate per reference from rf-count."""
    if not external.rfcount_available():
        print("  rf-count / rf-rctools not found; skipping RNAFramework", file=sys.stderr)
        return {}
    rc = external.run_rfcount(bam, fasta, workdir / "rf-count", MATCH)
    return external.parse_rfcount(rc)


def shapemapper_rates(
    bam: str, fasta: list[tuple[str, str]], workdir: Path
) -> dict[str, np.ndarray]:
    """Raw per-position mutation rate per reference from shapemapper2."""
    sm_dir = os.environ.get("SM2_DIR")
    if not sm_dir:
        print("  SM2_DIR unset; skipping shapemapper2", file=sys.stderr)
        return {}
    return external.run_shapemapper_rates(bam, fasta, sm_dir, MATCH, workdir=workdir / "sm")


def write_csv(
    path: str,
    names: list[str],
    cmuts: dict[str, np.ndarray],
    rf: dict[str, np.ndarray],
    sm: dict[str, np.ndarray],
) -> int:
    """Write per (reference, position) rates and |cmuts - other| differences."""
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "reference",
                "position",
                "cmuts",
                "rf_count",
                "shapemapper2",
                "abs_diff_rf",
                "abs_diff_sm",
            ]
        )
        for name in names:
            c = cmuts.get(name)
            if c is None:
                continue
            r = rf.get(name)
            s = sm.get(name)
            for pos in range(len(c)):
                cv = c[pos]
                if np.isnan(cv):
                    continue
                rv = r[pos] if r is not None and pos < len(r) else np.nan
                sv = s[pos] if s is not None and pos < len(s) else np.nan
                w.writerow(
                    [
                        name,
                        pos,
                        f"{cv:.6g}",
                        "" if np.isnan(rv) else f"{rv:.6g}",
                        "" if np.isnan(sv) else f"{sv:.6g}",
                        "" if np.isnan(rv) else f"{abs(cv - rv):.6g}",
                        "" if np.isnan(sv) else f"{abs(cv - sv):.6g}",
                    ]
                )
                n += 1
    return n


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bam", help="Aligned, indexed BAM file")
    p.add_argument("-f", "--fasta", required=True, help="Reference FASTA")
    p.add_argument("-o", "--out", default="correctness.csv", help="Output CSV")
    p.add_argument("--keep", action="store_true", help="Keep the work directory")
    args = p.parse_args()

    fasta = read_fasta(args.fasta)
    names = [n for n, _ in fasta]
    print(f"{len(names)} references (MIN_DEPTH={MIN_DEPTH}); BAM {args.bam}")

    workdir = Path(tempfile.mkdtemp(prefix="correctness-"))
    try:
        print("running cmuts (no-spread) ...")
        cmuts = cmuts_rates(args.bam, args.fasta, names, workdir)
        print(f"  cmuts: {sum(np.any(~np.isnan(v)) for v in cmuts.values())} refs with data")

        print("running rf-count ...")
        rf = rf_count_rates(args.bam, args.fasta, workdir)
        print(f"  rf-count: {len(rf)} refs")

        print("running shapemapper2 ...")
        sm = shapemapper_rates(args.bam, fasta, workdir)
        print(f"  shapemapper2: {len(sm)} refs")

        rows = write_csv(args.out, names, cmuts, rf, sm)
        print(f"wrote {rows} positions to {args.out}")
    finally:
        if not args.keep:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            print(f"work dir: {workdir}")


if __name__ == "__main__":
    main()
