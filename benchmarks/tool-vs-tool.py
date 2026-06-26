#!/usr/bin/env python3
"""Per-nucleotide agreement of cmuts (tuned to match each existing tool) vs that
tool, built from modified/untreated alignments via `external.reactivity`.

With deletion spreading off and its settings tuned to an existing pipeline, cmuts
should reproduce that pipeline's reactivity profile nucleotide-by-nucleotide --
the Fig 2a control that isolates deletion spreading (cmuts reproduces prior
tools), not a standalone accuracy result. Writes, per (comparison, condition,
reference, position), the two reactivities and their absolute difference.

Comparisons (cmuts-tuned-to-X vs X):
  rnaframework   cmuts-match-rf  vs rnaframework
  shapemapper2   cmuts-match-sm  vs shapemapper2   (skipped unless SM2_DIR is set)

Each tool's profile is built here from the alignments (no precomputed HDF5), so
this needs samtools, rf-count, and (for the shapemapper2 comparison) SM2_DIR.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import os
from pathlib import Path

import external
import numpy as np

CONDITIONS = ["DMS", "2A3"]

# One map-seq counting config shared by every tool (cmuts spread set per dataset):
# insertions off, deletions on and right-aligned, collapse 2, mapq/phred 10,
# eval-surrounding (cmuts window 1), duplicates kept.
MAP = external.CountParams(
    insertions=False,
    deletions=True,
    right_align_deletions=True,
    collapse=2,
    eval_surrounding=True,
    quality_window=1,
    cov_low_qual=True,
    discard_duplicates=False,
    min_mapq=10,
    min_phred=10,
    min_length=2,
    median_quality=0,
    max_edit_distance=1.0,
    max_internal_match=2,
)

# (cmuts dataset tuned to the tool, the tool dataset); labelled by the tool.
COMPARISONS = [
    ("cmuts-match-rf", "rnaframework"),
    ("cmuts-match-sm", "shapemapper2"),
]


def _datasets(threads: int) -> dict[str, tuple[str, external.CountParams, external.NormParams]]:
    """The four datasets these comparisons need: cmuts tuned to each tool (no
    spread), plus each tool's own reactivity."""
    count = dataclasses.replace(MAP, threads=threads)
    nospread = dataclasses.replace(count, cmuts_spread="nospread")
    return {
        "cmuts-match-rf": (
            "cmuts",
            nospread,
            external.NormParams(cmuts_norm="outlier", cmuts_per_reference=True),
        ),
        "cmuts-match-sm": ("cmuts", nospread, external.NormParams(cmuts_norm="sm")),
        "rnaframework": (
            "rnaframework",
            count,
            external.NormParams(rf_scoring=3, rf_norm_method=1),
        ),
        "shapemapper2": ("shapemapper2", count, external.NormParams()),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-f", "--fasta", required=True, type=Path)
    p.add_argument("--dms-mod", dest="mod_dms", required=True, type=Path)
    p.add_argument("--dms-nomod", dest="nomod_dms", required=True, type=Path)
    p.add_argument("--2a3-mod", dest="mod_2a3", required=True, type=Path)
    p.add_argument("--2a3-nomod", dest="nomod_2a3", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8, help="cmuts MPI ranks / rf-count workers")
    p.add_argument(
        "-o",
        "--out",
        type=Path,
        default=Path("outputs/tool-vs-tool.csv"),
        help="Output CSV (parent dir created); default is the gitignored outputs/.",
    )
    args = p.parse_args()

    refs = external.read_fasta(args.fasta)
    names = [n for n, _ in refs]
    length = max(len(s) for _, s in refs)
    inputs = {
        "DMS": external.Inputs(args.fasta, args.mod_dms, args.nomod_dms),
        "2A3": external.Inputs(args.fasta, args.mod_2a3, args.nomod_2a3),
    }
    built = external.build_profiles(
        _datasets(args.threads), inputs, sm_dir=os.environ.get("SM2_DIR")
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["comparison", "condition", "reference", "position", "cmuts", "tool", "abs_diff"]
        )
        for cmuts_ds, tool_ds in COMPARISONS:
            if cmuts_ds not in built or tool_ds not in built:
                print(f"  {tool_ds}: dataset unavailable, skipping")
                continue
            for cond in CONDITIONS:
                cmuts = external.to_array(built[cmuts_ds][cond], names, length)
                tool = external.to_array(built[tool_ds][cond], names, length)
                for i in range(len(names)):
                    cv, tv = cmuts[i], tool[i]
                    for pos in range(length):
                        c, t = cv[pos], tv[pos]
                        if np.isnan(c) or np.isnan(t):
                            continue
                        w.writerow(
                            [
                                tool_ds,
                                cond,
                                names[i],
                                pos,
                                f"{c:.6g}",
                                f"{t:.6g}",
                                f"{abs(c - t):.6g}",
                            ]
                        )
                        rows += 1
    print(f"wrote {rows} positions to {args.out}")


if __name__ == "__main__":
    main()
