#!/usr/bin/env python3
"""Per-nucleotide correctness comparison: cmuts (tuned to match each existing
tool) vs that tool, from the precomputed benchmark profiles.

With deletion spreading off and its settings tuned to an existing pipeline, cmuts
should reproduce that pipeline's reactivity profile nucleotide-by-nucleotide --
the Fig 2a control that isolates deletion spreading (cmuts reproduces prior
tools), not a standalone accuracy result. This reads the profiles HDF5 written by
`profiles.py` and writes, per (comparison, condition, reference, position), the
two reactivities and their absolute difference.

Comparisons (cmuts-tuned-to-X vs X):
  rnaframework   cmuts-match-rf  vs rnaframework
  shapemapper2   cmuts-match-sm  vs shapemapper2   (skipped until that dataset exists)

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np

CONDITIONS = ["DMS", "2A3"]

# (cmuts dataset tuned to the tool, the tool dataset); the comparison is labelled
# by the tool. cmuts-match-sm is a TODO in profiles.py and is skipped if absent.
COMPARISONS = [
    ("cmuts-match-rf", "rnaframework"),
    ("cmuts-match-sm", "shapemapper2"),
]


def read_fasta_names(path: str) -> list[str]:
    """Reference names (header up to the first space), in file order."""
    names: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    return names


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profiles", required=True, type=Path, help="Profiles HDF5 from profiles.py")
    p.add_argument("-f", "--fasta", required=True, help="Reference FASTA (for reference names)")
    p.add_argument(
        "-o",
        "--out",
        type=Path,
        default=Path("outputs/correctness.csv"),
        help="Output CSV (parent dir created); default is the gitignored outputs/.",
    )
    args = p.parse_args()

    names = read_fasta_names(args.fasta)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with h5py.File(args.profiles, "r") as f, open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["comparison", "condition", "reference", "position", "cmuts", "tool", "abs_diff"]
        )
        for cmuts_ds, tool_ds in COMPARISONS:
            for cond in CONDITIONS:
                ck, tk = f"{cmuts_ds}/{cond}", f"{tool_ds}/{cond}"
                if ck not in f or tk not in f:
                    print(f"  {tool_ds}/{cond}: dataset missing, skipping")
                    continue
                cmuts, tool = f[ck][:], f[tk][:]
                n = min(cmuts.shape[0], tool.shape[0], len(names))
                for i in range(n):
                    cv, tv = cmuts[i], tool[i]
                    for pos in range(min(len(cv), len(tv))):
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
