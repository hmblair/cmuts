#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys
from typing import Any

import h5py

import cmuts

NAME = "cmuts normalize"


parser = argparse.ArgumentParser()
parser.add_argument(
    "file",
    help="The HDF5 file containing the input data.",
)
parser.add_argument(
    "--mod",
    nargs="+",
    help=(
        "Single-group mode: HDF5 paths for the modified condition "
        "(replicates are summed). Mutually exclusive with --groups."
    ),
)
parser.add_argument(
    "--nomod",
    nargs="+",
    help=(
        "Single-group mode: HDF5 paths for the unmodified control. "
        "Mutually exclusive with --groups."
    ),
)
parser.add_argument(
    "--group",
    help=(
        "Single-group mode: name for the output HDF5 group. " "Mutually exclusive with --groups."
    ),
    default="",
)
parser.add_argument(
    "--groups",
    help=(
        "Multi-group mode: path to a TOML file describing the groups. "
        "See the documentation for the file format."
    ),
)
parser.add_argument(
    "--independent-norm",
    help=(
        "Normalize each group independently rather than pooling across "
        "groups. Has no effect in single-group mode."
    ),
    action="store_true",
)
parser.add_argument(
    "--fasta",
    help="FASTA file of reference sequences.",
    required=True,
)
parser.add_argument(
    "-o", "--out", help="The output HDF5 file to write to.", default="reactivity.h5"
)
parser.add_argument(
    "--overwrite",
    help="Overwrite an existing HDF5 file (the whole file, not just the group).",
    action="store_true",
)
parser.add_argument(
    "--clip-below",
    help="Clip reactivity values below this threshold (e.g. 0).",
    type=float,
    default=None,
    metavar="N",
)
parser.add_argument(
    "--clip-above",
    help="Clip reactivity values above this threshold (e.g. 1).",
    type=float,
    default=None,
    metavar="N",
)
parser.add_argument(
    "--blank-5p",
    help="NaN out this many bases on the 5' end.",
    type=int,
    default=0,
)
parser.add_argument(
    "--blank-3p",
    help="NaN out this many bases on the 3' end.",
    type=int,
    default=0,
)
parser.add_argument(
    "--blank-cutoff",
    help="NaN out any positions with less than this many reads.",
    type=int,
    default=10,
)
parser.add_argument(
    "--norm",
    help="Normalization method: ubr (90th percentile), raw (none), outlier (2-8%%).",
    choices=["ubr", "raw", "outlier"],
    default="ubr",
)
parser.add_argument(
    "--norm-cutoff",
    help="References with at least this many reads are used for normalization.",
    type=int,
    default=500,
)
parser.add_argument(
    "--norm-percentile",
    help="The reactivity percentile to use for normalization.",
    type=int,
    default=90,
)
parser.add_argument(
    "--no-insertions",
    help="Do not use insertions to compute the reactivity profile.",
    action="store_true",
)
parser.add_argument(
    "--no-deletions",
    help="Do not use deletions to compute the reactivity profile.",
    action="store_true",
)
parser.add_argument(
    "--sig",
    help="The significance level for correlation plots.",
    type=float,
    default=0.05,
)


def _load_toml(path: str) -> dict[str, Any]:
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ImportError:
            sys.exit(
                "Reading TOML group files requires Python 3.11+ or the 'tomli' "
                "package. Install with: pip install tomli"
            )
    with open(path, "rb") as f:
        return tomllib.load(f)


def _groups_from_toml(path: str) -> list[cmuts.Group]:
    """Parse a multi-group TOML config into a list of Group objects.

    Expected format::

        [[group]]
        name = "with_ligand"
        mod = ["mod1", "mod2"]
        nomod = ["nomod1"]  # optional
    """
    data = _load_toml(path)
    raw = data.get("group")
    if raw is None:
        sys.exit(f"{path}: no [[group]] entries found.")
    if not isinstance(raw, list):
        sys.exit(f"{path}: 'group' must be an array of tables ([[group]]).")

    groups: list[cmuts.Group] = []
    for i, entry in enumerate(raw):
        if not isinstance(entry, dict):
            sys.exit(f"{path}: [[group]] entry {i} is not a table.")
        name = entry.get("name")
        mod = entry.get("mod")
        nomod = entry.get("nomod")
        if not isinstance(name, str) or not name:
            sys.exit(f"{path}: [[group]] entry {i} is missing 'name'.")
        if not isinstance(mod, list) or not mod or not all(isinstance(m, str) for m in mod):
            sys.exit(f"{path}: group '{name}' needs a non-empty 'mod' list of strings.")
        if nomod is not None and (
            not isinstance(nomod, list) or not all(isinstance(n, str) for n in nomod)
        ):
            sys.exit(f"{path}: group '{name}' has invalid 'nomod' (must be list of strings).")
        groups.append(cmuts.Group(name=name, mod=list(mod), nomod=list(nomod) if nomod else None))

    return groups


def _remove_if_exists(path: str, overwrite: bool = False) -> None:
    if os.path.exists(path):
        if overwrite:
            os.remove(path)


def main():
    args = parser.parse_args()

    using_groups_file = args.groups is not None
    using_single = args.mod is not None
    if using_groups_file and using_single:
        parser.error("--groups is mutually exclusive with --mod/--nomod/--group.")
    if not using_groups_file and not using_single:
        parser.error("Provide either --mod (single group) or --groups (multi-group).")

    if using_groups_file:
        groups = _groups_from_toml(args.groups)
    else:
        groups = [
            cmuts.Group(
                name=args.group, mod=list(args.mod), nomod=list(args.nomod) if args.nomod else None
            )
        ]

    print()
    print(cmuts.title(NAME, cmuts.__version__))
    if len(groups) == 1:
        print(cmuts.subtitle(groups[0].name))
    else:
        print(cmuts.subtitle(", ".join(g.name for g in groups)))

    _remove_if_exists(args.out, args.overwrite)

    clip = (args.clip_below, args.clip_above)
    blank = (args.blank_5p, args.blank_3p)

    # opts.mod and opts.nomod are unused by compute_reactivities; pass empties.
    opts = cmuts.Opts(
        cmuts.DataGroups([]),
        cmuts.DataGroups(None),
        args.blank_cutoff,
        not args.no_insertions,
        not args.no_deletions,
        args.norm,
        blank,
        clip,
        args.sig,
    )

    with h5py.File(args.file, "r") as f:
        results = cmuts.compute_reactivities(
            f,
            args.fasta,
            groups,
            opts,
            shared_norm=not args.independent_norm,
        )

    # Save all groups to a single output file.

    cmuts.save_groups(args.out, [(r.group.name, r.combined) for r in results])

    # Print stats and plot per group.

    figdir = os.path.join(os.path.dirname(os.path.abspath(args.out)), "figures")
    for r in results:
        if len(results) > 1:
            print()
            print(cmuts.subtitle(r.group.name))
        cmuts.stats(r.mod, r.nomod, r.combined)
        cmuts.visualize.plot_all(r.combined, r.group.name, figdir)
        cmuts.visualize.plot_snr_scaling(r.mod, r.nomod, r.combined, r.group.name, figdir)


if __name__ == "__main__":
    main()
