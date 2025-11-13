#!/usr/bin/env python3

from __future__ import annotations
import os
import h5py
import argparse
import cmuts


NAME = 'cmuts normalize'


parser = argparse.ArgumentParser()
parser.add_argument(
    'file',
    help='The HDF5 file containing the input data.',
)
parser.add_argument(
    '--mod',
    nargs="+",
    help='The names of the datasets containing the mutation counts of the modified sequences.',
    required=True,
)
parser.add_argument(
    '--nomod',
    nargs="+",
    help='The names of the datasets containing the mutation counts of the non-modified sequences.',
)
parser.add_argument(
    '--fasta',
    help='FASTA file of reference sequences.',
    required=True,
)
parser.add_argument(
    '--group',
    help='The group in the output file to place the reads and reactivities.',
    default="",
)
parser.add_argument(
    '-o', '--out',
    help='The output HDF5 file to write to.',
    default='reactivity.h5'
)
parser.add_argument(
    '--overwrite',
    help='Overwrite an existing HDF5 file (the whole file, not just the group).',
    action='store_true',
)
parser.add_argument(
    '--clip-low',
    help='Clip negative reactivity values.',
    action='store_true',
)
parser.add_argument(
    '--clip-high',
    help='Clip reactivity values above 1.',
    action='store_true',
)
parser.add_argument(
    '--blank-5p',
    help='NaN out this many bases on the 5\' end.',
    type=int,
    default=0,
)
parser.add_argument(
    '--blank-3p',
    help='NaN out this many bases on the 3\' end.',
    type=int,
    default=0,
)
parser.add_argument(
    '--blank-cutoff',
    help='NaN out any positions with less than this many reads.',
    type=int,
    default=10,
)
parser.add_argument(
    '--norm-independent',
    help='Normalize each profile separately, rather than using experiment-wide statistics.',
    action='store_true',
)
parser.add_argument(
    '--raw',
    help='Do not normalize the reactivity values.',
    action='store_true',
)
parser.add_argument(
    '--norm-outlier',
    help='Use the 2-8 outlier-based normalization method.',
    action='store_true',
)
parser.add_argument(
    '--norm-cutoff',
    help='References with at least this many reads are used for normalization.',
    type=int,
    default=500,
)
parser.add_argument(
    '--norm-percentile',
    help='The reactivity percentile to use for normalization.',
    type=int,
    default=90,
)
parser.add_argument(
    '--no-insertions',
    help='Do not use insertions to compute the reactivity profile.',
    action='store_true',
)
parser.add_argument(
    '--no-deletions',
    help='Do not use deletions to compute the reactivity profile.',
    action='store_true',
)
parser.add_argument(
    '--sig',
    help='The significance level for correlation plots.',
    type=float,
    default=0.05,
)


# ANSI escape codes (for printing)


BOLD = '\033[1m'
RESET = '\033[0m'


def _subtitle(name: str) -> str:
    if name:
        return " " * 8 + f"{BOLD}Statistics for {name}:{RESET}"
    else:
        return " " * 8 + f"{BOLD}Statistics:{RESET}"


def _remove_if_exists(
    path: str,
    overwrite: bool = False
) -> None:
    if os.path.exists(path):
        if overwrite:
            os.remove(path)


if __name__ == '__main__':

    args = parser.parse_args()

    print()
    print(cmuts.title(NAME, cmuts.version))
    print(_subtitle(args.group))

    _remove_if_exists(args.out, args.overwrite)

    clip = (args.clip_low, args.clip_high)
    blank = (args.blank_5p, args.blank_3p)

    opts = cmuts.Opts(
        cmuts.DataGroups(args.mod),
        cmuts.DataGroups(args.nomod),
        args.blank_cutoff,
        not args.no_insertions,
        not args.no_deletions,
        cmuts.NormScheme.UBR,
        blank,
        clip,
        args.sig,
    )

    with h5py.File(args.file, 'r') as f:

        mod, nomod, combined = cmuts.normalize(f, args.fasta, opts)

    # Save to file

    combined.save(args.group, args.out)

    # Print stats

    cmuts.stats(mod, nomod, combined)

    # Plot figures

    cmuts.plotting.generate(combined, args.group)
