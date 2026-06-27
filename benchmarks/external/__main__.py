"""CLI: `python -m external reactivity ...` runs one tool's full pipeline and
writes the HDF5, so the profile benchmark can time it in a subprocess."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import reactivity
from .common import TOOLS, CountParams, Inputs, NormParams


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="python -m external")
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
    raise SystemExit(main())
