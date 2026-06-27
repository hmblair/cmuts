"""Shared types and helpers for the benchmark tool layer.

`CountParams` / `NormParams` describe the count and normalize stages and `Inputs`
holds the alignments. The per-tool modules (`cmuts`, `rnaframework`,
`shapemapper`) build on these; the package `__init__` ties them together behind
`reactivity()` / `build_profiles()`.

Tools are resolved from the environment (CMUTS, RF_COUNT, RF_RCTOOLS, RF_NORM,
SAMTOOLS); shapemapper2 has no PATH convention, so its install dir is passed in.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

# numpy is imported lazily in the functions that build arrays, so paths that only
# build/run commands do not pay its import cost; annotations resolve here.
if TYPE_CHECKING:
    import numpy as np

TOOLS = ("cmuts", "rnaframework", "shapemapper2")

# The full pipeline imposes no coverage floor of its own -- it only blanks
# zero-coverage positions (the most permissive each tool allows) and exposes
# per-reference depth as /reads, so consumers filter as they choose.
FLOOR = 1


@dataclass
class CountParams:
    """Count-stage knobs, mapped to each tool's count flags by the per-tool
    command builders. Fields defaulting to None omit their flag (the tool's own
    default applies); a few apply to only one tool (noted inline)."""

    insertions: bool = True
    deletions: bool = True
    collapse: int = 2  # >1 -> collapse mods within (collapse - 1) of each other
    eval_surrounding: bool = True  # rf --eval-surrounding / cmuts --quality-window 1
    cov_low_qual: bool = False  # count low-quality positions toward coverage
    threads: int = 1
    min_mapq: int | None = None
    min_phred: int | None = None
    min_length: int | None = None
    max_length: int | None = None  # cmuts --max-length (max read length)
    max_indel: int | None = None
    cmuts_spread: str = "default"  # cmuts core: default | nospread | uniform
    max_internal_match: int | None = None  # shapemapper parser only


@dataclass
class NormParams:
    """Normalize-stage knobs. cmuts uses `cmuts_norm`/`cmuts_per_reference`;
    rf-count uses `rf_scoring`/`rf_norm_method`; shapemapper2 normalizes internally
    (DMS vs SHAPE from the condition). `cmuts_norm == "sm"` -> sm-dms/sm-shape."""

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


def strip_ext(path) -> str:
    """Drop a trailing .bam/.cram/.sam extension (cmuts sample identifiers)."""
    s = str(path)
    for ext in (".bam", ".cram", ".sam"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


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


def covered_references(bam, *, min_reads: int, samtools: str | None = None) -> set[str]:
    """References in `bam` with at least `min_reads` mapped reads (from idxstats).

    A position's depth cannot exceed a reference's total reads, so references
    below this floor can be skipped. Requires an indexed BAM/CRAM.
    """
    idx = subprocess.run(
        [samtools or samtools_bin(), "idxstats", str(bam)], capture_output=True, text=True
    )
    return {
        cols[0]
        for cols in (ln.split("\t") for ln in idx.stdout.splitlines())
        if len(cols) >= 3 and cols[2].isdigit() and int(cols[2]) >= min_reads
    }


def max_combined(arrays) -> float:
    """Max over positions of the element-wise sum of per-position depth arrays --
    the combined (mod + untreated) reference depth (cmuts' `reads` measure)."""
    import numpy as np

    arrays = [np.asarray(a, dtype=float) for a in arrays if len(a)]
    if not arrays:
        return 0.0
    total = np.zeros(max(len(a) for a in arrays))
    for a in arrays:
        total[: len(a)] += a
    return float(total.max())


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


def ensure_index(path, samtools: str) -> None:
    """Index a coordinate-sorted BAM/CRAM if its companion index is missing."""
    path = Path(path)
    suffix = ".crai" if path.suffix.lower() == ".cram" else ".bai"
    if not Path(str(path) + suffix).exists():
        run_checked([samtools, "index", str(path)])
