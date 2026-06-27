"""How the benchmarks drive rf-count (RNAFramework): `rf-count` (count) then
`rf-norm` (normalize). rf-count reads only BAM, so CRAM/SAM are converted first.

`reactivity` at the bottom is the full per-tool pipeline. Several settings are
fixed for every benchmark and so are not carried in `CountParams` -- see
`count_command`.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import TYPE_CHECKING

from .common import CountParams, max_combined, run_checked, samtools_bin

if TYPE_CHECKING:
    import numpy as np


def available() -> bool:
    return (
        shutil.which(os.environ.get("RF_COUNT", "rf-count")) is not None
        and shutil.which(os.environ.get("RF_RCTOOLS", "rf-rctools")) is not None
    )


def count_command(bam, fasta, outdir, params: CountParams, *, overwrite: bool = True) -> list[str]:
    """Build the `rf-count` argv from `params`. The single home for the rf-count
    flag spelling the benchmarks rely on."""
    rf = os.environ.get("RF_COUNT", "rf-count")
    cmd = [rf, "-m", "-f", str(fasta), "-o", str(outdir), "-wt", str(params.threads)]
    if overwrite:
        cmd.append("--overwrite")
    # Settings fixed for every benchmark (not carried in CountParams): fast mode,
    # right-aligned deletions, duplicates kept, edit-distance 1.0, median-quality 0.
    cmd += ["--fast", "--right-deletion", "--no-discard-duplicates",
            "--max-edit-distance", "1.0", "--median-quality", "0"]  # fmt: skip
    if params.eval_surrounding:
        cmd.append("--eval-surrounding")
    if not params.insertions:
        cmd.append("--no-insertions")
    if not params.deletions:
        cmd.append("--no-deletions")
    if not params.cov_low_qual:
        cmd.append("--no-cov-low-qual")
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
    cmd.append(str(bam))
    return cmd


def run_count(bam, fasta, outdir, params: CountParams) -> Path:
    """Run rf-count and return the path to the resulting .rc file."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_checked(count_command(bam, fasta, outdir, params))
    base = Path(bam).name
    base = base[:-4] if base.endswith(".bam") else base
    return outdir / f"{base}.rc"


def run_norm(rc_path, outdir, *, scoring_method: int, norm_method: int, untreated=None) -> Path:
    """Normalize an rf-count .rc with rf-norm; returns the output directory.

    Coverage filters are off (`--mean-coverage 0 --median-coverage 0`). `untreated`
    is required for the Ding/Siegfried scoring methods.
    """
    rfn = os.environ.get("RF_NORM", "rf-norm")
    outdir = Path(outdir)
    cmd = [rfn, "--treated", str(rc_path), "-sm", str(scoring_method), "-nm", str(norm_method),
           "--mean-coverage", "0", "--median-coverage", "0", "--output-dir", str(outdir),
           "--overwrite"]  # fmt: skip
    if untreated is not None:
        cmd += ["--untreated", str(untreated)]
    run_checked(cmd)
    return outdir


def parse_norm(norm_dir) -> dict[str, np.ndarray]:
    """Parse an rf-norm output directory into {name: reactivity}. rf-norm writes
    one <id>.xml per transcript; the <reactivity> element is comma-separated."""
    import numpy as np

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


def coverage(rc_path) -> dict[str, np.ndarray]:
    """Per-position coverage per reference from a .rc (`rf-rctools view`); each
    transcript is four lines (id, sequence, counts, coverage)."""
    import numpy as np

    rc = os.environ.get("RF_RCTOOLS", "rf-rctools")
    proc = subprocess.run([rc, "view", str(rc_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"rf-rctools view failed: {proc.stderr[-2000:]}")
    out: dict[str, np.ndarray] = {}
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip() != ""]
    i = 0
    while i + 4 <= len(lines):
        out[lines[i].strip()] = np.array([float(x) for x in lines[i + 3].split(",") if x != ""])
        i += 4
    return out


def _to_bam(path, fasta, wd, samtools) -> Path:
    """rf-count needs a BAM; convert CRAM/SAM (passthrough if already BAM)."""
    path = Path(path)
    if path.suffix.lower() == ".bam":
        return path
    out = Path(wd) / f"{path.stem}.bam"
    run_checked([samtools, "view", "-b", "-T", str(fasta), "-o", str(out), str(path)])
    return out


def reactivity(count, norm, mod, nomod, fasta, wd):
    """rf-count full pipeline: convert to BAM -> rf-count -> rf-norm ->
    ({ref: reactivity}, {ref: depth}). Depth is the combined (mod + untreated) max
    per-position coverage."""
    wd = Path(wd)
    st = samtools_bin()
    mod = _to_bam(mod, fasta, wd, st)
    nomod = _to_bam(nomod, fasta, wd, st) if nomod is not None else None
    rc_mod = run_count(mod, fasta, wd / "rf-mod", count)
    rc_nomod = run_count(nomod, fasta, wd / "rf-nomod", count) if nomod is not None else None
    norm_dir = run_norm(
        rc_mod,
        wd / "rf-norm",
        scoring_method=norm.rf_scoring,
        norm_method=norm.rf_norm_method,
        untreated=rc_nomod,
    )
    cov_mod = coverage(rc_mod)
    cov_nomod = coverage(rc_nomod) if rc_nomod is not None else {}
    reads = {
        name: max_combined([c] + ([cov_nomod[name]] if name in cov_nomod else []))
        for name, c in cov_mod.items()
    }
    return parse_norm(norm_dir), reads
