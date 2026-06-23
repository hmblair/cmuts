#!/usr/bin/env python3
"""
Accuracy of reactivity profiles against 3D base-pairing, for cmuts vs rf-count
and shapemapper2 on the PDB130 dataset.

For each reference (named PDBID_chain_wt in the cmuts reference FASTA), match it
by sequence to an RNA chain in the deposited PDB structure, detect base-base
hydrogen bonds geometrically from the 3D coordinates, and derive a probe-matched
"unpaired" ground truth:

  - 2A3 / SHAPE : a nucleotide is "paired" if its base participates in ANY
                  base-base hydrogen bond (canonical or not). SHAPE reports
                  backbone flexibility, which any base contact constrains.
  - DMS         : scored only at A and C. A nucleotide is "protected" if its
                  Watson-Crick-face atom (A-N1 / C-N3, the atoms DMS methylates)
                  is specifically hydrogen bonded. G/U are not DMS-probed (NaN).

Base pairing is computed from geometry (heavy-atom N/O distances), not from
deposited struct_conn records, because many entries (especially cryo-EM
ribosomes) do not annotate base pairs.

All three tools' profiles are computed from BAMs and scored together. Each
condition (DMS, 2A3) takes a modified and an unmodified (nomod) sample. cmuts is
run for every spread mode (no-spread / mutation-informed default / uniform) via
`cmuts core` + `cmuts normalize --mod --nomod`; rf-count via rf-norm and
shapemapper2 via its make_reactivity_profiles/normalize_profiles step (each its
own native normalization, see `RF_PARAMS` / `SM_PARAMS`).
Alternatively, precomputed cmuts profiles can be read from an HDF5 (--profiles)
instead of recomputing them. Each (mode/tool x condition) profile is scored by
AUC / Pearson / Spearman of reactivity against the unpaired ground truth over
the core region.

Requires ciffy, biopython, and scipy in addition to the cmuts deps. BAMs must be
sorted and indexed (the per-reference shapemapper2 loop needs the index).

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import os
import tempfile
import warnings
from functools import lru_cache
from pathlib import Path

import ciffy
import external
import h5py
import numpy as np
from Bio import Align
from ciffy.biochemistry import Residue
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

MODES = ["nospread", "default", "uniform"]  # cmuts spread modes
TOOLS = ["rf-count", "shapemapper2"]  # external tools
CONDITIONS = ["DMS", "2A3"]
MIN_VALID = 10  # minimum scored positions per reference

# cmuts core flag per spread mode ("" = the mutation-informed default).
SPREAD_FLAG = {"nospread": "--no-spread", "default": "", "uniform": "--uniform-spread"}

# Library adapter blanking and per-reference read cap for cmuts (PDB130 defaults).
DEFAULT_BLANK_5P = 26
DEFAULT_BLANK_3P = 20
DEFAULT_DOWNSAMPLE = 50000

# rf-count counting tuned for accuracy, matching the reference comparison
# pipeline (papers/cmuts scratch/comparison): mutational profiling WITHOUT
# deletions (rf-count cannot spread them), right-aligned, duplicates kept,
# quality floor 10. Normalization is rf-norm (treated-only Zubradt, 2-8%).
RF_PARAMS = external.Params(
    insertions=True,
    deletions=False,
    right_align_deletions=True,
    discard_duplicates=False,
    collapse=2,
    eval_surrounding=True,
    cov_low_qual=True,
    min_mapq=10,
    min_phred=10,
    min_length=2,
    median_quality=0,
    max_edit_distance=1.0,
)

# shapemapper2 parser settings; its make_reactivity_profiles + normalize_profiles
# turn the counts into a (background-subtracted, normalized) reactivity profile.
SM_PARAMS = external.Params(
    right_align_deletions=True,
    min_mapq=10,
    min_phred=10,
    max_internal_match=2,
)

# Base N/O atoms that participate in base-base hydrogen bonds (ciffy atom codes).
BASE_ATOMS = {
    int(getattr(getattr(Residue, b), a))
    for b, atoms in {
        "A": ("N1", "N3", "N6", "N7"),
        "G": ("N1", "N2", "N3", "N7", "O6"),
        "C": ("N3", "N4", "O2"),
        "U": ("N3", "O2", "O4"),
    }.items()
    for a in atoms
}
# Watson-Crick-face atoms DMS methylates (protected only when H-bonded).
WC_FACE = {int(Residue.A.N1), int(Residue.C.N3)}
A_VAL, C_VAL = int(Residue.A.value), int(Residue.C.value)
HBOND_CUTOFF = 3.4  # heavy-atom N/O...N/O distance (Angstrom)


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(name, sequence)] in file order."""
    entries: list[tuple[str, str]] = []
    name, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                entries.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        entries.append((name, "".join(seq)))
    return entries


def parse_header(name: str) -> tuple[str, str]:
    """'5TZS_z_wt' -> ('5TZS', 'z')."""
    parts = name.split("_")
    return parts[0], parts[1]


# ---------------------------------------------------------------------------
# Ground-truth pairing labels from 3D structure
# ---------------------------------------------------------------------------


@lru_cache(maxsize=8)
def _load_structure(pdbid: str, struct_dir: str):
    cif = Path(struct_dir) / f"{pdbid}.cif"
    if not cif.exists():
        return None
    try:
        return ciffy.load(str(cif), skip=["descriptions"])
    except Exception:
        return None


def _rna_chains(s):
    """Yield (sequence, chain_polymer) for each RNA chain."""
    for c in s.chains():
        mt = np.atleast_1d(np.asarray(c.molecule_types)).ravel()
        if mt.size and int(mt[0]) == int(ciffy.RNA):
            seq = c.sequence_str()
            if seq:
                yield seq, c


def _geometry_labels(chain) -> tuple[np.ndarray, np.ndarray]:
    """Probe-matched unpaired labels over a chain's residues, from geometry.

    Returns (shape_unpaired, dms_unpaired):
      shape_unpaired[r] = 0 if residue r has any base-base H-bond else 1.
      dms_unpaired[r]   = 0/1 for A,C by whether the WC-face atom is H-bonded;
                          NaN for G,U (not DMS-probed).
    """
    coords = np.asarray(chain.coordinates)
    atype = np.asarray(chain.atoms)
    resid = np.asarray(chain.membership(ciffy.RESIDUE))
    seq = np.asarray(chain.sequence)
    n = len(seq)

    sel = np.isin(atype, list(BASE_ATOMS))
    P, A, R = coords[sel], atype[sel], resid[sel]

    base_paired = np.zeros(n, dtype=bool)
    wc_hbonded = np.zeros(n, dtype=bool)
    if len(P):
        for i, j in cKDTree(P).query_pairs(HBOND_CUTOFF):
            ri, rj = int(R[i]), int(R[j])
            if abs(ri - rj) < 2:  # skip self/adjacent (stacking, backbone)
                continue
            base_paired[ri] = base_paired[rj] = True
            if A[i] in WC_FACE:
                wc_hbonded[ri] = True
            if A[j] in WC_FACE:
                wc_hbonded[rj] = True

    shape_unpaired = (~base_paired).astype(float)
    dms_unpaired = np.full(n, np.nan)
    is_ac = np.isin(seq, [A_VAL, C_VAL])
    dms_unpaired[is_ac] = (~wc_hbonded[is_ac]).astype(float)
    return shape_unpaired, dms_unpaired


def structure_labels(pdbid: str, ref_seq: str, struct_dir: str):
    """(struct_seq, shape_unpaired, dms_unpaired) for the structure's RNA chain
    that best matches the reference by sequence, or None."""
    s = _load_structure(pdbid, struct_dir)
    if s is None:
        return None
    target = ref_seq.upper().replace("U", "T")
    best = None  # (score, seq, chain)
    for seq, c in _rna_chains(s):
        cand = seq.upper().replace("U", "T")
        if set(cand) <= {"N"}:
            continue
        try:
            sc = _ALIGNER.score(target, cand)
        except Exception:
            continue
        if best is None or sc > best[0]:
            best = (sc, seq, c)
    if best is None:
        return None
    _, seq, c = best
    try:
        shape_unp, dms_unp = _geometry_labels(c)
    except Exception:
        return None
    return seq, shape_unp, dms_unp


_ALIGNER = Align.PairwiseAligner()
_ALIGNER.mode = "global"
_ALIGNER.match_score = 2
_ALIGNER.mismatch_score = -1
_ALIGNER.open_gap_score = -5
_ALIGNER.extend_gap_score = -0.5


def ref_to_struct_map(ref_seq: str, struct_seq: str) -> dict[int, int]:
    """Map reference position -> structure residue index via global alignment."""
    a, b = ref_seq.upper().replace("U", "T"), struct_seq.upper().replace("U", "T")
    aln = _ALIGNER.align(a, b)[0]
    mapping: dict[int, int] = {}
    ri = si = 0
    for r, s in zip(aln[0], aln[1]):
        if r != "-" and s != "-":
            mapping[ri] = si
            ri += 1
            si += 1
        elif r == "-":
            si += 1
        else:
            ri += 1
    return mapping


def project(n_ref: int, mapping: dict[int, int], labels: np.ndarray) -> np.ndarray:
    out = np.full(n_ref, np.nan)
    for ri, si in mapping.items():
        if 0 <= si < len(labels):
            out[ri] = labels[si]
    return out


# ---------------------------------------------------------------------------
# Scoring: reactivity vs unpaired ground truth
# ---------------------------------------------------------------------------


def _auc(y: np.ndarray, x: np.ndarray) -> float:
    """ROC AUC = P(score of an unpaired position > score of a paired one), via
    the rank-based (Mann-Whitney) formula, tie-safe through rankdata."""
    r = stats.rankdata(x)
    n_pos = float(y.sum())
    n_neg = float(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.nan
    return (r[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def score_profile(react: np.ndarray, unpaired: np.ndarray) -> dict | None:
    mask = ~np.isnan(react) & ~np.isnan(unpaired)
    if mask.sum() < MIN_VALID:
        return None
    y, x = unpaired[mask], react[mask]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pear = stats.pearsonr(x, y)[0]
        spear = stats.spearmanr(x, y)[0]
    return {"n": int(mask.sum()), "auc": _auc(y, x), "pearson": pear, "spearman": spear}


# ---------------------------------------------------------------------------
# Reactivity profiles from the three tools
# ---------------------------------------------------------------------------


def _strip_ext(bam: Path) -> str:
    """cmuts identifies a sample by its path without the alignment extension."""
    s = str(bam)
    for ext in (".bam", ".cram", ".sam"):
        if s.endswith(ext):
            return s[: -len(ext)]
    return s


def cmuts_profiles(
    names: list[str],
    fasta: Path,
    cond_samples: dict[str, tuple[Path | None, Path | None]],
    *,
    blank_5p: int,
    blank_3p: int,
    downsample: int,
    workdir: Path,
) -> dict[tuple[str, str], dict[str, np.ndarray]]:
    """Compute cmuts reactivity for every spread mode x condition from BAMs.

    Mirrors the established pipeline: one `cmuts core` per spread mode over all
    samples (--no-insertions plus the mode's spread flag), then `cmuts normalize
    --mod --nomod` per condition. Returns {(mode, condition): {ref_name: react}}.
    """
    cmuts = os.environ.get("CMUTS", "cmuts")
    all_bams = [b for pair in cond_samples.values() for b in pair if b is not None]
    out: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for mode in MODES:
        counts = workdir / f"counts-{mode}.h5"
        prof = workdir / f"profiles-{mode}.h5"
        core = [
            cmuts,
            "core",
            "-f",
            str(fasta),
            "-o",
            str(counts),
            "--overwrite",
            "--no-insertions",
        ]
        if SPREAD_FLAG[mode]:
            core.append(SPREAD_FLAG[mode])
        if downsample:
            core += ["--downsample", str(downsample)]
        core += [str(b) for b in all_bams]
        print(f"  cmuts core [{mode}] ...")
        external.run_checked(core)

        for cond, (mod, nomod) in cond_samples.items():
            if mod is None:
                continue
            ncmd = [
                cmuts,
                "normalize",
                "-o",
                str(prof),
                "--fasta",
                str(fasta),
                "--mod",
                _strip_ext(mod),
            ]
            if nomod is not None:
                ncmd += ["--nomod", _strip_ext(nomod)]
            ncmd += [
                "--group",
                f"{mode}-{cond}",
                "--blank-5p",
                str(blank_5p),
                "--blank-3p",
                str(blank_3p),
                "--clip-below",
                "0",
                str(counts),
            ]
            print(f"  cmuts normalize [{mode}/{cond}] ...")
            external.run_checked(ncmd)
            with h5py.File(prof, "r") as f:
                react = f[f"{mode}-{cond}/reactivity"][:]
            out[(mode, cond)] = {name: react[i] for i, name in enumerate(names)}
    return out


def tool_profiles(
    refs: list[tuple[str, str]],
    fasta: Path,
    cond_samples: dict[str, tuple[Path | None, Path | None]],
) -> dict[tuple[str, str], dict[str, np.ndarray]]:
    """Native reactivity profiles from rf-count and shapemapper2 per condition,
    as {(tool, condition): {ref_name: reactivity}}.

    rf-count is normalized with rf-norm (modified sample only); shapemapper2 uses
    its own make_reactivity_profiles / normalize_profiles step (modified +
    untreated). Tools that are unavailable (rf-count/rf-norm not on PATH, SM2_DIR
    unset) are skipped, so this returns only what could be computed.
    """
    sm_dir = os.environ.get("SM2_DIR")
    out: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for cond, (mod, nomod) in cond_samples.items():
        if mod is None:
            continue
        workdir = Path(tempfile.mkdtemp(prefix=f"accuracy-{cond.lower()}-"))
        # rf-count -> rf-norm (native normalization; treated-only Zubradt, 2-8%).
        if external.rfcount_available() and external.rfnorm_available():
            print(f"  rf-count + rf-norm [{cond}] ...")
            rc = external.run_rfcount(mod, fasta, workdir / "rf", RF_PARAMS)
            norm_dir = external.run_rfnorm(rc, workdir / "rf_norm")
            out[("rf-count", cond)] = external.parse_rfnorm(norm_dir)
        else:
            print(f"  rf-count / rf-norm not found; skipping [{cond}]")
        # shapemapper2 native reactivity (modified + untreated control).
        if sm_dir:
            print(f"  shapemapper2 reactivity [{cond}] ...")
            out[("shapemapper2", cond)] = external.run_shapemapper_reactivity(
                mod,
                nomod,
                refs,
                sm_dir,
                SM_PARAMS,
                dms=(cond == "DMS"),
                workdir=workdir / "sm",
            )
        else:
            print(f"  SM2_DIR unset; skipping shapemapper2 [{cond}]")
    return out


# ---------------------------------------------------------------------------
# Entry point: build profiles, score against structure, write tables
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--structures", required=True, type=Path)
    ap.add_argument(
        "--dms-mod",
        dest="mod_dms",
        type=Path,
        default=None,
        help="Modified BAM for the DMS condition",
    )
    ap.add_argument(
        "--dms-nomod",
        dest="nomod_dms",
        type=Path,
        default=None,
        help="Unmodified (nomod) BAM for the DMS condition",
    )
    ap.add_argument(
        "--2a3-mod",
        dest="mod_2a3",
        type=Path,
        default=None,
        help="Modified BAM for the 2A3 condition",
    )
    ap.add_argument(
        "--2a3-nomod",
        dest="nomod_2a3",
        type=Path,
        default=None,
        help="Unmodified (nomod) BAM for the 2A3 condition",
    )
    ap.add_argument(
        "--profiles",
        type=Path,
        default=None,
        help="Precomputed cmuts profiles HDF5 (skips recomputing cmuts from BAMs)",
    )
    ap.add_argument("--blank-5p", type=int, default=DEFAULT_BLANK_5P)
    ap.add_argument("--blank-3p", type=int, default=DEFAULT_BLANK_3P)
    ap.add_argument(
        "--downsample",
        type=int,
        default=DEFAULT_DOWNSAMPLE,
        help="Per-reference read cap for cmuts core (0 to disable)",
    )
    ap.add_argument("--per-ref", type=Path, default=Path("per_reference.tsv"))
    ap.add_argument("--summary", type=Path, default=Path("summary.tsv"))
    args = ap.parse_args()

    refs = load_fasta(args.fasta)
    names = [n for n, _ in refs]
    cond_samples: dict[str, tuple[Path | None, Path | None]] = {
        "DMS": (args.mod_dms, args.nomod_dms),
        "2A3": (args.mod_2a3, args.nomod_2a3),
    }
    have_bams = any(mod is not None for mod, _ in cond_samples.values())
    if not have_bams and args.profiles is None:
        ap.error("provide per-condition BAMs (--dms-mod/--2a3-mod ...) or --profiles")

    # Every profile, keyed by (label, condition) -> {ref_name: reactivity}.
    profiles: dict[tuple[str, str], dict[str, np.ndarray]] = {}

    # cmuts: compute every spread mode from the BAMs, or load a precomputed h5.
    if have_bams:
        print("computing cmuts profiles ...")
        cmuts_workdir = Path(tempfile.mkdtemp(prefix="accuracy-cmuts-"))
        profiles.update(
            cmuts_profiles(
                names,
                args.fasta,
                cond_samples,
                blank_5p=args.blank_5p,
                blank_3p=args.blank_3p,
                downsample=args.downsample,
                workdir=cmuts_workdir,
            )
        )
    elif args.profiles is not None:
        with h5py.File(args.profiles, "r") as f:
            for mode in MODES:
                for cond in CONDITIONS:
                    key = f"{mode}-{cond}/reactivity"
                    if key in f:
                        arr = f[key][:]
                        profiles[(mode, cond)] = {
                            names[i]: arr[i] for i in range(min(len(names), arr.shape[0]))
                        }

    # rf-count + shapemapper2 from the modified sample of each condition.
    if have_bams:
        print("running external tools ...")
        profiles.update(tool_profiles(refs, args.fasta, cond_samples))

    rows: list[dict] = []
    skipped = {"no_structure": 0, "low_coverage": 0}
    for name, ref_seq in refs:
        pdbid, chain = parse_header(name)
        labels = structure_labels(pdbid, ref_seq, str(args.structures))
        if labels is None:
            skipped["no_structure"] += 1
            continue
        struct_seq, shape_unp, dms_unp = labels
        mapping = ref_to_struct_map(ref_seq, struct_seq)
        gt = {
            "2A3": project(len(ref_seq), mapping, shape_unp),
            "DMS": project(len(ref_seq), mapping, dms_unp),
        }

        any_scored = False
        for (label, cond), per_ref in profiles.items():
            arr = per_ref.get(name)
            if arr is None:
                continue
            res = score_profile(np.asarray(arr, dtype=float), gt[cond])
            if res is None:
                continue
            any_scored = True
            rows.append(
                {
                    "reference": name,
                    "pdb": pdbid,
                    "chain": chain,
                    "mode": label,
                    "condition": cond,
                    **res,
                }
            )
        if not any_scored:
            skipped["low_coverage"] += 1

    cols = ["reference", "pdb", "chain", "mode", "condition", "n", "auc", "pearson", "spearman"]
    with open(args.per_ref, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    def agg(metric: str, mode: str, cond: str) -> tuple[int, float, float]:
        vals = np.array(
            [
                r[metric]
                for r in rows
                if r["mode"] == mode and r["condition"] == cond and not np.isnan(r[metric])
            ]
        )
        if len(vals) == 0:
            return 0, np.nan, np.nan
        return len(vals), float(np.mean(vals)), float(np.median(vals))

    order = MODES + TOOLS
    present = [m for m in order if any(r["mode"] == m for r in rows)]
    with open(args.summary, "w") as f:
        f.write("condition\tmetric\tmode\tn\tmean\tmedian\tmean_minus_default\n")
        for cond in CONDITIONS:
            for metric in ("auc", "pearson", "spearman"):
                _, def_mean, _ = agg(metric, "default", cond)
                for mode in present:
                    n, mean, median = agg(metric, mode, cond)
                    delta = (
                        mean - def_mean if not (np.isnan(mean) or np.isnan(def_mean)) else np.nan
                    )
                    f.write(
                        f"{cond}\t{metric}\t{mode}\t{n}\t{mean:.4f}\t{median:.4f}\t{delta:+.4f}\n"
                    )

    n_scored = len({r["reference"] for r in rows})
    print(
        f"references: {len(refs)}  scored: {n_scored}  "
        f"no_structure: {skipped['no_structure']}  low_coverage: {skipped['low_coverage']}"
    )
    print(f"per-reference: {args.per_ref}")
    print(f"summary:       {args.summary}\n")
    print(Path(args.summary).read_text())


if __name__ == "__main__":
    main()
