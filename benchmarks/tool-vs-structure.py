#!/usr/bin/env python3
"""
Accuracy of reactivity profiles against 3D base-pairing, for cmuts vs rf-count
and shapemapper2 on the PDB130 dataset.

For each reference (named PDBID_chain_wt in the cmuts reference FASTA), match it
by sequence to an RNA chain in the deposited PDB structure, read its base-base
hydrogen bonds from the deposited annotations, and derive a probe-matched
"unpaired" ground truth:

  - 2A3 / SHAPE : a nucleotide is "paired" if its base participates in ANY
                  base-base hydrogen bond (canonical or not). SHAPE reports
                  backbone flexibility, which any base contact constrains.
  - DMS         : scored only at A and C. A nucleotide is "protected" if its
                  Watson-Crick-face atom (A-N1 / C-N3, the atoms DMS methylates)
                  is specifically hydrogen bonded. G/U are not DMS-probed (NaN).

Base pairing is read from the deposited hydrogen-bond annotations (mmCIF
struct_conn 'hydrog' records, via ciffy `connections`), which are present for
~all of PDB130 (incl. cryo-EM ribosomes). References whose structure carries no
such annotations are skipped.

The reactivity profiles are built here from the modified/untreated alignments via
`external.reactivity` (one dataset per cmuts spread mode and per existing tool, x
condition). Each (dataset x condition) profile is scored by AUC / Pearson /
Spearman of reactivity against the unpaired ground truth over the core region,
and the summary reports each dataset's delta against cmuts' default spread mode,
averaged over the references that every dataset scored (so coverage filters that
drop different references per tool do not bias the comparison).

Builds each tool's profile from the alignments, so it needs samtools, rf-count,
and (for shapemapper2) SM2_DIR. Requires ciffy, biopython, and scipy too.

Hamish M. Blair, 2026
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import tempfile
import warnings
from functools import lru_cache
from pathlib import Path

import ciffy
import external
import numpy as np
from Bio import Align
from ciffy.biochemistry import Residue
from scipy import stats

warnings.filterwarnings("ignore")

CONDITIONS = ["DMS", "2A3"]
MIN_VALID = 10  # minimum scored positions per reference

# Datasets this benchmark scores: cmuts at its three spread modes plus the two
# existing tools. Built from the alignments by `_datasets` below.
DATASETS = ["cmuts-nospread", "cmuts-uniform", "cmuts-default", "rnaframework", "shapemapper2"]
BASELINE = "cmuts-default"  # the summary reports each dataset's delta against this

# One map-seq counting config shared by every tool (cmuts spread set per dataset).
MAP = external.CountParams(
    insertions=False,
    deletions=True,
    collapse=2,
    eval_surrounding=True,
    quality_window=1,
    cov_low_qual=True,
    min_mapq=10,
    min_phred=10,
    min_length=2,
    max_internal_match=2,
)


def _datasets(threads: int) -> dict[str, tuple[str, external.CountParams, external.NormParams]]:
    """The five datasets scored here: cmuts at each spread mode (ubr norm), plus
    rf-count and shapemapper2."""
    count = dataclasses.replace(MAP, threads=threads)

    def cmuts(spread: str) -> tuple[str, external.CountParams, external.NormParams]:
        return (
            "cmuts",
            dataclasses.replace(count, cmuts_spread=spread),
            external.NormParams(cmuts_norm="ubr"),
        )

    return {
        "cmuts-nospread": cmuts("nospread"),
        "cmuts-uniform": cmuts("uniform"),
        "cmuts-default": cmuts("default"),
        "rnaframework": (
            "rnaframework",
            count,
            external.NormParams(rf_scoring=3, rf_norm_method=1),
        ),
        "shapemapper2": ("shapemapper2", count, external.NormParams()),
    }


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


def load_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(name, sequence)] in file order."""
    entries: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
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


def _pairing_labels(chain) -> tuple[np.ndarray, np.ndarray] | None:
    """Probe-matched unpaired labels from the deposited hydrogen-bond annotations.

    Reads ciffy `connections` (the mmCIF struct_conn 'hydrog' records, atom-level)
    rather than re-deriving contacts geometrically: the deposited annotations are
    present for ~all of PDB130 (incl. cryo-EM ribosomes) and avoid the
    over-calling of a bare distance cutoff. Returns None for the rare structure
    with no connections, so the caller skips that reference.

    Returns (shape_unpaired, dms_unpaired):
      shape_unpaired[r] = 0 if residue r's base participates in any base-base
                          hydrogen bond else 1.
      dms_unpaired[r]   = 0/1 for A,C by whether the WC-face atom (A-N1/C-N3) is
                          hydrogen bonded; NaN for G,U (not DMS-probed).
    """
    conns = np.asarray(chain.connections)
    if conns.ndim != 2 or conns.shape[0] == 0:
        return None  # no deposited annotations -> skip this reference

    atype = np.asarray(chain.atoms)
    resid = np.asarray(chain.membership(ciffy.RESIDUE))
    seq = np.asarray(chain.sequence)
    n = len(seq)

    base_paired = np.zeros(n, dtype=bool)
    wc_hbonded = np.zeros(n, dtype=bool)
    for a, b in conns:
        # base-base hydrogen bonds only (both partners are base N/O atoms)
        if atype[a] in BASE_ATOMS and atype[b] in BASE_ATOMS:
            ra, rb = int(resid[a]), int(resid[b])
            base_paired[ra] = base_paired[rb] = True
            if atype[a] in WC_FACE:
                wc_hbonded[ra] = True
            if atype[b] in WC_FACE:
                wc_hbonded[rb] = True

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
        labels = _pairing_labels(c)
    except Exception:
        return None
    if labels is None:  # structure carries no hydrogen-bond annotations
        return None
    shape_unp, dms_unp = labels
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
# Entry point: load profiles, score against structure, write tables
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True, type=Path)
    ap.add_argument("--structures", required=True, type=Path)
    ap.add_argument("--dms-mod", dest="mod_dms", required=True, type=Path)
    ap.add_argument("--dms-nomod", dest="nomod_dms", required=True, type=Path)
    ap.add_argument("--2a3-mod", dest="mod_2a3", required=True, type=Path)
    ap.add_argument("--2a3-nomod", dest="nomod_2a3", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8, help="cmuts MPI ranks / rf-count workers")
    ap.add_argument("--per-ref", type=Path, default=Path("outputs/per_reference.tsv"))
    ap.add_argument("--summary", type=Path, default=Path("outputs/summary.tsv"))
    args = ap.parse_args()

    args.per_ref.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    refs = load_fasta(args.fasta)

    # Build each dataset's reactivity from the alignments, keyed (dataset, cond).
    inputs = {
        "DMS": external.Inputs(args.fasta, args.mod_dms, args.nomod_dms),
        "2A3": external.Inputs(args.fasta, args.mod_2a3, args.nomod_2a3),
    }
    built = external.build_profiles(
        _datasets(args.threads),
        inputs,
        Path(tempfile.mkdtemp(prefix="tool-vs-structure-")),
        sm_dir=os.environ.get("SM2_DIR"),
    )
    profiles: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for dataset in DATASETS:
        if dataset not in built:
            continue
        for cond in CONDITIONS:
            if cond not in built[dataset]:
                continue
            ref_names, react, _ = external.read_profiles(built[dataset][cond])
            profiles[(dataset, cond)] = {ref_names[i]: react[i] for i in range(len(ref_names))}
    if not profiles:
        ap.error("no datasets built; need samtools, rf-count, and SM2_DIR for shapemapper2")

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
        for (dataset, cond), per_ref in profiles.items():
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
                    "dataset": dataset,
                    "condition": cond,
                    **res,
                }
            )
        if not any_scored:
            skipped["low_coverage"] += 1

    cols = ["reference", "pdb", "chain", "dataset", "condition", "n", "auc", "pearson", "spearman"]
    with open(args.per_ref, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    present = [d for d in DATASETS if any(r["dataset"] == d for r in rows)]

    # Apples-to-apples: average each dataset only over references that EVERY
    # present dataset scored. Coverage/quality filters drop different references
    # per tool, so the per-dataset valid sets differ; comparing on their
    # intersection removes that selection bias from the means. Computed per
    # condition (DMS and 2A3 score different references).
    scored: dict[tuple[str, str], set[str]] = {
        (d, cond): {r["reference"] for r in rows if r["dataset"] == d and r["condition"] == cond}
        for d in present
        for cond in CONDITIONS
    }
    common: dict[str, set[str]] = {
        cond: set.intersection(*(scored[(d, cond)] for d in present)) if present else set()
        for cond in CONDITIONS
    }

    def agg(metric: str, dataset: str, cond: str) -> tuple[int, float, float]:
        vals = np.array(
            [
                r[metric]
                for r in rows
                if r["dataset"] == dataset
                and r["condition"] == cond
                and r["reference"] in common[cond]
                and not np.isnan(r[metric])
            ]
        )
        if len(vals) == 0:
            return 0, np.nan, np.nan
        return len(vals), float(np.mean(vals)), float(np.median(vals))

    with open(args.summary, "w") as f:
        for cond in CONDITIONS:
            f.write(
                f"# {cond}: {len(common[cond])} references scored by all {len(present)} datasets\n"
            )
        f.write(f"condition\tmetric\tdataset\tn\tmean\tmedian\tmean_minus_{BASELINE}\n")
        for cond in CONDITIONS:
            for metric in ("auc", "pearson", "spearman"):
                _, base_mean, _ = agg(metric, BASELINE, cond)
                for dataset in present:
                    n, mean, median = agg(metric, dataset, cond)
                    delta = (
                        mean - base_mean if not (np.isnan(mean) or np.isnan(base_mean)) else np.nan
                    )
                    f.write(
                        f"{cond}\t{metric}\t{dataset}\t{n}\t{mean:.4f}\t{median:.4f}\t{delta:+.4f}\n"
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
