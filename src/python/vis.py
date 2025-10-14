import numpy as np
import matplotlib.pyplot as plt
import subprocess
from Bio import Align, PDB

FIGURES = "figures"


def _seq_from_fasta(file: str, n: int = 0) -> str:

    seq = ""
    count = 0

    with open(file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
                continue
            if count == n + 1:
                seq += line.strip()

    return seq.upper().replace("U", "T")


def _seq_from_cif(filename: str, chain_id: str | None = None) -> str:

    if chain_id is None:
        chain_id = 'A'

    if filename.endswith('.cif'):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)

    rna_dict = {
        'A': 'A', 'ADE': 'A', 
        'U': 'U', 'URA': 'U', 'URI': 'U',
        'G': 'G', 'GUA': 'G', 'GTP': 'G',
        'C': 'C', 'CYT': 'C',
        'PSU': 'U',
        '5MU': 'U',
        'H2U': 'U',
        'M2G': 'G',
        'M7G': 'G',
        '1MA': 'A',
        '5MC': 'C',
    }

    try:
        structure = parser.get_structure("rna", filename)

        # Get the specified chain
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                sequence = ""

                for residue in chain:
                    resname = residue.get_resname().strip()

                    # Check if it's an RNA nucleotide
                    if resname in rna_dict:
                        sequence += rna_dict[resname]
                    elif PDB.is_aa(residue):
                        # Skip amino acids
                        continue
                    else:
                        # Unknown nucleotide, use 'N'
                        sequence += 'N'

                return sequence.replace('U', 'T')

        print(f"Chain {chain_id} not found in structure")
        return ""

    except Exception as e:
        print(f"Error parsing file: {e}")
        return ""


def _seq_align(seq1: str, seq2: str) -> tuple[str, str]:

    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]

    return best_alignment[0], best_alignment[1]


def _data_aln(data: np.ndarray, aln1: str, aln2: str):
    """
    Map data values to align with seq2, using zeros for gaps.
    """

    if data.shape[0] != len([c for c in aln1 if c != '-']):
        raise ValueError("Data length must match number of non-gap characters in aln1")

    mapped_data = []
    data_idx = 0

    for i, (char1, char2) in enumerate(zip(aln1, aln2)):
        if char2 == '-':
            # Gap in seq2, skip this position entirely
            if char1 != '-':
                data_idx += 1  # Still advance data index for seq1
            continue

        if char1 == '-':
            # Gap in seq1, insert zero for seq2
            mapped_data.append(0)
        else:
            # Both sequences have characters, use data value
            mapped_data.append(data[data_idx])
            data_idx += 1

    return np.array(mapped_data)


def _plot_single_profile(
    reactivity: np.ndarray,
    name: str,
    dir: str = FIGURES
) -> None:

    prefix = f"{name}-" if name else ""
    x = range(len(reactivity))

    plt.grid(axis="y", alpha=0.5)
    plt.fill_between(x, reactivity, alpha=0.5, color=plt.cm.RdPu(0.3))
    plt.plot(reactivity, color=plt.cm.RdPu(0.8), linewidth=1)

    plt.xlabel("Residue", fontsize=14)
    plt.ylabel("Reactivity", fontsize=14)
    if name:
        plt.title(f"Profile of {name}", fontsize=14)
    plt.tick_params(axis='both', labelsize=13)

    plt.savefig(f"{dir}/{prefix}profile.png", dpi=300, bbox_inches="tight")
    plt.close()


def _to_defattr(
    values: np.ndarray,
    out: str,
    start_resnum: int = 1,
    attr_name: str = "value",
    pad5: int = 0,
    pad3: int = 0,
    chain: str | None = None
):
    """
    Convert numpy array to ChimeraX defattr format.
    """

    values = np.concatenate([np.zeros(pad5), values, np.zeros(pad3)])

    with open(out, 'w') as f:

        f.write(f"attribute: {attr_name}\n")
        f.write("recipient: residues\n\n")

        for i, value in enumerate(values):
            resnum = start_resnum + i
            if chain is not None:
                f.write(f"\t/{chain}:{resnum}\t{value}\n")
            else:
                f.write(f"\t:{resnum}\t{value}\n")


def _color_by_defattr(bin: str, cif: str, defattr: str, color: str) -> None:
    """
    Color a .cif file with the given .defattr file.
    """

    chmx_cmd = (
        f"open {cif}; " +
        "color grey; " +
        "graphics quality 5; " +
        "renumber start 1 relative false; " +
        f"open {defattr}; " +
        f"color byattribute value palette white:{color} range 0,1; " +
        "hide cartoons; " +
        "nucleotides atoms; " +
        "style sphere; " +
        "lighting soft; " +
        "lighting ambientIntensity 1.3"
    )

    cmd = [bin, "--cmd", chmx_cmd]
    subprocess.run(cmd)


def _color_by_reactivity(
    bin: str,
    cif: str,
    reactivity: np.ndarray,
    chain: str | None,
    color: str,
) -> None:
    """
    Color a .cif file with the given reactivity.
    """

    defattr = ".cmuts-visualize.defattr"
    _to_defattr(reactivity, defattr, chain=chain)
    _color_by_defattr(bin, cif, defattr, color)


def _color_structure(
    bin: str,
    cif: str,
    reactivity: np.ndarray,
    seq1: str,
    seq2: str,
    chain: str | None,
    color: str,
) -> None:

    aln1, aln2 = _seq_align(seq1, seq2)
    aln_data = _data_aln(reactivity, aln1, aln2)
    _color_by_reactivity(bin, cif, aln_data, chain, color)
