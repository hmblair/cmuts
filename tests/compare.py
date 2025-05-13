import sys
import numpy as np
import h5py


def _open_rf(filename: str) -> np.ndarray:

    with open(filename, 'r') as file:
        content = file.read()
    blocks = content.strip().split('\n\n')

    data = []
    labels = []
    for block in blocks:
        lines = block.split('\n')
        labels.append(lines[0])
        block_data = [
            list(map(int, line.split('\t')[1].split(',')))
            for line in lines[1:]
        ]
        data.append(block_data)

    return np.array(data)


CMUTS_FILE = sys.argv[1]
EXPECTED_FILE = sys.argv[2]
RF_FILE = "rf_count/raw_counts/aln.txt"
TOLERANCE = 1E-6

with h5py.File(CMUTS_FILE, "r") as cmuts_f, h5py.File(EXPECTED_FILE, "r") as exp_f:

    cmuts = cmuts_f['aln'][:]
    exp = exp_f['aln'][:]
    diff_with_exp = cmuts - exp

    cmuts_matches = (
        cmuts[..., 0, 0] +
        cmuts[..., 1, 1] +
        cmuts[..., 2, 2] +
        cmuts[..., 3, 3]
    )
    cmuts_cov = cmuts[..., :-1].sum((2, 3))
    cmuts_muts = cmuts.sum((2, 3)) - cmuts_matches

    exp_matches = (
        exp[..., 0, 0] +
        exp[..., 1, 1] +
        exp[..., 2, 2] +
        exp[..., 3, 3]
    )
    exp_cov = exp[..., :-1].sum((2, 3))
    exp_muts = exp.sum((2, 3)) - exp_matches

    if (np.abs(diff_with_exp).sum() < TOLERANCE):
        if (exp_cov.sum() < TOLERANCE):
            print("   STATUS: (TRIVIALLY) PASSED")
        else:
            print("   STATUS: PASSED")
        print()
        sys.exit(0)
    else:
        print("   STATUS: FAILED")
        print("   MEAN DIFFERENCE: ", np.abs(diff_with_exp).mean())
        print("   MAX DIFFERENCE:  ", np.abs(diff_with_exp).max())
        print(f"   EXP COV:   {exp_cov.sum().item()}")
        print(f"   EXP MUT:   {exp_muts.sum().item()}")
        print(f"   CMUTS COV: {cmuts_cov.sum().item()}")
        print(f"   CMUTS MUT: {cmuts_muts.sum().item()}")
        print("   MUTATION TABLE:")
        print(diff_with_exp.sum((0, 1)))
        print()
        sys.exit(1)
