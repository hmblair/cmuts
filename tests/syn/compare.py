import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


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

RF_MUTS_FILE = "rf_muts.txt"
RF_COV_FILE = "rf_cov.txt"

TOLERANCE = 1E-2

with h5py.File(CMUTS_FILE, "r") as cmuts_f, h5py.File(EXPECTED_FILE, "r") as exp_f:

    cmuts = cmuts_f['aln'][:]
    exp = exp_f['aln'][:]
    diff_with_exp = (cmuts - exp)

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

    rf_muts = np.loadtxt(RF_MUTS_FILE, delimiter=",")
    rf_cov = np.loadtxt(RF_COV_FILE, delimiter=",")

    cmuts_profile = np.nan_to_num(cmuts_muts / cmuts_cov, 0.0)
    rf_profile = np.nan_to_num(rf_muts / rf_cov, 0.0)
    exp_profile = np.nan_to_num(exp_muts / exp_cov, 0.0)

    cmuts_profile[cmuts_profile > 1] = 1
    rf_profile[rf_profile > 1] = 1

    exp_mae = np.mean(np.abs(cmuts_profile - exp_profile))
    rf_mae = np.mean(np.abs(cmuts_profile - rf_profile))

    print()
    print("   TEST RESULTS (SYNTHETIC)")
    print(" ─────────────────────────────")

    if (exp_mae < TOLERANCE):
        if (exp_cov.sum() < TOLERANCE):
            print("   STATUS: (TRIVIALLY) PASSED")
        else:
            print("   STATUS: PASSED")
    else:
        print("   STATUS: FAILED")

    print(f"   MAE: {exp_mae}")

    print()
    print("    TEST RESULTS (RF-COUNT)")
    print(" ─────────────────────────────")

    if (rf_mae < TOLERANCE):
        if (exp_cov.sum() < TOLERANCE):
            print("   STATUS: (TRIVIALLY) PASSED")
        else:
            print("   STATUS: PASSED")
    else:
        print("   STATUS: FAILED")

    print(f"   MAE: {rf_mae}")
    
    plt.plot(rf_profile[0], label="rf-count", color="blue")
    plt.plot(cmuts_profile[0], label="cmuts", color="red")
    plt.legend()
    plt.show()
