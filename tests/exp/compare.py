import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


CMUTS_FILE = sys.argv[1]
RF_MUTS_FILE = sys.argv[2]
RF_COV_FILE = sys.argv[3]

TOLERANCE = 1E-2

with h5py.File(CMUTS_FILE, "r") as f:

    cmuts = f["../data/pk50/aln"][:]

cmuts_matches = (
    cmuts[..., 0, 0] +
    cmuts[..., 1, 1] +
    cmuts[..., 2, 2] +
    cmuts[..., 3, 3]
)
cmuts_cov = cmuts[..., :-1].sum((2, 3))
cmuts_muts = cmuts.sum((2, 3)) - cmuts_matches

rf_muts = np.loadtxt(RF_MUTS_FILE, delimiter=",")
rf_cov = np.loadtxt(RF_COV_FILE, delimiter=",")

cmuts_profile = np.nan_to_num(cmuts_muts / cmuts_cov, 0.0)
rf_profile = np.nan_to_num(rf_muts / rf_cov, 0.0)

mae = np.mean(np.abs(cmuts_profile - rf_profile))

print()
print("         TEST RESULTS")
print(" ─────────────────────────────")

if (mae < TOLERANCE):
    print("   STATUS: PASSED")
else:
    print("   STATUS: FAILED")

print(f"   MAE: {mae:.4f}")

ix = cmuts_cov.max(-1).argmax()

plt.plot(rf_profile[ix], label="rf-count", color="blue", alpha=0.5, linewidth=1.0)
plt.plot(cmuts_profile[ix], label="cmuts", color="red", alpha=0.5, linewidth=1.0)
plt.xlabel("Residue", fontsize=11)
plt.ylabel("Reactivity Rate", fontsize=11)
plt.title(f"MAE: {mae:.4f}", fontsize=12)
plt.legend()
plt.savefig("../figures/rf-vs-cmuts.png", dpi=300, bbox_inches="tight")
