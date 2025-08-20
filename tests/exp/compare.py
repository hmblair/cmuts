import sys
import xml.etree.ElementTree as ET
import numpy as np
import h5py
import matplotlib.pyplot as plt


CMUTS_FILE = sys.argv[1]
RF_FILE = sys.argv[2]
NAME = sys.argv[3]

TOLERANCE = 1E-2

tree = ET.parse(RF_FILE)
root = tree.getroot()

reactivity_values = []
for transcript in root.findall('transcript'):
    reactivity = transcript.find('reactivity')
    if reactivity is not None:
        reactivity_values = reactivity.text.split(',')

rf = np.array([x.strip() for x in reactivity_values], dtype=float)

with h5py.File(CMUTS_FILE, "r") as f:
    cmuts = f['onp/reactivity'][:]
    reads = f['onp/reads'][:]

if rf.ndim == 1:
    rf = rf[None, ...]

mae = np.mean(np.abs(cmuts - rf))
mae_indiv = np.mean(np.abs(cmuts - rf), axis=1)

print()
print("         TEST RESULTS")
print(" ─────────────────────────────")

if (mae < TOLERANCE):
    print("   STATUS: PASSED")
else:
    print("   STATUS: FAILED")

print(f"   MAE: {mae:.4f}")

ix = reads.argmax()

plt.plot(rf[ix], label="rf-count", color="blue", alpha=0.5, linewidth=1.0)
plt.plot(cmuts[ix], label="cmuts", color="red", alpha=0.5, linewidth=1.0)
plt.xlabel("Residue", fontsize=11)
plt.ylabel("Reactivity Rate", fontsize=11)
plt.title(f"MAE: {mae_indiv[ix]:.4f}", fontsize=12)
plt.legend()
plt.savefig("../figures/rf-vs-cmuts.png", dpi=300, bbox_inches="tight")
