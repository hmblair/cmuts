import sys
import numpy as np
import h5py


CMUTS_FILE = sys.argv[1]
EXPECTED_FILE = sys.argv[2]
TOLERANCE = 1E-6

with h5py.File(CMUTS_FILE, "r") as cmuts_f, h5py.File(EXPECTED_FILE, "r") as exp_f:
    cmuts = cmuts_f['aln'][:]
    exp = exp_f['aln'][:]
    diff = np.abs(cmuts - exp)

    if (diff.sum() < TOLERANCE):
        print("   STATUS: PASSED")
        print()
    else:
        print("   STATUS: FAILED")
        print("   MEAN DIFFERENCE: ", diff.mean())
        print("   MAX DIFFERENCE:  ", diff.max())
        print()
