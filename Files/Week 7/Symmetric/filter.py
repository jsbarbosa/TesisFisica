import sys
sys.path += ["../.."]

import smbh
import numpy as np

from glob import glob

format = "Data/return_%d_%d.dat"

format_1 = format.replace("_%d_", "_0_").replace("_%d", "*")
format_2 = format.replace("%d", "*")
n = len(glob(format_1))
m = len(glob(format_2)) // n

results = np.zeros((m, n), dtype = object)

for i in range(m):
    for j in range(n):
        results[i, j] = smbh.Results(format%(i, j))

for i in range(m):
    v = results[i, 0].speeds[0, 0]
    fs = [r.STELLAR_FRACTION for r in results[i]]
    ts = [r.getReturnTime() for r in results[i]]
    ms = [r.masses[r.return_index] for r in results[i]]
    txt = "fs\tTr(Gyr)\tm(1e5M0)\n"
    txt_ = ["%f\t%f\t%f" % data for data in zip(fs, ts, ms)]
    txt += "\n".join(txt_)
    with open("s%d_resume.dat"%v, "w") as f: f.write(txt)
