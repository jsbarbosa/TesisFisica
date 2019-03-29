import sys
sys.path += [".."]

import smbh
import numpy as np

from glob import glob

format = "Symmetric/return_%d_%d.dat"
#format = "Triaxial/return_%d_%d.dat"


format_1 = format.replace("_%d_", "_0_").replace("_%d", "*")
format_2 = format.replace("%d", "*")
n = len(glob(format_1))
m = len(glob(format_2)) // n

print(n, m)

# results = np.zeros((m, n), dtype = object)
speeds = np.zeros(m)
results = np.zeros((m, n, 3))

for i in range(0, m):
    for j in range(n):
        r = smbh.Results(format % (i, j), header_only = True)
        try:
            t = r.RETURN_TIME
            m_ = r.RETURN_MASS
        except AttributeError:
            print(i, j)
            t, m_ = -1, -1
        f = r.STELLAR_FRACTION
        v = r.getInitialSpeed()
        results[i, j] = [f, t, m_]
        del r
        
    speeds[i] = v

for i in range(m):
    txt = "fs\tTr(Gyr)\tm(1e5M0)\n"
    txt_ = ["%f\t%f\t%f" % tuple(data) for data in results[i]]
    txt += "\n".join(txt_)
    with open("s%d_resume.dat" % speeds[i], "w") as f: f.write(txt)
    f.close()
