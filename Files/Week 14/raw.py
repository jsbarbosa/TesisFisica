import sys
sys.path += ['..']

import smbh
import numpy as np

n = 1000

data = np.zeros((n, 2))

for i in range(n):
    r = smbh.Results("Results/%d.dat" % i, header_only = True)
    data[i] = r.RETURN_TIME, r.RETURN_MASS

np.savetxt("results.txt", data)
