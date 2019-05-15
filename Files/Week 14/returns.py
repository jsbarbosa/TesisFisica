import sys
sys.path += ["../"]

import smbh
import numpy as np

fileFormat = "../Week 13/Results/%d_%d.dat"

n = 30
m = 1000

masses = np.zeros((n, m))

for i in range(n):
    for j in range(m):
        file = fileFormat % (i, j)
        results = smbh.Results(file)
        masses[i, j] = results.getMassAt(0.512)
        del results
        print(i, j, masses[i, j])


np.savetxt("masses_at_time.txt", masses)
