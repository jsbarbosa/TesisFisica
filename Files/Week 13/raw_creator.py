import sys
sys.path += ["../"]

import smbh
import numpy as np

n = 30
m = 1000

files = [["Results/%d_%d.dat" % (i, j) for j in range(m)] for i in range(n)]

results = np.zeros((n, m), dtype = object)
for i in range(n):
    for j in range(m):
        results[i, j] = smbh.Results(files[i][j], header_only = True)
    print(i)

return_times = np.zeros_like(results).astype(float)
return_masses = np.zeros_like(return_times)

for i in range(results.shape[0]):
    for j in range(results.shape[1]):
        return_times[i, j] = results[i, j].RETURN_TIME
        return_masses[i, j] = results[i, j].RETURN_MASS


lyapunov = np.genfromtxt("lyapunov.txt", skip_header = 1, delimiter = ', ')
lyapunov_ = np.zeros((n, m))

for l in lyapunov:
    i, j, l = l
    i = int(i)
    j = int(j)
    try:
        lyapunov_[i, j] = l
    except IndexError:
        print(i, j, 'ignored')
np.savetxt('return_times.txt', return_times)
np.savetxt('return_masses.txt', return_masses)
np.savetxt('lyapunov_exponents.txt', lyapunov_)
