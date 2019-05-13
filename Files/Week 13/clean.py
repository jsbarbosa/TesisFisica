import numpy as np

AXIS = 'axis.txt'
MASSES = 'return_masses.txt'
TIMES = 'return_times.txt'
LYAPUNOV = 'lyapunov_exponents.txt'

axis = np.genfromtxt(AXIS, skip_header = 1, delimiter = ', ')
masses = np.genfromtxt(MASSES)
times = np.genfromtxt(TIMES)
lyapunov = np.genfromtxt(LYAPUNOV)

ignore = []

for i, row in enumerate(masses):
    pos = np.where(row > 1e1)[0]
    smaller = np.where(row < 1)[0]
    masses[i, smaller] = 1 + np.random.random(len(smaller))
    if len(pos): ignore.append(i)
    
save = [i for i in range(len(masses)) if i not in ignore]
print(save)

arrays = [axis, masses, times, lyapunov]
files = [AXIS, MASSES, TIMES, LYAPUNOV]

for (array, file) in zip(arrays, files):
    array = array[save]
    np.savetxt(file, array)
    
