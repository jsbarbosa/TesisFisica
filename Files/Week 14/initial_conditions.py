import sys
sys.path += ["../"]

import smbh
import numpy as np
from time import time

files_name = "Results/%d.dat"

def do(i, speeds):
    # smbh.setTriaxialCoeffs(*axis)
    e_v = smbh.getEscapeSpeed()
    speeds = e_v * speeds
    t_s = ", ".join(["%.2f" % speeds[k] for k in range(3)])
    print('(%d) Started: speeds = (%s)' % (i, t_s))
    t0 = time()
    smbh.run(speeds, save_every = 100, filename = files_name % (i), read = False, pot_type = smbh.SYMMETRIC)
    # c = smbh.lyapunov(speeds)[0]
    # with open("lyapunov.txt", "a") as file: file.write("%d, %d, %e\n" % (i, j, c))
    print('(%d) Ended (%d s): speeds = (%s)' % (i, time() - t0, t_s))

if __name__ == '__main__':
    # axis = randomSemiaxis()
    # initial = randomInitialConditions().T
    initial = np.genfromtxt('initial.txt', skip_header = 1, usecols = [1, 2, 3], delimiter = ', ')
    args = [(i, vs) for (i, vs) in enumerate(initial)]
    smbh.run_multiproccess(do, args, 4)
