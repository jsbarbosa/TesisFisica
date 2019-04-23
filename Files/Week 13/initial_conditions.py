import sys
sys.path += ["../"]

import smbh
import numpy as np
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

files_name = "Results/%d_%d.dat"

def randomInitialConditions(n_random = 500):
    v_min = 0.45
    v_max = 0.87

    thetas = 0.5 * np.pi * np.random.random(n_random)
    phis = 0.5 * np.pi * np.random.random(n_random)
    vs = (v_max - v_min) * np.random.random(n_random) + v_min

    xs = smbh.sphericalToCartesian(vs, thetas, phis)
    return np.array(xs)

def randomSemiaxis(n_random = 500):
    a_2 = np.random.random(n_random)
    a_3 = a_2 * np.random.random(n_random)

    axis = np.ones((n_random, 3))
    axis[:, 1] = a_2
    axis[:, 2] = a_3
    return axis

def do(i, j, axis, speeds):
    smbh.setTriaxialCoeffs(*axis)
    e_v = smbh.getEscapeSpeedOrthogonalTriaxial()
    speeds = e_v * speeds
    t_a = ", ".join(["%.2f" % axis[k] for k in range(3)])
    t_s = ", ".join(["%.2f" % speeds[k] for k in range(3)])
    print('(%d, %d) Started: axis = (%s), speeds = (%s)' % (i, j, t_a, t_s))
    t0 = time()
    smbh.run(speeds, save_every = 100, filename = files_name % (i, j), read = False, pot_type = smbh.TRIAXIAL)
    c = smbh.lyapunov(speeds)[0]
    with open("lyapunov.txt", "a") as file: file.write("%d, %d, %e\n" % (i, j, c))
    print('(%d, %d) Ended (%d s): axis = (%s), speeds = (%s)' % (i, j, time() - t0, t_a, t_s))

if __name__ == '__main__':
    axis = randomSemiaxis()
    initial = randomInitialConditions().T
    args = [(i, j, ax, vs) for (i, ax) in enumerate(axis) for (j, vs) in enumerate(initial)]

    with open("axis.txt", "w") as file:
        file.write("index, a_1, a_2, a_3\n")
        for (i, ax) in enumerate(axis):
            ax = ", ".join(["%e" % a for a in ax])
            txt = "%d, %s\n" % (i, ax)
            file.write(txt)

    with open("initial.txt", "w") as file:
        file.write("index, rv_x, r_vy, r_vz\n")
        for (i, v) in enumerate(initial):
            ax = ", ".join(["%e" % a for a in v])
            txt = "%d, %s\n" % (i, ax)
            file.write(txt)

    with open("lyapunov.txt", "w") as file:
        file.write("i, j, lyapunov\n")
    smbh.run_multiproccess(do, args, 4)
