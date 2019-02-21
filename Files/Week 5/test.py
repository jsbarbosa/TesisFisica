import sys
sys.path += [".."]

import smbh
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

r_vir = 9.79
r0 = 0
v0 = 70 * 1.023 # to kpc / gyr

POS = smbh.sphericalToCartesian(r0, np.pi / 4, np.pi / 4)
SPEEDS = smbh.sphericalToCartesian(v0, np.pi / 4, np.pi / 4)

dt = 1e-6
smbh_m = 1
until = 0.065 # gyr
n_points = until // dt
filename = "Results.dat"

smbh.setR_vir(r_vir)
smbh.printConstants()
t, pos, speeds, masses = smbh.run(POS, SPEEDS, smbh_m, dt, n_points, n_points // 10, filename, delete_file = False)

# t, pos, speeds, masses = smbh.readFile(filename%r_vir)
smbh.plotProperties(t, pos, speeds, masses, r_vir)

plt.show()
