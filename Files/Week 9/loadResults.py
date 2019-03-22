import sys
sys.path += [".."]

import smbh
import numpy as np
import matplotlib.pyplot as plt

positions = 1e-3 * np.ones(3) / (3 ** 0.5)
speeds = [70, 0, 0]
smbh_mass = 1
dt = 1e-6
d_q0 = [-1e3, 0, 0]
d_p0 = [0, 0, 0]

T = 1e-5
ls = np.logspace(0, 3, 10).astype(int)
coeffs = []

for l in ls:
    lya = smbh.lyapunov(positions, speeds, d_q0, d_p0, smbh_mass, l = l, T = T, triaxial = False)
    coeffs.append(lya)

plt.loglog(ls * T, coeffs)

print(coeffs)
plt.xlabel('lT')
plt.ylabel('L')

plt.show()
