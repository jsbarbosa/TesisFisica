import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# smbh.setTriaxialCoeffs(1, 0.99, 0.95)

a1, a2, a3 = smbh.getTriaxialCoeffs()

x_vir = smbh.R_VIR_z20
# y_vir = z_vir = x_vir
y_vir = (a2 / a1) * smbh.R_VIR_z20
z_vir = (a3 / a1) * smbh.R_VIR_z20

# poss = [[x_vir, 0, 0],
#         [0, y_vir, 0],
#         [0, 0, z_vir]]

coeffs = np.array([[1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])

r_0 = np.ones(3) / 3 ** 0.5
phi_0 = smbh.darkMatterPotential_triaxial(*r_0) \
        + smbh.stellarPotential_triaxial(*r_0) \
        + smbh.gasPotential_triaxial(*r_0)

rs = np.logspace(-6, 2)
pots = np.zeros((len(rs), 3))
for (i, r) in enumerate(rs):
    for (j, coeff) in enumerate(coeffs):
        pos = coeff * r
        pot = smbh.darkMatterPotential_triaxial(*pos) + \
                smbh.stellarPotential_triaxial(*pos) + \
                smbh.gasPotential_triaxial(*pos)
        pots[i, j] = pot

#     print(abs(2 * (pot - phi_0) ** 0.5))
#
# for pos in poss:
#     pos = smbh.R_VIR_z20
#     pot = smbh.darkMatterPotential(pos)
#     pot += smbh.stellarPotential(pos)
#     pot += smbh.gasPotential(pos)
#     print(abs(2 * pot) ** 0.5)

labels = ["x", "y", "z"]

for j in range(3):
    plt.plot(rs, pots[:, j], label = labels[j])

plt.xscale('log')
# plt.yscale('log')
plt.legend()
plt.show()


# n_random = 5000
#
# v_min = 30
# v_max = 85
#
# thetas = 0.5 * np.pi * np.random.random(n_random)
# phis = 0.5 * np.pi * np.random.random(n_random)
# vs = (v_max - v_min) * np.random.random(n_random) + v_min
#
# x, y, z = smbh.sphericalToCartesian(vs, thetas, phis)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection = '3d')
#
# ax.scatter(x, y, z, alpha = 0.5, s = 5)
#
# ax.set_xlim(0, v_max)
# ax.set_ylim(0, v_max)
# ax.set_zlim(0, v_max)
#
# plt.show()
