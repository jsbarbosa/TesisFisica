import sys
sys.path += [".."]

import smbh
import numpy as np
import matplotlib.pyplot as plt

def getPotential(func, x, y, z):
    int_ = func(x, y, z)
    int_ = np.array([int_[i] for i in range(3)])
    return sum(int_ ** 2) ** 0.5

n = 100
xs = np.logspace(-3, 5, n)
y = 0
z = 0

distances = np.zeros((n, 3))
errors = np.zeros((n, 3))

smbh.setTriaxalCoeffs(1, 1, 1)

for i in range(n):
    x = xs[i]
    r = (x * x + y * y + z * z) ** 0.5

    m_dm = smbh.darkMatterMass(r)
    m_s = smbh.stellarMassHernquist(r)
    p_dm = m_dm * smbh.G / (r ** 2)
    p_s = m_s * smbh.G / (r ** 2)

    int_dm = getPotential(smbh.triaxial_gravDM, x, y, z)
    int_s = getPotential(smbh.triaxial_gravS, x, y, z)

    e_dm = 100 * abs(p_dm - int_dm) / p_dm
    e_s = 100 * abs(p_s - int_s) / p_s

    distances[i] = r
    errors[i] = [e_dm, e_s, 0]

fig, ax = plt.subplots(figsize = (8, 4.5))

distances = distances[:, 0] / smbh.R_VIR_z20

ax.plot(distances, errors[:, 0], label = 'NFW')
ax.plot(distances, errors[:, 1], label = 'Hernquist')

ax.legend()

ax.grid()

# ax.set_ylim(0, 1)
ax.set_xscale('log')
ax.set_xlabel("$R/R_{vir}$")
ax.set_ylabel("Error (%)")

fig.savefig("error.png")

plt.show()
