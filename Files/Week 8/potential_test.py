import sys
sys.path += [".."]

import smbh
import numpy as np
import matplotlib.pyplot as plt

def getPotential(func, x, y, z, gamma = 0.2):
    int_ = func(x, y, z, gamma)
    int_ = np.array([int_[i] for i in range(3)])
    return sum(int_ ** 2) ** 0.5

n = 100
xs = np.logspace(-4, 1, n)
y = 0
z = 0

distances = np.zeros(n)
errors = np.zeros((n, 4))

smbh.setTriaxalCoeffs(1, 1, 1)

for i in range(n):
    x = xs[i]
    r = (x * x + y * y + z * z) ** 0.5

    m_dm = smbh.darkMatterMass(r)
    m_s = smbh.stellarMassHernquist(r)
    m_g = smbh.gasMass(r)
    p_dm = m_dm * smbh.G / (r ** 2)
    p_s = m_s * smbh.G / (r ** 2)
    p_g = m_g * smbh.G / (r ** 2)

    p_t = (m_dm + m_s + m_g) * smbh.G / (r ** 2)

    int_dm = getPotential(smbh.triaxial_gravDM, x, y, z)
    int_s = getPotential(smbh.triaxial_gravS, x, y, z)
    int_g = getPotential(smbh.triaxial_gravG, x, y, z)

    int_t = int_dm + int_s + int_g

    e_dm = 100 * abs(p_dm - int_dm) / p_dm
    e_s = 100 * abs(p_s - int_s) / p_s
    e_g = 100 * abs(p_g - int_g) / p_g
    e_t = 100 * abs(p_t - int_t) / p_t

    distances[i] = r
    errors[i] = [e_dm, e_s, e_g, e_t]

fig, ax = plt.subplots(figsize = (5, 4.5))

distances = distances / smbh.R_VIR_z20

ax.plot(distances, errors[:, 0], label = 'NFW')
ax.plot(distances, errors[:, 1], label = 'Hernquist')
ax.plot(distances, errors[:, 2], label = 'Gas')
ax.plot(distances, errors[:, 3], label = 'Cumulative')

ax.legend()

ax.grid()

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("$R/R_{vir}$")
ax.set_ylabel("Error (%)")

# ax.set_ylim(0, 1)
fig.tight_layout()

fig.savefig("error.png", dpi = 300)

plt.show()
