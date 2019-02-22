import sys
sys.path += [".."]

import smbh
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0 = 67.8, Om0 = 0.309, Ode0 = 0.73)

# z = np.linspace(0, 30)
z = 20
H = cosmo.H(z).value
# zs = np.linspace(30, 0)
# ts = cosmo.age(zs)
#
# plt.plot(ts, zs)
# plt.show()

# a, b, c = np.polyfit(z, H, 2)
# H_ = a * z ** 2 + b * z + c
# plt.plot(z, H)
# plt.plot(z, H_, '--')
# plt.show()


H = smbh.HubbleToNaturalUnits(H)
r_vir = smbh.calculateR_vir(smbh.G, H)
r0 = 0
v0 = 90 * 1.023 # to kpc / gyr

POS = smbh.sphericalToCartesian(r0, np.pi / 4, np.pi / 4)
SPEEDS = smbh.sphericalToCartesian(v0, np.pi / 4, np.pi / 4)

dt = 1e-6
smbh_m = 1
until = 0.057 # gyr
n_points = until // dt
filename = "Results.dat"

smbh.setR_vir(r_vir)
smbh.printConstants()
t, pos, speeds, masses = smbh.run(POS, SPEEDS, smbh_m, dt, n_points, n_points // 10, filename, delete_file = False)

# t, pos, speeds, masses = smbh.readFile(filename%r_vir)
smbh.plotProperties(t, pos, speeds, masses, r_vir)
plt.show()

r = smbh.magnitude(pos)
d = smbh.darkMatterDensity(r)
m = smbh.darkMatterMass(r)

smbh.plotDensityMass(r, d, m)
plt.show()
