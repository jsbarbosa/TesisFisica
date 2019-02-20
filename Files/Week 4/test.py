import smbh
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt

r_vir = 9.79
G = 0.449338

r0 = 1e-4
t0 = np.pi / 4
p0 = np.pi / 4

v0 = 70 * 1.023 # to kpc / gyr
vt0 = np.pi / 4
vp0 = np.pi / 4

POS = smbh.sphericalToCartesian(r0, t0, p0)
SPEEDS = smbh.sphericalToCartesian(v0, vt0, vp0)

dt = 1e-6
smbh_m = 1
until = 1 # gyr
n_points = until // dt
filename = "Data/Results%.2f.dat"

smbh.setR_vir(r_vir)
t, pos, speeds, masses = smbh.run(POS, SPEEDS, smbh_m, dt, n_points, n_points // 10, filename%r_vir, delete_file = False)

# t, pos, speeds, masses = smbh.readFile(filename%r_vir)

accr = smbh.SMBHAccretion(pos, speeds)
dm = smbh.dynamicalFrictionDM(pos, speeds)

t *= 1000
r = smbh.magnitude(pos)
v = smbh.magnitude(speeds)
dm = smbh.magnitude(dm)

mass = 1e3
sigma = (0.5 * G * mass / r_vir) ** 0.5
x = v / (np.sqrt(2) * sigma)

grav = np.array([smbh.gravitationalForce(r_) for r_ in r])

rho_dm = np.array([smbh.darkMatterDensity(r_) for r_ in r])
rho_h = np.array([smbh.baryonicDensityHernquist(r_) for r_ in r])
m_dm = np.array([smbh.darkMatterMass(r_) for r_ in r])

erf_ = erf(x)
exp = 2 * x  * np.exp(-x**2) / np.sqrt(np.pi)

factor = erf_ - exp

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True)

# ax1.plot(t, m_dm, label = '$m_{dm}$')


ax1.plot(t, r/r_vir, '-o', ms = 1, c = 'k')
ax3.plot(t, rho_h, label = r'$\rho_h$')
ax3.plot(t, abs(grav), label = 'grav')

ax3.plot(t, dm, label = 'DF')
ax3.plot(t, rho_dm, label = r'$\rho_{dm}$')

# ax3.plot(t, v, label = 'v')
# ax3.plot(t, x, label = '$x$')

# ax1.plot(r/r_vir, abs(grav / dm))

ax2.plot(t, x, label = 'x')
ax2.plot(t, erf_, label = 'erf')
ax2.plot(t, exp, label = 'exp')
ax2.plot(t, factor, label = 'factor')


# ax1.plot(t, 1/v**3)
ax1.grid()
ax2.grid()
ax3.grid()
# ax1.set_yscale('log')
# ax2.set_yscale('log')
ax3.set_yscale('log')


ax1.legend()
ax2.legend()
ax3.legend()
# ax1.set_xlim(0, 10)

plt.show()
