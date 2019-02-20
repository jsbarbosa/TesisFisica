import os
import smbh
import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt

r0 = 0
t0 = np.pi / 4
p0 = np.pi / 4

v0 = 70 * 1.023 # to kpc / gyr
vt0 = np.pi / 4
vp0 = np.pi / 4

POS = np.array(smbh.sphericalToCartesian(r0, t0, p0))
SPEEDS = smbh.sphericalToCartesian(v0, vt0, vp0)

dt = 1e-6
smbh_m = 1
until = 0.06 # gyr
n_points = until // dt
filename = "Data/Results%.3f.dat"

n = 50
label_every = n // 10
# r_virs = np.linspace(0.47, 0.48, n)
r_virs = np.logspace(-2, 3, n)

colors = cm.jet(np.linspace(0, 1, n))

# for r_vir in r_virs:
r_vir = 9.79
smbh.setR_vir(r_vir)
smbh.run(POS, SPEEDS, smbh_m, dt, n_points, n_points // 100, filename%r_vir, delete_file = False)

r_virs = [r_vir]

fig, ax = plt.subplots()
for (i, r_vir) in enumerate(r_virs):
    points =  until // dt
    smbh.setR_vir(r_vir)
    t, pos, speeds, mass = smbh.readFile(filename% r_vir)
    r = np.linalg.norm(pos, axis = 1) / r_vir
    if i%label_every == 0:
        ax.plot(1000 * t, r, label = "$R_{vir} = %.4f$ kpc"%r_vir, c = colors[i])
    else:
        ax.plot(1000 * t, r, c = colors[i])

ax.set_xlabel('Time since recoil (Myr)')
ax.set_ylabel(r'$R / R_{vir}$')
# ax.set_yscale('log')
ax.set_ylim(0, 1)
ax.legend()
fig.savefig('recoils.png', dpi = 300)
plt.show()
