import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt

pos = 1e-3 * np.ones(3) / 3 ** 0.5

v = 60
theta = np.pi / 4
phi = np.pi / 4

x = v * np.sin(theta) * np.cos(phi)
y = v * np.sin(theta) * np.sin(phi)
z = v * np.cos(theta)

speeds = [x, y, z]
ls = np.logspace(1, 3, 10).astype(int)
lya_s = np.zeros((len(ls), 7))
lya_t = np.zeros((len(ls), 7))

T = 1e-5
dq = 1e-8

with open('lyapunov_s.dat', 'w') as file_s: pass
with open('lyapunov_t.dat', 'w') as file_t: pass
for (i, l) in enumerate(ls):
    lya_s[i] = l, *smbh.lyapunov(pos, speeds, dq, T = T, l = l, triaxial = False)
    lya_t[i] = l, *smbh.lyapunov(pos, speeds, dq, T = T, l = l, triaxial = True)

    with open('lyapunov_s.dat', 'a') as file_s:
        txt = "%d %f %f %f %f %f %f\n" % tuple(lya_s[i])
        file_s.write(txt)
    with open('lyapunov_t.dat', 'a') as file_t:
        txt = "%d %f %f %f %f %f %f\n" % tuple(lya_t[i])
        file_t.write(txt)

# lya_s = np.genfromtxt('lyapunov_s.dat')
# lya_t = np.genfromtxt('lyapunov_t.dat')

fig, ax = plt.subplots(figsize = (8, 4.5))
ax.plot(lya_s[:, 0] * T, lya_s[:, 1], label = "Spherical")
ax.plot(lya_t[:, 0] * T, lya_t[:, 1], label = "Triaxial")
# fig, axes = plt.subplots(1, 3, figsize = (8, 4.5))
# for i in range(3):
#     axes[i].plot(lya_s[:, 0] * T, lya_s[:, i * 2 + 1], label = 'Spherical +')
#     axes[i].plot(lya_t[:, 0] * T, lya_t[:, i * 2 + 1], label = 'Triaxial +')
#     axes[i].plot(lya_s[:, 0] * T, lya_s[:, i * 2 + 2], '--o', label = 'Spherical -')
#     axes[i].plot(lya_t[:, 0] * T, lya_t[:, i * 2 + 2], '--o', label = 'Triaxial -')
#
#     axes[i].set_xlabel('$lT$')
#     # axes[i].set_yscale('log')
#     axes[i].set_xscale('log')
# print(lya_s[-1, 0] * T)

# axes[0].set_ylabel(r'$\vec{x}$')

# ax.grid()
# axes[0].legend()

ax.set_ylabel(r'$\mathcal{L}$')
ax.set_xlabel('$l$ $T$ (Gyr)')

ax.set_xscale('log')
ax.legend()
ax.grid()

fig.tight_layout()
fig.savefig('lyapunov.png', dpi = 300)

plt.show()

# results = smbh.run(speeds, pot_type = smbh.SYMMETRIC)
# print(results.RETURN_TIME)
