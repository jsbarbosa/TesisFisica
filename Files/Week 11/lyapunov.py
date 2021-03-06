import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt

pos = 1e-3 * np.ones(3) / 3 ** 0.5
dq = np.ones(3) * 1e-3

v = 25
theta = np.pi / 4
phi = np.pi / 4

x = v * np.sin(theta) * np.cos(phi)
y = v * np.sin(theta) * np.sin(phi)
z = v * np.cos(theta)

T = 1e-5

speeds = [x, y, z]
ls = np.linspace(1e1, 1e3).astype(int)
lya_s = np.zeros(len(ls))
lya_t = np.zeros(len(ls))

file_s = open('lyapunov_s.dat', 'w')
file_t = open('lyapunov_t.dat', 'w')

for (i, l) in enumerate(ls):
    print('(%d) l = %d' % (i, l))
    lya_s[i] = smbh.lyapunov(speeds, dq, T = T, l = l, triaxial = False)[0]
    lya_t[i] = smbh.lyapunov(speeds, dq, T = T, l = l, triaxial = True)[0]

    file_s.write("%d %f\n" % (l, lya_s[i]))
    file_t.write("%d %f\n" % (l, lya_t[i]))

file_s.close()
file_t.close()

# lya_s = np.genfromtxt("lyapunov_s.dat")[:, 1]
# lya_t = np.genfromtxt("lyapunov_t.dat")[:, 1]

fig, ax = plt.subplots(figsize = (5, 4.5))
ax.plot(ls * T, lya_s, label = 'Spherical')
ax.plot(ls * T, lya_t, label = 'Triaxial')

ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel('$l$ $T$ (Gyr)')
ax.set_ylabel(r'$\mathcal{L}$')

ax.grid()
ax.legend()

fig.tight_layout()

fig.savefig('lyapunov_.png', dpi = 300)
# fig.savefig('lyapunov_triaxial.png', dpi = 300)

plt.show()

# for l in ls:
#     a = smbh.lyapunov(pos, speeds, d_q0 = dq, T = 1e-5, l = l, triaxial = False)
    # print(l, a)
