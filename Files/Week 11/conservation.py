import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt

# speeds = [70, 0, 0]
#
# results_s = smbh.run(speeds, filename = 's_conservative.dat', end_time = 6, save_every = 10, pot_type = smbh.C_SYMMETRIC)
#
# print('Symmetric')
#
# smbh.setTriaxialCoeffs(1, 1, 1)
# results_t = smbh.run(speeds, filename = 't_conservative.dat', end_time = 6, save_every = 10, pot_type = smbh.C_TRIAXIAL)
# print('Triaxial')

# results_s = smbh.Results('s_conservative.dat')
# results_t = smbh.Results('t_conservative.dat')

# results = [results_s, results_t]

files = ['s_energy.dat', 't_energy.dat']

# for i in range(2):
#     t = results[i].times
#     v = results[i].speed
#     r = results[i].distance
#
#     pot = smbh.darkMatterPotential(r) + smbh.stellarPotential(r) + smbh.gasPotential(r)
#
#     k = 0.5 * v ** 2
#     p = pot
#
#     e = k + p
#
#     text = ["Time (Gyr) Kinetic Potential Energy"] + ["%e %e %e %e" % (t[i], k[i], p[i], e[i]) for i in range(len(k))]
#
#     text = '\n'.join(text)
#
#     with open(files[i], 'w') as file: file.write(text)

# fig1, axes = plt.subplots(1, 2, sharey = True, figsize = (8, 4.5))
fig2, ax = plt.subplots(figsize = (8, 4.5))

labels = ['Symmetric', 'Triaxial ($a_1 = a_2 = a_3 = 1$)']

for i in range(2):
    data = np.genfromtxt(files[i], skip_header = 1)
    t, k, p, e = data.T

    max_ = e.max()
    min_ = e.min()
    mean = e.mean()

    print("Max: %.2f, Min: %.2f, Fluctuacion: %.3f %%" % (max_, min_, 100 * abs((max_ - min_) / mean)))

    # axes[i].plot(t, k, label = 'K')
    # axes[i].plot(t, p, label = 'V')
    # axes[i].plot(t, e, label = 'K + V')

    ax.plot(t, 100 * (e - mean) / mean, label = labels[i])

    # axes[i].set_xlabel('Time (Gyr)')
    # axes[i].legend()
    # axes[i].grid()

# axes[0].set_ylabel('Energy ($10^5M_\odot$ kpc$^2$Gry$^{-2}$)')

ax.set_xlabel('Time (Gyr)')
ax.set_ylabel('Energy variation (%)')

ax.legend()
ax.grid()

# fig1.tight_layout()
fig2.tight_layout()

# fig1.subplots_adjust(wspace = 0.05)

# fig1.savefig('Energies.png', dpi = 300)
fig2.savefig('Comparison.png', dpi = 300)
plt.show()
