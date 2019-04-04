import sys
sys.path += [".."]

import smbh

import matplotlib.pyplot as plt

colors = smbh.getColors(2)
r0 = smbh.Results("tOrbit0.dat")
fig, ax = smbh.make3dPlot(r0.positions, color = colors[0], alpha = 0.8)

r1 = smbh.Results("tOrbit1.dat")
fig, ax = smbh.make3dPlot(r1.positions, fig_axes = (fig, ax), color = colors[1], alpha = 0.8)

fig.tight_layout()
fig.savefig('3dorbit.png', dpi = 300)

fig, ax = plt.subplots(1)
for (i, results) in enumerate([r0, r1]):
    t = results.times * 1000
    r = results.distance / results.R_VIR

    results.setDistance(1, 0.5, 0.25)
    m = results.distance / results.R_VIR

    ax.plot(t, r, '-', lw = 1, color = colors[i])

ax.set_xlabel('Time (Myr)')
ax.set_ylabel('Distance / $R_{vir}$')
ax.grid()

fig.tight_layout()
fig.savefig('3dorbit_distances.png', dpi = 300)

plt.show()
