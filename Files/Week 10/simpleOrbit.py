import sys
sys.path += [".."]

import smbh

import matplotlib.pyplot as plt

results = smbh.Results("tOrbit.dat")
# x = results.times
# y = results.distance / results.R_VIR
x, y, z = results.positions.T

fig, ax = smbh.make3dPlot(results.positions)
fig.savefig('3dorbit.png', dpi = 300)

t = results.times * 1000
r = results.distance / results.R_VIR

results.setDistance(1, 0.5, 0.25)
m = results.distance / results.R_VIR

fig, ax = plt.subplots(1)

ax.plot(t, r, "-", label = '$r$')
ax.plot(t, m, '-', label = '$m$')

ax.set_xlabel('Time (Myr)')
ax.set_ylabel('Distance / $R_{vir}$')
ax.legend()
ax.grid()

fig.savefig('3dorbit_distances.png', dpi = 300)

plt.show()
