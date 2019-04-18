import sys
sys.path += [".."]

import smbh
import numpy as np

import matplotlib.pyplot as plt

smbh.setTriaxialCoeffs(1, 0.99, 0.95)

results_x = smbh.run([30, 0, 0], save_every = 10, filename = "x.dat")
print("X done")

results_y = smbh.run([0, 30.3, 0], save_every = 10, filename = "y.dat")
print("Y done")

results_z = smbh.run([0, 0, 31.46], save_every = 10, filename = "z.dat")
print("Z done")

# results_xyz = smbh.run([70 / 3 ** 0.5, 70 / 3 ** 0.5, 70 / 3 ** 0.5], save_every = 10, filename = "xyz.dat")
# print("Z done")


# results_x = smbh.Results("x.dat")
# results_y = smbh.Results("y.dat")
# results_z = smbh.Results("z.dat")
# results_xyz = smbh.Results("xyz.dat")

results = [results_z, results_y, results_x]
labels = [r"$\vec{v}_0 = 70 \hat{%s}$"%i for i in ["k", "j", "i"]]

# results = [results_z, results_xyz, results_y, results_x]
# labels = [r"$\vec{v}_0 = 70 \hat{%s}$"%i for i in ["k", "ijk", "j", "i"]]


for result in results:
    result.setDistance(1, 0.99, 0.95)

fig, ax = smbh.plotOrbits(results, figsize = (5, 4.5))

ax.set_ylabel(r"$m / R_{vir}$")
ax.legend(labels)

fig.tight_layout()
fig.savefig("orthogonal_triaxial.png", dpi = 300)

plt.show()
