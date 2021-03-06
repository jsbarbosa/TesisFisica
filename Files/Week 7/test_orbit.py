import sys
sys.path += [".."]

import smbh
import numpy as np

import matplotlib.pyplot as plt

smbh.setTriaxalCoeffs(1, 1, 1)

results_0 = smbh.run([70, 0, 0], save_every = 10, filename = "symmetric.dat", triaxial = False)
print("Symmetric done")

results_1 = smbh.run([70, 0, 0], save_every = 10, filename = "triaxial.dat")
print("Triaxial done")
# results_0 = smbh.Results("symmetric.dat")
# results_1 = smbh.Results("triaxial.dat")

fig, ax = plt.subplots(figsize = (5, 4.5))

for result in [results_0, results_1]:
    r = result.distance / result.R_VIR
    t = result.times * 1000

    ax.plot(t, r, lw = 1)

ax.grid()
# smbh.plotProperties(results)
ax.set_xlabel("Time (Myr)")
ax.set_ylabel(r"$R / R_{vir}$")

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.legend(["Symmetric", "Triaxial $(a_1 = a_2 = a_3 = 1)$"])

fig.tight_layout()
fig.savefig("symmetric_triaxial.png", dpi = 300)

plt.show()
