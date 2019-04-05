import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt

pos = 1e-3 * np.ones(3) / 3 ** 0.5

v = 25
theta = np.pi / 4
phi = np.pi / 4

x = v * np.sin(theta) * np.cos(phi)
y = v * np.sin(theta) * np.sin(phi)
z = v * np.cos(theta)

speeds = [x, y, z]

pos_ = pos.copy()

dq = 1e-4

coeffs = [[1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1, 0],
            [0, 0, 1],
            [0, 0, -1]]

coeffs = np.array(coeffs)

results_0 = smbh.run(speeds, pos, pot_type = 1)
results_0.setDistance(1, 0.5, 0.25)

rvir = results_0.R_VIR
fig, ax = plt.subplots(figsize = (8, 4.5))

ax.plot(results_0.times * 1000, results_0.distance / rvir)

labels = ["$x_+$", "$x_-$",
            "$y_+$", "$y_-$",
            "$z_+$", "$z_-$",]

for (i, d) in enumerate(coeffs):
    results_1 = smbh.run(speeds, pos + d * dq, pot_type = 1)
    results_1.setDistance(1, 0.5, 0.25)
    ax.plot(results_1.times * 1000, results_1.distance / rvir, label = labels[i])

ax.legend()
ax.grid()

ax.set_xlabel('Time (Myr)')
ax.set_ylabel('$R / R_{vir}$')

fig.savefig('lyapunovTrayectories.png', dpi = 300)
plt.show()
