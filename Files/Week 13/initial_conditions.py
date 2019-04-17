import sys
sys.path += ["../"]

import smbh
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n_random = 500

v_min = 30
v_max = 85

thetas = 0.5 * np.pi * np.random.random(n_random)
phis = 0.5 * np.pi * np.random.random(n_random)
vs = (v_max - v_min) * np.random.random(n_random) + v_min

x, y, z = smbh.sphericalToCartesian(vs, thetas, phis)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.scatter(x, y, z, alpha = 0.5, s = 5)

ax.set_xlim(0, v_max)
ax.set_ylim(0, v_max)
ax.set_zlim(0, v_max)

# ax.set_xlim(-v_max, v_max)
# ax.set_ylim(-v_max, v_max)
# ax.set_zlim(-v_max, v_max)
plt.show()
