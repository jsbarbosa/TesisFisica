import sys
sys.path += ["../"]

import smbh
import numpy as np

v0 = 30
theta = np.pi / 4
phi = theta

x = v0 * np.sin(theta) * np.cos(phi)
y = v0 * np.sin(theta) * np.sin(phi)
z = v0 * np.cos(theta)

smbh.run([x, y, z], save_every = 10, filename = "tOrbit.dat", read = False, triaxial = True)
