import sys
sys.path += ["../"]

import smbh
import nolds
import numpy as np
import matplotlib.pyplot as plt

pos = 1e-0 * np.ones(3) / 3 ** 0.5

v = 50
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

# results = smbh.run(speeds, pot_type = smbh.SYMMETRIC)
# r = results.distance
#
# l = nolds.lyap_e(r)
# print(l)
#
# results = smbh.run(speeds, pot_type = smbh.TRIAXIAL)
# r = results.distance
# 
# l = nolds.lyap_e(r)
# print(l)
