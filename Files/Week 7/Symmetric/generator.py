import sys
sys.path += ["../.."]

import smbh
import numpy as np

import matplotlib.pyplot as plt

ratios = np.linspace(0.01, 0.1, 50)
speeds = np.linspace(95, 55, 9).astype(int)

file = "Data/return_%d_%d.dat"

for (i, speed) in enumerate(speeds):
    for (j, ratio) in enumerate(ratios):
        smbh.setStellarRatio(ratio)
        smbh.run([speed, 0, 0], save_every = 10, filename = file%(i, j), read = False, triaxial = False)
