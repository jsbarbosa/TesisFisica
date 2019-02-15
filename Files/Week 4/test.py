import sys
sys.path += [".."]

import rebound
import numpy as np
from scipy.special import erf

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

from smbh import *

p = 1e-6
v = 60

ps = sphericalToCartesian(p, np.pi / 4, np.pi / 4)
vs = sphericalToCartesian(v, np.pi / 4, np.pi / 4)

sim = setupSimulation(1, ps, vs, gravitational_DF)

sim.t = 0
sim.dt = 5e-7
sim.integrator = 'leapfrog'

n = int(0.1 // sim.dt)
print(n)
times_df, positions_df, speeds_df = runSimulation(sim, n)

fig, ax = plotDistance(times_df, positions_df)

# ax.set_ylim(0, 1)
plt.show()
