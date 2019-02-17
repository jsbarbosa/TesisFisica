import os
import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt

R_VIR = 1
def makeArgs(name, r0, t0, p0, v0, vt0, vp0):
    return name + ' ' + ("%e " * 7)%(r0, t0, p0, v0, vt0, vp0, R_VIR)

def runSimulation(name, r0, t0, p0, v0, vt0, vp0):
    args = makeArgs(name, r0, t0, p0, v0, vt0, vp0)
    try: os.remove(name)
    except FileNotFoundError: pass
    os.system("./rebound.run %s" % args)

    data = np.genfromtxt(name)
    t = data[:, 0]
    pos = data[:, 1:4]
    speeds = data[:, 4:7]
    mass = data[:, -1]

    return t, pos, speeds, mass

r0 = 1e-4
t0 = np.pi / 4
p0 = np.pi / 4

v0 = 70 * 1.023 # to kpc / gyr

vt0 = np.pi / 4
vp0 = np.pi / 4

n = 10
label_every = n // 10
r_virs = np.linspace(5, 7, n)
colors = cm.jet(np.linspace(0, 1, n))

rs = []

fig, ax = plt.subplots()

filename = "Results.dat"
for (i, R_VIR) in enumerate(r_virs):
    t, pos, speeds, mass = runSimulation(filename, r0, t0, p0, v0, vt0, vp0)

    r = np.linalg.norm(pos, axis = 1)  / R_VIR
    rs.append(r)
    if i%label_every == 0:
        ax.plot(1000 * t, r, c = colors[i], label = "$R_{vir} = %.2f$ kpc"%R_VIR)
    else:
        ax.plot(1000 * t, r, c = colors[i])

ax.set_xlabel('Time since recoil (Myr)')
ax.set_ylabel(r'$R / R_{vir}$')
# ax.set_ylim(0, 2)
ax.legend()
# fig.savefig('recoil%f.png'%R_VIR, dpi = 300)
plt.show()
