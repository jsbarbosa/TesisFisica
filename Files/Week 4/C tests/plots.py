import os
import numpy as np
import matplotlib.pyplot as plt

def makeArgs(name, r0, t0, p0, v0, vt0, vp0):
    return name + ' ' + ("%e " * 6)%(r0, t0, p0, v0, vt0, vp0)

def runSimulation(name, r0, t0, p0, v0, vt0, vp0):
    args = makeArgs(name, r0, t0, p0, v0, vt0, vp0)
    try: os.remove(name)
    except FileNotFoundError: pass

    os.system("./rebound.run %s" % args)

    data = np.genfromtxt(name)
    t = data[:, 0]
    pos = data[:, 1:4]
    speeds = data[:, 4:]

    return t, pos, speeds

def determineStop(time, pos, speeds):
    r = np.linalg.norm(pos, axis = 1)
    diff = abs(np.diff(r))

    # np.integrate()
    # p = np.where(diff < 1e-5)[0][0]
    # l = len(time)
    # if p * 1.25 < l:
    #     p = int(p * 1.25)
    # return time[:p], pos[:p], speeds[:p]

r0 = 1e-6
t0 = np.pi / 4
p0 = np.pi / 4

v0 = 70
vt0 = np.pi / 4
vp0 = np.pi / 4

filename = "Results.dat"

t, pos, speeds = runSimulation(filename, r0, t0, p0, v0, vt0, vp0)
# determineStop(t, pos, speeds)
# t, pos, speeds = determineStop(t, pos, speeds)

r = np.linalg.norm(pos, axis = 1)  / 0.254111
v = np.linalg.norm(speeds, axis = 1)

fig, ax = plt.subplots()
ax.plot(1000 * t, r)
ax.set_xlabel('Time since recoil (Myr)')
ax.set_ylabel(r'$R / R_{vir}$')

fig.savefig('recoil.png', dpi = 300)
plt.show()
