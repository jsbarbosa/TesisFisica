import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from smbh.lib import *
from smbh.useful import magnitude

mpl.rcParams['grid.color'] = 'k'
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['grid.linewidth'] = 0.5

def slicePlot(data, axes = None):
    if axes == None:
        gs = gridspec.GridSpec(2, 2)
        ax1 = plt.subplot(gs[1, 0])
        ax2 = plt.subplot(gs[0, 0])
        ax3 = plt.subplot(gs[1, 1])
    else:
        ax1, ax2, ax3 = axes

    c = ax1.plot(data[:, 0], data[:, 1])[0].get_color()
    ax1.plot(data[0, 0], data[0, 1], 'o', c = c, alpha = 0.5)
    ax1.plot(data[-1, 0], data[-1, 1], 'o', c = c, alpha = 0.7)
    ax1.plot([0], [0], 'o')
    ax2.plot(data[:, 0], data[:, 2])
    ax3.plot(data[:, 2], data[:, 1])

    ax2.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)

    ax1.set_xlabel('$x (kpc)$')
    ax1.set_ylabel('$y (kpc)$')
    ax2.set_ylabel('$z (kpc)$')
    ax3.set_xlabel('$z (kpc)$')

    return plt.gcf(), (ax1, ax2, ax3)

def plotDensityMass(distance, densities, masses, figsize = (8, 4.5)):
    fig, ax1 = plt.subplots(figsize = figsize)
    ax2 = ax1.twinx()

    ax1.set_xlabel('Distance from center (kpc)')
    ax1.set_ylabel(r'Density ($1\times10^5$ $M_\odot$kpc$^{-3}$)', color = 'b')
    ax1.tick_params('y', colors = 'b')

    ax2.set_ylabel(r'Mass ($1\times10^5$ $M_\odot$)', color = 'g')
    ax2.tick_params('y', colors = 'g')

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax2.set_yscale('log')

    n = len(densities)
    if n <= 3:
        lines = ["-", "--", ":"]
        for i in range(n):
            ax1.plot(distance, densities[i], lines[i], c = 'b')
            ax2.plot(distance, masses[i], lines[i], c = 'g')
    else:
        ax1.plot(distance, densities, c = 'b')
        ax2.plot(distance, masses, c = 'g')

    r_vir = getR_vir()

    x = [r_vir, r_vir]
    y = ax1.get_ylim()

    ax1.plot(x, y, lw = 0.5, c = 'r')

    ax1.set_ylim(y)

    fig.tight_layout()

    ax1.grid()

    return fig, (ax1, ax2)

def plotDensityMassForAll(r, r_vir = None, baryonic_fraction = None, stellar_ratio = None, figsize = (8, 4.5)):
    if r_vir != None: setR_vir(r_vir)
    if baryonic_fraction != None: setBaryonicFraction(baryonic_fraction)
    if stellar_ratio != None: setStellarRatio(stellar_ratio)

    d_dm = darkMatterDensity(r)
    m_dm = darkMatterMass(r)

    d_s = stellarDensityHernquist(r)
    m_s = stellarMassHernquist(r)

    d_g = gasDensity(r)
    m_g = gasMass(r)

    ds = [d_dm, d_s, d_g]
    ms = [m_dm, m_s, m_g]

    fig, (ax1, ax2) = plotDensityMass(r, ds, ms, figsize)

    ax1.legend(["Dark matter", "Stars", "Gas"])

    return fig, (ax1, ax2)

def plotDistance(times, positions, ax = None, figsize = (8, 4.5)):
    rs = np.linalg.norm(positions, axis = 1)

    if ax == None:
        fig, ax = plt.subplots(figsize = figsize)
        ax.set_xlabel("Time (Myr)")
        ax.set_ylabel("$R / R_{vir}$")
    else:
        fig = plt.gcf()

    y = rs / constants.R_VIR
    c = ax.plot(times * 1000, y)[0].get_color()

    return fig, ax

def plotProperties(results, figsize = (6, 8.5)):
    r_vir = results.R_VIR
    setR_vir(r_vir)
    pos = results.positions
    speeds = results.speeds
    accr = SMBHAccretion(pos, speeds)

    t = results.times * 1000
    r = results.distance
    v = results.speed
    masses = results.masses

    dm = dynamicalFrictionDM(r, v)
    dg = dynamicalFrictionGas(r, v)
    factor = dampingFactor(r, v)

    grav = gravitationalForce(r)
    rho_dm = darkMatterDensity(r)
    rho_h = stellarDensityHernquist(r)
    rho_g = gasDensity(r)
    m_dm = darkMatterMass(r)

    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize = figsize, sharex = True)

    ax11 = ax1.twinx()

    ax1.plot(t, r / r_vir, c = 'k')
    ax11.plot(t, v, c = 'b')

    ax2.plot(t, masses)

    ax3.plot(t, abs(grav), label = '$a_{grav}$')
    ax3.plot(t, abs(dm), label = '$a_{DF}$')
    ax3.plot(t, abs(dg), label = '$a_{gas}$')
    ax3.plot(t, rho_g, label = r'$\rho_{gas}$')
    ax3.plot(t, rho_h, label = r'$\rho_h$')
    ax3.plot(t, rho_dm, label = r'$\rho_{dm}$')

    ax1.grid()
    ax2.grid()
    ax3.grid()

    ax3.set_yscale('log')

    # ax3.legend(prop={'size': 8})
    ax3.legend(loc = 'center', bbox_to_anchor = (1.15, 0.5),
          ncol = 1, fancybox = True, shadow = True, prop = {'size': 12})

    ax1.set_ylabel('$R / R_{vir}$')
    ax11.set_ylabel('Speed (kpc/Gyr)', color = 'b')
    ax11.tick_params('y', colors='b')

    ax2.set_ylabel(r'$M_\bullet$ ($10^5M_\odot$)')
    ax3.set_xlabel('Time (Myr)')

    fig.tight_layout()
    fig.subplots_adjust(hspace = 0)

    return fig, (ax1, ax2, ax3)

def make3dPlot(positions):
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2])
    ax.set_xlabel("$x$ (kpc)")
    ax.set_ylabel("$y$ (kpc)")
    ax.set_zlabel("$z$ (kpc)")

    return fig, ax

def makeAnimation(positions, points = 100):
    def update_lines(i):
        x, y, z = data[:i].T
        if i > 0:
	        line.set_xdata(x)
	        line.set_ydata(y)
	        line.set_3d_properties(z)
	        dot.set_xdata(x[-1])
	        dot.set_ydata(y[-1])
	        dot.set_3d_properties(z[-1])

        return line, dot

    fig = plt.figure()
    ax = Axes3D(fig)

    line = ax.plot([0], [0], [0])[0]
    dot = ax.plot([0], [0], [0], "o")[0]

    steps = len(positions) // points

    data = positions[:: steps]

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    ax.set_xlim3d([data[:, 0].min(), data[:, 0].max()])
    ax.set_ylim3d([data[:, 1].min(), data[:, 1].max()])
    ax.set_zlim3d([data[:, 2].min(), data[:, 2].max()])

    return fig, ax, animation.FuncAnimation(fig, update_lines, frames = len(data), interval = 50, blit = False)
