import numpy as np
from copy import copy
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib import cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

from smbh.lib import *
from smbh.useful import magnitude

mpl.rcParams['grid.color'] = 'k'
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['grid.linewidth'] = 0.5

def getColors(n_points):
    return cm.jet(np.linspace(0, 1, n_points))

def rgbToHex(color):
    if len(color.shape) == 1:
        c = tuple((255 * color[:3]).astype(int))
        return '#%02x%02x%02x' % c
    else:
        return [rgbToHex(c) for c in color]

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

def plotDensityMass(distance, densities, masses, lines = [":", "--", "-.", "-"], figsize = (8, 4.5)):
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
    if n > 1:
        if lines != None:
            for i in range(n):
                ax1.plot(distance, densities[i], lines[i], c = 'b')
                ax2.plot(distance, masses[i], lines[i], c = 'g')
        else:
            for i in range(n):
                ax1.plot(distance, densities[i], c = 'b')
                ax2.plot(distance, masses[i], c = 'g')
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

    d = d_dm + d_s + d_g
    ds = [d_dm, d_s, d_g, d]

    m = m_dm + m_s + m_g
    ms = [m_dm, m_s, m_g, m]

    fig, (ax1, ax2) = plotDensityMass(r, ds, ms, figsize = figsize)

    ax1.legend(["Dark matter", "Stars", "Gas", "Cumulative"])

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

    t = results.times * 1000
    r = results.distance
    v = results.speed
    masses = results.masses

    dm = dynamicalFrictionDM(r, v)
    dg = dynamicalFrictionGas(r, v)
    factor = dampingFactor(r, v)
    accr = SMBHAccretion(r, v)

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

def plotOrbits(results, figsize = (5, 4)):
    t_max = 0
    fig, ax = plt.subplots(figsize = figsize)
    for result in results:
        t = 1000 * result.times
        if t[-1] > t_max: t_max = t[-1]
        ratio = result.STELLAR_FRACTION
        fb = result.FB
        y = result.distance / result.R_VIR
        ax.fill_between(t, 0, y, alpha = 0.7, zorder = 1, linewidth = 1, edgecolor = 'k')

    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel(r'$R / R_{vir}$')

    ax.set_xlim(0, t_max)
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.grid()
    return fig, ax

def coolLegend(ax, labels = None, loc = 'upper center', bbox_to_anchor = (0.5, 1.15), ncol = 4, **kwargs):
    if labels == None:
        ax.legend(loc = loc, bbox_to_anchor = bbox_to_anchor,
                  ncol = ncol, fancybox = True, shadow = True, **kwargs)
    else:
        ax.legend(labels, loc = loc, bbox_to_anchor = bbox_to_anchor,
                  ncol = ncol, fancybox = True, shadow = True, **kwargs)

def changeFontSize(ax, fs = 25, xbins = None, ybins = None):
    ax.xaxis.label.set_fontsize(fs)
    ax.yaxis.label.set_fontsize(fs)

    [tick.label.set_fontsize(fs) for tick in ax.yaxis.get_major_ticks()]
    [tick.label.set_fontsize(fs) for tick in ax.xaxis.get_major_ticks()]

    if xbins != None:
        ax.xaxis.set_major_locator(plt.MaxNLocator(xbins))

    if ybins != None:
        ax.yaxis.set_major_locator(plt.MaxNLocator(ybins))

def evaluateMesh(function, array, symmetric):
    n = len(array)
    if symmetric: f = lambda x, y, z: function((x * x + y * y + z * z)**0.5)
    else: f = function

    all = []
    for axis in range(3):
        answer = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if axis == 0: answer[j, i] = f(array[i], array[j], 0)
                elif axis == 1: answer[i, j] = f(0, array[j], array[i])
                elif axis == 2: answer[i, j] = f(array[j], 0, array[i])
        all.append(answer)
    return all

def slice3d(from_, to_, function, symmetric = True, points = 35, alpha = 1, label = r'$1\times10^5M_\odot$ kpc$^{-3}$'):
    fig = plt.figure()

    ax = fig.add_subplot(111, projection = '3d')

    distance = np.linspace(from_, to_, points)
    X, Y = np.meshgrid(distance, distance)
    PLANE =  0 * np.ones(X.shape)

    data_XY, data_XZ, data_YZ = evaluateMesh(function, distance, symmetric)

    norm = colors.LogNorm(vmin = data_XY.min(), vmax = data_XY.max())

    color_XY = cm.jet(norm(data_XY))
    color_XZ = cm.jet(norm(data_XZ))
    color_YZ = cm.jet(norm(data_YZ))

    ax.plot_surface(X, Y, PLANE, rstride = 1, cstride = 1, facecolors = color_XY,
                    alpha = alpha, shade = False)

    X, Z = np.meshgrid(distance, distance)
    ax.plot_surface(PLANE, X, Z, rstride = 1, cstride = 1, facecolors = color_XZ,
                    alpha = alpha, shade = False)

    Y, Z = np.meshgrid(distance, distance)
    surf = ax.plot_surface(Y, PLANE, Z, rstride = 1, cstride = 1, facecolors = color_YZ,
                    alpha = alpha, shade = False)

    ax.invert_yaxis()

    ax.set_xlabel("$x$ (kpc)")
    ax.set_ylabel("$y$ (kpc)")
    ax.set_zlabel("$z$ (kpc)")

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.grid(False)

    m = cm.ScalarMappable(cmap = plt.cm.jet, norm = norm)
    m.set_array([])
    cbar = fig.colorbar(m, shrink = 0.5, aspect = 5)

    cbar.ax.set_ylabel(label, rotation = 90)
    ax.set_aspect('equal')

    fig.tight_layout()

    return fig, ax, cbar

def generateCMap(base_r, base_g, base_b, alpha_start = 0.1, alpha_stop = 1):
    cdict = {'red':   ((0.0, base_r, base_r),
                   (1.0, base_r, base_r)),
         'green': ((0.0, base_g, base_g),
                   (1.0, base_g, base_g)),
         'blue':  ((0.0, base_b, base_b),
                   (1.0, base_b, base_b)),
        'alpha':((0.0, alpha_start, alpha_start),
                (1.0, alpha_stop, alpha_stop)),
            }
    return LinearSegmentedColormap('map', cdict)

def plotPhaseSpace(result, fig_axes = None, colors = 'Greys', linewidth = 1, figsize = (8, 4.5)):
    def make_segments(x, y):
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis = 1)
        return segments

    pos = result.positions.T
    speeds = result.speeds.T
    momentums = [v * result.masses for v in speeds]

    ymin = 0
    ymax = 0

    n = len(pos[0])
    x_labels = ["x", "y", "z"]

    array = np.linspace(0, 1, n)

    if fig_axes == None:
        fig, axes = plt.subplots(1, 3, sharey = True, figsize = figsize)
        fig_axes = True
    else:
        fig, axes = fig_axes
        fig_axes = False

    for i in range(3):
        x, y = pos[i], momentums[i]
        segments = make_segments(x, y)
        lc = LineCollection(segments, cmap = colors, array = array,
            linewidth = linewidth)

        axes[i].add_collection(lc)
        if fig_axes:
            axes[i].set_xlim(1.1 * x.min(), 1.1 * x.max())
            axes[i].set_xlabel("$%s$ (kpc)" % x_labels[i])

            axes[i].grid()

            if min(y) < ymin: ymin = min(y)
            if max(y) > ymax: ymax = max(y)

    if fig_axes:
        axes[0].set_ylim(1.1 * ymin, 1.1 * ymax)
        axes[0].set_ylabel(r'$\vec{P}$ (10$^5M_\odot$kpc/Gyr)')

    return fig, axes

def plotPhaseSpaces(results, max_points = 200, linewidth = 1, figsize = (8, 4.5)):
    m = len(results)
    mins = np.zeros((m, 6))
    maxs = np.zeros_like(mins)

    colors = plt.cm.jet(np.linspace(0, 1, m))[:, :3]

    fig, axes = plt.subplots(1, 3, sharey = True, linewidth = 1, figsize = figsize)

    for (i, result) in enumerate(results):
        color = generateCMap(*colors[i])
        skip = int(round(result.positions.shape[0] / max_points, 0))
        if skip > 0:
            result = copy(result)
            result.positions = result.positions[:: skip]
            result.speeds = result.speeds[:: skip]
            result.masses = result.masses[:: skip]
        mins[i, : 3] = result.positions.min(axis = 0)
        mins[i, 3 :] = result.speeds.min(axis = 0)
        maxs[i, : 3] = result.positions.max(axis = 0)
        maxs[i, 3 :] = result.speeds.max(axis = 0)
        fig, axes = plotPhaseSpace(result, colors = color, fig_axes = [fig, axes], linewidth = linewidth, figsize = (8, 4.5))

    x_labels = ['x', 'y', 'z']

    mins = mins.min(axis = 0)
    maxs = maxs.max(axis = 0)
    for (i, ax) in enumerate(axes):
        ax.grid()
        ax.set_xlabel('$%s$ (kpc)' % x_labels[i])
        ax.set_xlim(1.1 * mins[i], 1.1 * maxs[i])

    axes[0].set_ylabel(r'$\vec{P}$ (10$^5M_\odot$kpc/Gyr)')
    axes[0].set_ylim(1.1 * mins[3:].min(), 1.1 * maxs[3:].max())
    return fig, axes

def make3dPlot(positions, fig_axes = None, color = None, alpha = 1, lw = 1):
    if fig_axes == None:
        fig = plt.figure()
        ax = Axes3D(fig)
    else:
        fig, ax = fig_axes

    points = positions.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis = 1)

    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], color = color, alpha = alpha, lw = lw)
    ax.set_xlabel("$x$ (kpc)")
    ax.set_ylabel("$y$ (kpc)")
    ax.set_zlabel("$z$ (kpc)")

    return fig, ax

def makePolarPlot(r, theta, phi, colors = None, cmap = 'jet', alpha = 0.45):
    fig = plt.figure()
    ax1 = fig.add_subplot(121, polar = True)
    ax2 = fig.add_subplot(122, polar = True)
    ax1.scatter(theta, r, c = colors, cmap = 'jet', alpha = alpha)
    ax2.scatter(phi, r, c = colors, cmap = 'jet', alpha = alpha)

    ax1.set_theta_zero_location('N')

    ax1.set_thetamin(0)
    ax1.set_thetamax(90)

    ax2.set_thetamin(0)
    ax2.set_thetamax(90)

    # ax1.set_rorigin(-0.6)
    # ax2.set_rorigin(-0.6)

    theta_labels = [(r'$\theta = %d^\circ$' + '\t') % i for i in np.arange(0, 90 + 1, 20)]
    phi_labels = ['\t' + r'$\phi = %d^\circ$' % i for i in np.arange(0, 90 + 1, 20)]
    # r_labels = ['%.1f' % i for i in np.arange(0.4, 0.9 + 1, 0.1)]

    ax1.set_xticklabels(theta_labels)
    ax2.set_xticklabels(phi_labels)
    # ax2.set_yticklabels(r_labels, rotation = 90)
    ax2.set_rlabel_position(22.5)  # get radial labels away from plotted line

    ax1.set_yticklabels([])

    fig.subplots_adjust(hspace = 0, wspace = 0.025)

    return fig, (ax1, ax2)

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
