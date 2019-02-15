import rebound
from scipy.special import erf

import numpy as np

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

import smbh.constants as constants

def darkMatterDensity(r):
    r += constants.SOFTENING_RADIUS
    factor = r / constants.DARK_MATTER_SCALE_RADIUS
    return constants.DARK_MATTER_DENSITY_0 / (factor * (1 + factor) ** 2)

def darkMatterMass(r):
    factor = np.log(1 + r / constants.DARK_MATTER_SCALE_RADIUS) - r / (constants.DARK_MATTER_SCALE_RADIUS + r)
    return 4 * np.pi * constants.DARK_MATTER_DENSITY_0 * factor * constants.DARK_MATTER_SCALE_RADIUS ** 3

def darkMatterConcentration(mh, z):
    c_0 = (4.58 / 2) * (((1 + z) / 2.24) ** 0.107 + ((1 + z) / 2.24) ** -1.29)
    alpha = -0.0965 * np.exp(- z / 4.06)
    factor = mh / (1e13 * constants.SOLAR_MASS)
    return c_0 * factor ** alpha

def baryonicDensityOld(r):
    radio = 1e-3
    def f1(r):
        return constants.BARIONIC_DENSITY_0
    def f2(r):
        return constants.BARIONIC_DENSITY_0 * (radio / r) ** (2.2)
    try:
        if r < radio:
            return f1(r)
        else:
            return f2(r)
    except ValueError:
        answer = np.zeros_like(r)
        answer[r < radio] = f1(r[r < radio])
        answer[r >= radio] = f2(r[r >= radio])
        return answer

def baryonicMassOld(r):
    def f1(r):
        return (4 / 3) * np.pi * constants.BARIONIC_DENSITY_0 * r ** 3
    def f2(r):
        return 5 * np.pi * constants.BARIONIC_DENSITY_0 * radio ** 2.2 * (r ** 0.8 - radio ** 0.8) + f1(radio)
    radio = 1e-3
    try:
        if r < radio:
            return f1(r)
        else:
            return f2(r)
    except ValueError:
        answer = np.zeros_like(r)
        answer[r < radio] = f1(r[r < radio])
        answer[r >= radio] = f2(r[r >= radio])
        return answer

def baryonicDensityHernquist(r):
    r += constants.SOFTENING_RADIUS
    return constants.TOTAL_MASS * constants.SCALE_LENGTH / (2 * np.pi * r * (r + constants.SCALE_LENGTH)**3)

def baryonicMassHernquist(r):
    return constants.TOTAL_MASS * (r / (r + constants.SCALE_LENGTH))**2

def sphericalToCartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z

def setupSimulation(mass, position, speed, additional_force, velocity_dependent = 1):

    sim = rebound.Simulation()
    sim.G = constants.G

    sim.add(m = mass, x = position[0], y = position[1], z = position[2],
            vx = speed[0], vy = speed[1], vz = speed[2])

    constants.SMBH_MASS = mass
    constants.particle = sim.particles[0]
    sim.additional_forces = additional_force
    sim.force_is_velocity_dependent = velocity_dependent

    return sim

def runSimulation(sim, n_points, save_points = 1000):
    if save_points > n_points: save_points = n_points

    every = n_points // save_points

    j = 0
    times = np.zeros(save_points)
    positions = np.zeros((save_points, 3))
    speeds = np.zeros((save_points, 3))

    for i in range(n_points):
        if (i % every == 0) & (j < save_points):
            times[j] = i * sim.dt
            positions[j] = constants.particle.xyz
            speeds[j] = constants.particle.vxyz
            j += 1
        constants.SIM_DT = sim.dt
        sim.integrate(i * sim.dt, exact_finish_time = 0) # integrate to the nearest timestep so WHFast's timestep stays constant

    return times, positions, speeds

def getLocalSoundSpeed(z):
    # TODO: change of units
    return 1.8 * (1 + z) ** 0.5 * (constants.HALO_MASS / (1e7 * constants.SOLAR_MASS)) ** (1/3) * ((constants.MATTER_DENSITY_PARAMETER * constants.h ** 2) / 0.14)

def gasDensity(r):
    return baryonicDensityHernquist(r) / 1000

def dynamicalFrictionDM(position, speed):
    r = np.linalg.norm(position)
    v = np.linalg.norm(speed) + constants.SOFTENING_RADIUS

    mass = constants.HALO_MASS
    # mass = darkMatterMass(r)
    # print()
    sigma = (0.5 * constants.G * mass / constants.R_VIR) ** 0.5
    x = v / ((2 ** 0.5) * sigma)

    rho = darkMatterDensity(r) + BARYONIC_DENSITY_FUNCTION(r)
    factor = -4 * np.pi * (constants.G ** 2) * constants.SMBH_MASS * rho * constants.LN_LAMBDA \
                * (erf(x) - 2 * (np.pi ** -0.5) * x * np.exp(-x**2))

    return factor * speed / (v ** 3)

def dynamicalFrictionGas(position, speed):
    cs = getLocalSoundSpeed(0)
    v = np.linalg.norm(speed) + constants.SOFTENING_SPEED
    match = v / cs
    if match <= 1.7:
        factor = erf(match / (2 ** 0.5)) - (2 / np.pi) ** 0.5 * match * np.exp(-0.5 * match**2)
        if match <= 0.8:
            f = 0.5 * constants.LN_LAMBDA * factor
        else:
            f = 1.5 * constants.LN_LAMBDA * factor
    else:
        f = 0.5 * np.log(1 - match ** -2) + constants.LN_LAMBDA

    rho = gasDensity(np.linalg.norm(position))
    return -4 * np.pi * constants.G ** 2 * constants.SMBH_MASS * rho * f * speed / (v ** 3)
    # return -4 * np.pi * (constants.G / (v + constants.SOFTENING_SPEED)) ** 2 * mass * rho * f * (speed / v)

def SMBHAccretion(position, speed):
    r = np.linalg.norm(position)
    v = np.linalg.norm(speed)
    return 4 * np.pi * (constants.G * constants.SMBH_MASS) ** 2 * baryonicDensityHernquist(r) / (getLocalSoundSpeed(0) ** 2 + v ** 2) ** (1.5)

def gravitationalForce(r):
    return -constants.G * (darkMatterMass(r) + baryonicMassHernquist(r)) / ((r + constants.SOFTENING_RADIUS) ** 2)

def baseCase(pos, speed):
    r = np.linalg.norm(pos)
    dir_ = pos / r

    grav = gravitationalForce(r)
    df_g = dynamicalFrictionGas(pos, speed)
    df_dm = dynamicalFrictionDM(pos, speed)
    df = df_g + df_dm
    # df = df_g

    m_change = SMBHAccretion(pos, speed)

    constants.SMBH_MASS += m_change * constants.SIM_DT

    v = np.linalg.norm(speed)
    accretion = v * m_change / constants.SMBH_MASS

    pos_dependent = grav + accretion

    return [pos_dependent * dir_[0] + df[0],
            pos_dependent * dir_[1] + df[1],
            pos_dependent * dir_[2] + df[2]]

"""
Sims
"""
def gravitational_only(sim):
    pos = np.array(constants.particle.xyz)
    r = np.linalg.norm(pos)
    coeff = gravitationalForce(r)
    dir_ = pos / r

    constants.particle.ax = coeff * dir_[0]
    constants.particle.ay = coeff * dir_[1]
    constants.particle.az = coeff * dir_[2]

def gravitational_DF(sim):
    pos = np.array(constants.particle.xyz)
    speed = np.array(constants.particle.vxyz)

    ax, ay, az = baseCase(pos, speed)

    constants.particle.ax = ax
    constants.particle.ay = ay
    constants.particle.az = az

"""
Plots
"""

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

def getZ(t):
    z1 = -1 + (1 -2*(1 - 1 / (constants.H * t)))**0.5
    return z1

def plotDensityMass(distance, density, mass, figsize = (8, 4.5)):
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

    ax1.plot(distance, density, c = 'b')
    ax2.plot(distance, mass, c = 'g')

    fig.tight_layout()

    return fig, (ax1, ax2)

# def plotDistance(times, positions, axes = None, figsize = (8, 4.5)):
#     rs = np.linalg.norm(positions, axis = 1)
#
#     if axes == None:
#         fig, ax1 = plt.subplots(figsize = figsize)
#         ax1.set_xlabel("Time (Myr)")
#         ax1.set_ylabel("$R / R_{vir}$")
#         ax2 = ax1.twiny()
#         ax2.invert_xaxis()
#     else:
#         fig = plt.gcf()
#         ax1, ax2 = axes
#
#     zs = getZ(times + 0.06335)
#     y = rs / constants.R_VIR
#     c = ax1.plot(times * 1000, y)[0].get_color()
#
#     ticks = ax1.get_xticks() / 1000
#     lims = ax1.get_xlim()
#     ticks[0] = lims[0] / 1000
#     ticks[-1] = lims[-1] / 1000
#
#     zs = getZ(0.04 + ticks)
#     zs_int = list(set(zs[1:-1].astype(int)))[::-1]
#
#     zs = [zs[0]] + zs_int + [zs[-1]]
#     labels = ["%d"%z for z in zs[1:-1]]
#     zs = zs[::-1]
#
#     ax2.set_xticks(zs[1:-1])
#     ax2.set_xticklabels(labels)
#     ax2.set_xlim(zs[0], zs[-1])
#     return fig, (ax1)

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

BARYONIC_DENSITY_FUNCTION = baryonicDensityHernquist
