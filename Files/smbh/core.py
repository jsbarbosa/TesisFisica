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
    return constants.TOTAL_MASS * constants.SCALE_LENGTH / (2 * np.pi * r * (r + constants.SCALE_LENGTH)**3)

def baryonicMassHernquist(r):
    return constants.TOTAL_MASS * (r / (r + constants.SCALE_LENGTH))**2

def dynamicalFrictionDM(position, speed, mass):    
    r = np.linalg.norm(position)
    v = np.linalg.norm(speed)
    
    factor = -4 * np.pi * (constants.G / v) ** 2 * mass * (darkMatterDensity(r) + BARYONIC_DENSITY_FUNCTION(r)) * constants.LN_LAMBDA
    
    sigma = (0.5 * constants.G * darkMatterMass(r) / constants.R_VIR) ** 0.5
    x = v / (2 ** 0.5 * sigma)
    
    factor *= erf(x) - 2 * np.pi ** -0.5 * x * np.exp(-x**2)
    
    return factor * speed / v


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
    
    particle = sim.particles[0]
    
    sim.additional_forces = additional_force
    sim.force_is_velocity_dependent = velocity_dependent
    
    return sim, particle

def runSimulation(sim, particle, n_points):
    times = np.zeros(n_points)
    positions = np.zeros((n_points, 3))
    speeds = np.zeros((n_points, 3))

    for i in range(n_points):
        times[i] = times[i - 1] + sim.dt
        sim.integrate(times[i], exact_finish_time = 0) # integrate to the nearest timestep so WHFast's timestep stays constant
        positions[i] = particle.xyz
        speeds[i] = particle.vxyz
    return times, positions, speeds


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

def plotDistance(times, positions, ax = None, figsize = (8, 4.5)):
    rs = np.linalg.norm(positions, axis = 1)
    
    if ax == None:
        fig, ax = plt.subplots(figsize = figsize)
        ax.set_xlabel("Time (Myr)")
        ax.set_ylabel("$R / R_{vir}$")
    else:
        fig = plt.gcf()
        
    ax.plot(times * 1000, rs / constants.R_VIR)
    
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
