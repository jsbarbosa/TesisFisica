import math
import numpy as np

from scipy.constants import G as G_0 # m3s-2kg-1

KPC = 3.0856776e19 # in meters
GYR = 60 * 60 * 24 * 365.25 * 1e9 # in seconds
SOLAR_MASSES = 1.98847e30 # in kg

def GToNaturalUnits():
    global KPC, GYR, SOLAR_MASSES, G_0
    distance = KPC ** -3
    time = GYR ** 2
    mass = SOLAR_MASSES * 1e5 # 1e5 solar masses base unit

    return G_0 * distance * time * mass

def HubbleToNaturalUnits(H_0):
    global KPC, GYR
    kpc2 = 1000 / KPC # 1 kpc in km
    return H_0 * kpc2 * GYR / 1000

def calculateR_vir(G, H):
    return ((G * 1e3) / (100 * H**2)) ** (1/3.)

def sphericalToCartesian(r, theta, phi):
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)

    return x, y, z

def magnitude(matrix):
    return np.linalg.norm(matrix, axis = 1)

def readFile(filename):
    data = np.genfromtxt(filename)
    t = data[:, 0]
    pos = data[:, 1:4]
    speeds = data[:, 4:7]
    mass = data[:, -1]

    return t, pos, speeds, mass
