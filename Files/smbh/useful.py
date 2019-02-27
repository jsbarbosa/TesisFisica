import math
import numpy as np
from scipy.signal import argrelextrema
from scipy.constants import G as G_0 # m3s-2kg-1

KPC = 3.0856776e19 # in meters
GYR = 60 * 60 * 24 * 365.25 * 1e9 # in seconds
SOLAR_MASSES = 1.98847e30 # in kg

R_VIR_THRESHOLD = 0.1

def kmsTokpcGyr(kms):
    global GYR, KPC
    return kms * 1000 * GYR / KPC

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

def findLocalMaxima(array):
    return argrelextrema(array, np.greater)[0]

def findLocalMinima(array):
    return argrelextrema(array, np.less)[0]

class Results(object):
    def __init__(self, filename):
        with open(filename, "r") as file:
            for (i, line) in enumerate(file):
                line = line.replace("\n", "").replace(" ", "")
                try:
                    conf, val = line.split("=")
                    setattr(self, conf, float(val))
                except ValueError:
                    break

        data = np.genfromtxt(filename, skip_header = i + 1)
        self.times = data[:, 0]
        self.positions = data[:, 1:4]
        self.speeds = data[:, 4:7]
        self.masses = data[:, -2]
        self.lyapunov = data[:, -1]

        self.distance = magnitude(self.positions)
        self.speed = magnitude(self.speeds)

        self.return_index = -1

    def getReturnTime(self, threshold = R_VIR_THRESHOLD):
        try:
            return self.return_time
        except AttributeError:
            r = self.distance / self.R_VIR
            indexes = findLocalMaxima(r)
            try:
                first_in = np.where(r[indexes] < threshold)[0][0]
                first_in = indexes[first_in]
                for i in range(first_in, 0, -1):
                    if r[i] >= threshold: break
            except Exception as e:
                i = -1
            self.return_time = self.times[i]
            self.return_index = i
            return self.return_time

    def getMassAtReturnTime(self, threshold = R_VIR_THRESHOLD):
        if self.return_index < 0:
            getReturnTime(threshold)
        return self.masses[self.return_index]
