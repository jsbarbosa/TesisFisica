import math
import numpy as np
from scipy.signal import argrelextrema
from scipy.constants import G as G_0 # m3s-2kg-1

KPC = 3.0856776e19 # in meters
GYR = 60 * 60 * 24 * 365.25 * 1e9 # in seconds
SOLAR_MASSES = 1.98847e30 # in kg

R_VIR_THRESHOLD = 0.01

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
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z

def magnitude(matrix):
    return np.linalg.norm(matrix, axis = 1)

def findLocalMaxima(array):
    return argrelextrema(array, np.greater)[0]

def findLocalMinima(array):
    return argrelextrema(array, np.less)[0]

class Results(object):
    def __init__(self, filename, header_only = False):
        with open(filename, "r") as file:
            for (i, line) in enumerate(file):
                line = line.replace("\n", "").replace(" ", "")
                try:
                    conf, val = line.split("=")
                    setattr(self, conf, float(val))
                except ValueError:
                    if "mass\t" in line:
                        i = line.find("mass\t")
                        line = line[i + 5:]
                    else:
                        line = file.readline()
                    vars = line.split('\t')
                    self.t0, self.x0, self.y0, self.z0, self.vx0, self.vy0, self.vz0, self.mass0 = [float(var) for var in vars]
                    break

        if not header_only:
            data = np.genfromtxt(filename, skip_header = i + 1)
            self.times = data[:, 0]
            self.positions = data[:, 1:4]
            self.speeds = data[:, 4:7]
            self.masses = data[:, 7]
            # self.lyapunov = data[:, -1]

            self.speed = magnitude(self.speeds)
            self.setDistance()

            self.return_index = -1

    def setDistance(self, a_1 = 1, a_2 = 1, a_3 = 1):
        x, y, z = self.positions.T
        self.distance = (x ** 2 + (a_1 * y / a_2) ** 2 + (a_1 * z / a_3) ** 2) ** 0.5

    def getMassAt(self, time):
        if time > self.RETURN_TIME:
            m, b = np.polyfit(self.times, np.log(self.masses), 1)
            mass = np.exp(m * time + b)
        else:
            p0, p1 = abs(self.times - time).argsort()[:2]
            m = (self.masses[p1] - self.masses[p0]) / (self.times[p1] - self.times[p0])
            mass = self.masses[p0] + (time - self.times[p0]) * m
        return mass

    def getReturnTime(self, threshold = R_VIR_THRESHOLD):
        try:
            return self.return_time
        except AttributeError:
            try:
                self.return_index = np.where(self.distance >= self.R_VIR * R_VIR_THRESHOLD)[0][-1]
            except IndexError: self.return_index = 0
            self.return_time = self.times[self.return_index]
            return self.return_time

    def getMassAtReturnTime(self, threshold = R_VIR_THRESHOLD):
        if self.return_index < 0:
            getReturnTime(threshold)
        return self.masses[self.return_index]

    def getInitialSpeed(self):
        s0 = [self.vx0, self.vy0, self.vz0]
        s = sum([v ** 2 for v in s0])
        return s ** 0.5
