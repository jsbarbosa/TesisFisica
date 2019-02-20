import math
import numpy as np

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
