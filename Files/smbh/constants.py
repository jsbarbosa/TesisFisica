from numpy import pi, log

def setDarkMatterScaleRadius(radius):
    global DARK_MATTER_SCALE_RADIUS, DARK_MATTER_DENSITY_0
    DARK_MATTER_SCALE_RADIUS = radius
    DARK_MATTER_DENSITY_0 = HALO_MASS / (4 * pi * DARK_MATTER_SCALE_RADIUS ** 3 * (log((DARK_MATTER_SCALE_RADIUS + R_VIR) / DARK_MATTER_SCALE_RADIUS) - R_VIR / (DARK_MATTER_SCALE_RADIUS + R_VIR)))

"""
Universal constants
"""
G0 = 6.67408e-11
kpc = 1 / (3.0857e19 ** 3)
gyr = (3.154e16) ** 2
m0 = (1.98847e35)
G = G0 * kpc * gyr * m0

H_0 = 67.66
kpc2 = 1 / (3.0857e16) # 1 kpc in km
gyr2 = 3.154e16
H = H_0 * kpc2 * gyr2 / 1000

"""
Simulation constants
"""
HALO_MASS = 1e3
LN_LAMBDA = 2.3
SOLAR_MASS = 1e-5

BARIONIC_DENSITY_0 = 100
TOTAL_MASS = 0.158 * HALO_MASS # fb * Mh

R_VIR = 200  * 3 * H ** 2 / (8 * pi * G)
SCALE_LENGTH = 0.01 * R_VIR / (1 + 2**0.5)

DARK_MATTER_SCALE_RADIUS = 0
DARK_MATTER_DENSITY_0 = 0

MATTER_DENSITY_PARAMETER = 0.309
h = 0.678

setDarkMatterScaleRadius(1)
