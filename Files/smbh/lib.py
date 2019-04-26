import os
import numpy as np
from scipy.special import erf
from functools import partial

from ctypes import cdll, CDLL, c_double, POINTER, c_int, c_char_p, pointer, ArgumentError

from .useful import Results

G = 0.449338
R_VIR_z20 = 0.6899734690284971
DEFAULT_STELLAR_FRACTION = 0.01
DEFAULT_GAS_POWER = -2.2

INT_LEAPFROG = 0
INT_IAS15 = 1
INT_WHFAST = 2
INT_SEI = 3
INT_JANUS = 4
INT_MERCURIUS = 5

SYMMETRIC = 0
TRIAXIAL = 1
C_SYMMETRIC = 2
C_TRIAXIAL = 3

dir = os.path.dirname(__file__)
os.environ['PATH'] = dir + ';' + os.environ['PATH']
path = os.path.abspath(os.path.join(dir, "rebound.so"))

try:
    lib = CDLL(path)
except OSError:
    cwd = os.getcwd()
    os.chdir(dir)
    lib = CDLL(path)
    os.chdir(cwd)

def setResArgs(func, res, args):
    f = getattr(lib, func)
    f.restype = res
    try:
        f.argtypes = tuple(args)
    except TypeError:
        if args != None: f.argtypes = args,
    setattr(lib, func, f)

def function(func, *args):
    if lib.getR_vir() == 0:
        lib.setR_vir(R_VIR_z20)
    try:
        return func(*args)
    except ArgumentError:
        values = np.array(args).T
        ls = [func(*row) for row in values]
        return np.array(ls)

func_info = [('getG', c_double, None),
            ('getR_vir', c_double, None),
            ('setR_vir', None, c_double),
            ('printConstants', None, None),
            ('gasDensity', c_double, c_double),
            ('gasMass', c_double, c_double),
            ('getNorm', c_double, POINTER(c_double)),
            ('darkMatterMass', c_double, c_double),
            ('darkMatterDensity', c_double, c_double),
            ('darkMatterDensity0', c_double, c_double),
            ('getLocalSoundSpeed', c_double, c_double),
            ('gravitationalForce', c_double, c_double),
            ('stellarMassHernquist', c_double, c_double),
            ('stellarDensityHernquist', c_double, c_double),
            ('SMBHAccretion', c_double, (c_double, c_double)),
            ('dynamicalFrictionDM', c_double, (c_double, c_double)),
            ('dynamicalFrictionGas', c_double, (c_double, c_double)),

            ('darkMatterPotential', c_double, c_double),
            ('stellarPotential', c_double, c_double),
            ('gasPotential', c_double, c_double),

            ('run', None, (2 * [POINTER(c_double)] + 3 * [c_double] + 3 * [c_int] + [c_char_p])),
            ('getRedshift', c_double, c_double),
            ('getHubbleParameter', c_double, c_double),
            ('calculateR_vir', c_double, (c_double, c_double)),
            ('dampingFactor', c_double, (c_double, c_double)),
            ('setBaryonicFraction', None, c_double),
            ('setStellarRatio', None, c_double),
            ('darkMatterVelocityDispersion', c_double, None),
            ('machFunction', c_double, c_double),
            ('setGasPower', None, c_double),
            ('darkMatterDensityTriaxial', c_double, (c_double, c_double, c_double)),
            ('getM', c_double, (c_double, c_double, c_double)),
            ('setTriaxialCoeffs', None, 3 * [c_double]),
            ('triaxial_gravDM', POINTER(c_double), 4 * [c_double]),
            ('triaxial_gravS', POINTER(c_double), 4 * [c_double]),
            ('triaxial_gravG', POINTER(c_double), 4 * [c_double]),
            ('testLoad', None, c_char_p),
            ('getReturnFraction', c_double, None),
            ('lyapunov', POINTER(c_double), 2 * [POINTER(c_double)] + 4 * [c_double] + [c_int, c_int]),

            ('darkMatterPotential_triaxial', c_double, 3 * [c_double]),
            ('stellarPotential_triaxial', c_double, 3 * [c_double]),
            ('gasPotential_triaxial', c_double, 3 * [c_double]),
            ('getTriaxialCoeffs', POINTER(c_double), None),

            ('getReturnProperties', None, c_char_p),
            ]

for func in func_info:
    name = func[0]
    setResArgs(*func)
    globals()[name] = partial(function, getattr(lib, name))

def pointerFunction(func, pos, speeds):
    pos = (c_double * 3)(*pos)
    speeds = (c_double * 3)(*speeds)
    return func(pos, speeds)

def pointerReturn(pointer):
    values = [pointer[i] for i in range(3)]
    lib.free(pointer)
    return np.array(values)

def run(speeds, pos = [1e-3 / (3**0.5)] * 3, smbh_mass = 1, dt = 1e-6, end_time = 0, pot_type = TRIAXIAL, \
        integrator = INT_LEAPFROG, save_every = 10, filename = None, read = True,
        header_only = False):
    global lib
    delete_file = False
    pos = (c_double * 3)(*pos)
    speeds = (c_double * len(speeds))(*speeds)
    if filename == None:
        filename = "temp.dat"
        delete_file = True
    lib.run(pos, speeds, smbh_mass, dt, end_time, pot_type, integrator, int(save_every), filename.encode())

    if read: data = Results(filename, header_only)
    else: data = None
    if delete_file: os.remove(filename)

    return data

def lyapunov(speeds, positions = [1e-3 / (3**0.5)] * 3, d_q0 = 1e-8, smbh_mass = 1, T = 1e-5, dt = 1e-6, l = 100, triaxial = True):
    pos = (c_double * 3)(*positions)
    speeds = (c_double * 3)(*speeds)

    pointer = lib.lyapunov(pos, speeds, d_q0, smbh_mass, T, dt, l, triaxial)
    values = [pointer[i] for i in range(6)]
    lib.free(pointer)
    return np.array(values)

def getTriaxialCoeffs():
    return pointerReturn(lib.getTriaxialCoeffs())

def getPotential(r):
    return darkMatterPotential(r) + stellarPotential(r) + gasPotential(r)

def getPotential_triaxial(x, y, z):
    return darkMatterPotential_triaxial(x, y, z) + \
            stellarPotential_triaxial(x, y, z) + \
            gasPotential_triaxial(x, y, z)

def getEscapeSpeed(r = R_VIR_z20, r0 = 1e-3):
    pot_1 = getPotential(r)
    # pot_0 = 0
    pot_0 = getPotential(r0)
    return abs(2 * (pot_1 - pot_0)) ** 0.5

def getEscapeSpeed_triaxial(r, r0 = [1e-3 / 3**0.5] * 3):
    pot_1 = getPotential_triaxial(*r)
    pot_0 = getPotential_triaxial(*r0)
    return abs(2 * (pot_1 - pot_0)) ** 0.5
# def triaxial_gravDM(x, y, z, gamma = 0.2):
#     return pointerReturn(lib.triaxial_gravDM(x, y, z, gamma))
#
# def triaxial_gravS(x, y, z, gamma = 0.2):
#     return pointerReturn(lib.triaxial_gravS(x, y, z, gamma))
#
# def triaxial_gravG(x, y, z, gamma = 0.2):
#     return pointerReturn(lib.triaxial_gravG(x, y, z, gamma))

def getEscapeSpeedOrthogonalTriaxial():
    a1, a2, a3 = getTriaxialCoeffs()

    x_vir = R_VIR_z20
    y_vir = (a2 / a1) * R_VIR_z20
    z_vir = (a3 / a1) * R_VIR_z20

    poss = [[x_vir, 0, 0],
            [0, y_vir, 0],
            [0, 0, z_vir]]

    e_v = [getEscapeSpeed_triaxial(pos) for pos in poss]
    return np.array(e_v)

G = getG()
setR_vir(R_VIR_z20)
setGasPower(2.2)
setStellarRatio(0.01)
