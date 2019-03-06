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

func_info = [('getR_vir', c_double, None),
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
            ('run', None, (POINTER(c_double), POINTER(c_double), c_double, c_double, c_int, c_int, c_int, c_char_p)),
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
            ('setTriaxalCoeffs', None, (c_double, c_double, c_double)),
            ('triaxial_gravDM', POINTER(c_double), (c_double, c_double, c_double)),
            ('triaxial_gravS', POINTER(c_double), (c_double, c_double, c_double)),
            ('triaxial_gravG', POINTER(c_double), (c_double, c_double, c_double)),
            ('triaxial_gravitationalDarkMatter', POINTER(c_double), (c_double, c_double, c_double, c_double)),
            ('triaxial_gravitationalStellar', c_double, (c_double, c_double, c_double, c_double, c_int))]

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

def run(speeds, smbh_mass = 1, dt = 1e-6, triaxial = True, integrator = INT_LEAPFROG, save_every = 10, filename = None, read = True):
    global lib
    delete_file = False
    pos = (1e-3 / (3**0.5)) * np.ones(3)
    pos = (c_double * 3)(*pos)
    speeds = (c_double * len(speeds))(*speeds)
    if filename == None:
        filename = "temp.dat"
        delete_file = True
    lib.run(pos, speeds, smbh_mass, dt, triaxial, integrator, int(save_every), filename.encode())

    if read: data = Results(filename)
    else: data = None
    if delete_file: os.remove(filename)

    return data

def SMBHAccretion(pos, speeds):
    if type(pos) is np.ndarray:
        if len(pos.shape) > 1:
            return [pointerFunction(lib.SMBHAccretion, pos[i], speeds[i]) for i in range(len(pos))]
        else:
            return pointerFunction(lib.SMBHAccretion, pos, speeds)
    return pointerFunction(lib.SMBHAccretion, pos, speeds)

setR_vir(R_VIR_z20)
setGasPower(-2.2)
setStellarRatio(0.01)
