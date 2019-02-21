import os
import numpy as np
from scipy.special import erf
from functools import partial

from ctypes import cdll, CDLL, c_double, POINTER, c_int, c_char_p, pointer, ArgumentError

from .useful import readFile

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

G = 0.449338

def setResArgs(func, res, args):
    f = getattr(lib, func)
    f.restype = res
    try:
        f.argtypes = tuple(args)
    except TypeError:
        if args != None:f.argtypes = args,
    setattr(lib, func, f)

def function(func, *args):
    try:
        return func(*args)
    except ArgumentError:
        if len(args) == 1:
            ls = [func(arg) for arg in args[0]]
            return np.array(ls)
        else:
            raise('Put argument one by one')

func_info = [('setR_vir', None, c_double),
            ('printConstants', None, None),
            ('gasDensity', c_double, c_double),
            ('getNorm', c_double, POINTER(c_double)),
            ('darkMatterMass', c_double, c_double),
            ('darkMatterDensity', c_double, c_double),
            ('darkMatterDensity0', c_double, c_double),
            ('getLocalSoundSpeed', c_double, c_double),
            ('gravitationalForce', c_double, c_double),
            ('stellarMassHernquist', c_double, c_double),
            ('stellarDensityHernquist', c_double, c_double),
            ('SMBHAccretion', c_double, (POINTER(c_double), POINTER(c_double))),
            ('dynamicalFrictionDM', POINTER(c_double), (POINTER(c_double), POINTER(c_double))),
            ('dynamicalFrictionGas', POINTER(c_double), (POINTER(c_double), POINTER(c_double))),
            ('run', None, (POINTER(c_double), POINTER(c_double), c_double, c_double, c_int, c_int, c_char_p)),
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

def run(positions, speeds, smbh_mass, dt, n_points, n_saves, filename, delete_file = True):
    global lib
    pos = (c_double * len(positions))(*positions)
    speeds = (c_double * len(speeds))(*speeds)
    try: os.remove(filename)
    except FileNotFoundError: pass
    try:
        lib.run(pos, speeds, smbh_mass, dt, int(n_points), int(n_saves), filename.encode())
    except Exception as e:
        print(e)

    data = readFile(filename)

    if delete_file: os.remove(filename)

    return data

def SMBHAccretion(pos, speeds):
    if type(pos) is np.ndarray:
        if len(pos.shape) > 1:
            return [pointerFunction(lib.SMBHAccretion, pos[i], speeds[i]) for i in range(len(pos))]
        else:
            return pointerFunction(lib.SMBHAccretion, pos, speeds)
    return pointerFunction(lib.SMBHAccretion, pos, speeds)

def dynamicalFrictionDM(pos, speeds):
    if type(pos) is np.ndarray:
        if len(pos.shape) > 1:
            ans = []
            for i in range(len(pos)):
                pointer = pointerFunction(lib.dynamicalFrictionDM, pos[i], speeds[i])
                values = pointerReturn(pointer)
                ans.append(values)
            return np.array(ans)
        else:
            return pointerReturn(pointerFunction(lib.dynamicalFrictionDM, pos, speeds))
    return pointerReturn(pointerFunction(lib.dynamicalFrictionDM, pos, speeds))

def dynamicalFrictionGas(pos, speeds):
    if type(pos) is np.ndarray:
        if len(pos.shape) > 1:
            return np.array([pointerReturn(pointerFunction(lib.dynamicalFrictionGas, pos[i], speeds[i])) for i in range(len(pos))])
        else:
            return pointerReturn(pointerFunction(lib.dynamicalFrictionGas, pos, speeds))
    return pointerReturn(pointerFunction(lib.dynamicalFrictionGas, pos, speeds))
