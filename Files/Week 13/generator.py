import sys
sys.path += ["../"]

import smbh
import numpy as np
from time import sleep, time
import matplotlib.pyplot as plt

from multiprocessing import Process

ratios = np.linspace(0.01, 0.1, 50)
speeds = np.linspace(90, 55, 36).astype(int)

file = "Triaxial/return_%d_%d.dat"

def do(speed, i, j):
    global ratios, file
    smbh.setStellarRatio(ratios[j])
    smbh.setTriaxialCoeffs(1, 1, 1)
    print('Process %d, %d started' % (i, j))
    t0 = time()
    smbh.run([speed, 0, 0], save_every = 50, filename = file%(i, j), read = False, triaxial = True)
    print('Process %d, %d ended (%d s)' % (i, j, time() - t0))

def startProcess(args):
    proc = Process(target = do, args = args, daemon = True)
    proc.start()
    return proc

if __name__ == '__main__':
    procs = []
    s_time = time()
    args = []
    for (i, speed) in enumerate(speeds):
        if i <= 5:
            args += [(speed, i, j) for j in range(len(ratios))]

    p = 4
    p1 = p
    started = []
    for i, arg in enumerate(args[: p1]):
        proc = startProcess(arg)
        started.append(proc)

    while p1 < len(args):
        temp = []
        for i in range(p):
            if started[i].is_alive(): temp.append(started[i])
            else: started[i].terminate()

        news = len(started) - len(temp)
        started = temp

        for i in range(news):
            if p1 < len(args):
                proc = startProcess(args[p1])
                started.append(proc)
                p1 += 1
            else:
                break
        sleep(0.5)

    while sum([proc.is_alive() for proc in started]) != 0:
        sleep(0.5)

    e_time = time()
    print('DONE: %d s'% (e_time - s_time))
