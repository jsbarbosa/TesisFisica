from time import time, sleep
from multiprocessing import Process

def startProcess(target, args):
    proc = Process(target = target, args = args, daemon = True)
    proc.start()
    return proc

def run_multiproccess(target, list_of_args, n_proccesors = 4):
    procs = []
    s_time = time()

    current = n_proccesors
    started = []
    for i, arg in enumerate(list_of_args[: current]):
        proc = startProcess(target, list_of_args[i])
        started.append(proc)

    while current < len(list_of_args):
        temp = []
        for i in range(n_proccesors):
            if started[i].is_alive(): temp.append(started[i])
            else: started[i].terminate()

        news = len(started) - len(temp)
        started = temp

        for i in range(news):
            if current < len(list_of_args):
                proc = startProcess(target, list_of_args[current])
                started.append(proc)
                current += 1
            else:
                break
        sleep(0.5)

    while sum([proc.is_alive() for proc in started]) != 0:
        sleep(0.5)

    e_time = time()
    print('DONE: %d s'% (e_time - s_time))
