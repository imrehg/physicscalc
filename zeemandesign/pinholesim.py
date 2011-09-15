"""
Calculating and simulating the angular distribution after a collimator
in the Zeeman slower
"""
import numpy as np
import pylab as pl
import sys

try:
    import multiprocessing as processing
except ImportError:
    import processing

def oneShot(args):
    """ 
    Do one simulation of atoms escaping
    """
    r1, r2, l, random = args
    while True:
        rval = random.uniform(0, 1, 4)
        x1, y1 = r1*rval[0], r1*rval[1]
        if (x1*x1 + y1*y1) <= r1*r1:

            th = np.arccos(rval[2])
            th = np.arcsin(rval[2])
            phi = 2*np.pi*rval[3]
            dl = np.sin(th)*l
            x2 = x1 + np.cos(phi)*dl
            y2 = y1 + np.sin(phi)*dl
            break    
    # return x1, y1, th, phi, x2, y2
    if (x2*x2 + y2*y2) <= r2*r2:
        res = True
    else:
        res = False
    return (th, phi, res)

def runsimu(repeats, r1=1, r2=1, l=1, pool=None, seed=None):
    
    random = np.random.mtrand.RandomState(seed)
    
    inargs = (r1, r2, l, random)

    TASKS = [(inargs)]*repeats
    if pool:
        try:
            result = pool.map_async(oneShot, TASKS).get(9999999)
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(0)
    else:
        result = map(oneShot, TASKS)
    return result

if __name__ == "__main__":
    import time

    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)
    pool = processing.Pool(processes=4)
    # pool = None

    # # Big series
    # for r1 in [0.5, 1, 2.0]:
    #     for d in [1, 5, 10, 20]:
    #         sim = {'r1': r1,
    #                'r2': 1,
    #                'd': d,
    #                'repeats': 4000000}
    #         filename = "pinhole_%g_%d" %(r1, d)
    #         print sim
    #         start = time.time()
    #         result = runsimu(sim['repeats'],
    #                          sim['r1'],
    #                          sim['r2'],
    #                          sim['d'],
    #                          pool=pool,
    #                          seed=None,
    #                          )
    #         elapsed = time.time()-start
    #         print "Elapsed time: %.3fs (%.1f repeats/s)" %(elapsed, sim['repeats']/elapsed)
    #         np.savez(filename, sim=sim, result=result)

    r1, d = 1, 1
    sim = {'r1': 1,
           'r2': 1,
           'd': 1,
           'repeats': 4000000}
    filename = "pinhole_%g_%d" %(r1, d)
    print sim
    start = time.time()
    result = runsimu(sim['repeats'],
                     sim['r1'],
                     sim['r2'],
                     sim['d'],
                     pool=pool,
                     seed=None,
                     )
    elapsed = time.time()-start
    print "Elapsed time: %.3fs (%.1f repeats/s)" %(elapsed, sim['repeats']/elapsed)
    np.savez(filename, sim=sim, result=result)
