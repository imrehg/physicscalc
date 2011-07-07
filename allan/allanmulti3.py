"""
Quick Allan variance calculation for the "allan_cont_long" type of data
(big stream of frequency records with constant time spacing)
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy as sp
import re
import os
import scipy.odr as odr
import ourgui
import time as t
import sys

try:
    import multiprocessing as processing
except:
    import processing

def apiter(base, limit):
    """
    Iterator over different orders of magnitudes
    base: list of base values (e.g. [1, 2, 5])
    limit: maximum value

    So apiter(base=[1, 2, 4],limit=300) will return an
    iterator through these values:
    1, 2, 5, 10, 20, 50, 100, 200
    """
    x, mul, i = base, 0, 0
    n = x[0] * 10**mul
    while n < limit:
        yield n
        i += 1
        if i >= len(x):
            i = 0
            mul += 1
        n = x[i] * 10**mul

def partallan(args):
    """ The core of the parallel Allan variation calculation code: a single set """
    y, num, n = args
    print "Region length: ", n
    divisor = 2 * n**2 * (num - 2 * n + 1)
    pp = [ (y[j+n] - y[j])  for j in xrange(0, num - n) ]
    part2 = np.array([np.sum(pp[i:(i+n)]) for i in xrange(0, num - 2*n+1)])
    result = np.sum(part2**2) / divisor
    return result


def allanover(freq, base):
    """ Allan variation with overlapping windows, parallel edition """
    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)  

    # # overlapping series
    y = (freq - base) / base
    num, s = len(y), []
    # iter = apiter([1, 2, 5], num/2)
    # iter = apiter(range(1, 10), num/2)
    iter = apiter([1], num/2)
    TASKS = [(y, num, n) for n in iter]
    try:
        out = pool.map_async(partallan, TASKS).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)
    s = [[TASKS[i][2], out[i]] for i in xrange(len(TASKS))]
    return s

def main(files, gatetime=0.1, marker='k.',title=None, highprec=True):
    """ Allan variance calculation and plotting """
    pl.figure(1)
    maxtime = 0
    maxgate = 0

    for filename in files:
        print filename
        data = np.loadtxt(filename, comments="#")
        filename = filename.split("/")[-1]
        base = np.mean(data)

        # Filter data
        if highprec:
            data = np.array([x for x in data if abs(x - base) < 5])

        n = len(data)
        time = np.array(range(0, n))*gatetime
        out = []

        variation = np.array(allanover(data, np.mean(data)))
        ogatetime = gatetime

        pl.subplot(211)
        pl.plot(time, data-np.mean(data), label=filename)
        pl.xlabel("Time (s)")
        pl.ylabel("Frequency deviation (Hz)")
        if time[-1] > maxtime:
            maxtime = time[-1]
        pl.xlim([0, maxtime])
        pl.legend(loc="best")

        pl.subplot(212)
        avgtime = variation[:,0]*ogatetime
        allandev = np.sqrt(variation[:,1])
        pl.loglog(avgtime, allandev, '.-', label=filename)
        pl.xlabel("Averaging time (s)")
        pl.ylabel("Allan deviation")
        if avgtime[-1] > maxgate:
            maxgate = avgtime[-1]
            pl.xlim([gatetime, maxgate])
        pl.legend(loc="best")
    pl.show()

if __name__ == "__main__":
    # Select multiple files
    files = []
    filename  = ourgui.openFile(type="log")
    while filename:
        files += [filename]
        filename  = ourgui.openFile(type="log")
    main(files, highprec=True)
