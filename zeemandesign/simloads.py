"""
Basic calculation about the Zeeman slower geometry
because math is hard.
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ

if __name__ == "__main__":
    seed = 0
    random = np.random.mtrand.RandomState()

    simdata = []
    # nrep = 1000000
    print "Simulate"
    temp = raw_input("Number of repetitions (default 100000): ")
    nrep = 100000 if temp == '' else int(temp)
    
    for i in xrange(nrep):
        x1, y1, a, b = random.rand(4)
        r = np.sqrt(x1*x1 + y1*y1)
        while r > 1:
            x1, y1 = random.rand(2)
            r = np.sqrt(x1*x1 + y1*y1)
        th = np.arcsin(np.sqrt(a))
        phi = b * 2 * np.pi
        simdata += [(x1, y1, th, phi)]

    np.savez("simloads.npz", simdata=simdata)
