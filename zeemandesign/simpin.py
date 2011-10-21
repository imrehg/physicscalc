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
    
    L = int(raw_input("Length (if r=1): "))

    for i in xrange(nrep):
        xp, yp, a, b = random.rand(4)

        x1 = 2*xp-1
        y1 = 2*yp-1
        r = np.sqrt(x1*x1 + y1*y1)
        while r > 1:
            x1, y1 = random.rand(2)
            r = np.sqrt(x1*x1 + y1*y1)
        th = np.arcsin(np.sqrt(a))
        phi = b * 2 * np.pi
        x2 = x1 + L * np.tan(th) * np.cos(phi)
        y2 = y1 + L * np.tan(th) * np.sin(phi)
        if np.sqrt(x2**2 + y2**2) <= 1:
            result = True
        else:
            result = False
        simdata += [(x1, y1, th, phi, result)]

    np.savez("simdata%d.npz" %(L), L=L, simdata=simdata)
