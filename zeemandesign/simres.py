from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ
import sys

if __name__ == "__main__":
    print "Analyze"
    Lin = int(raw_input("Length (if r=1): "))
    filename = "simdata%d.npz" %(Lin)

    try:
        data = np.load(filename)
    except (IOError):
        print "No such data"
        sys.exit(1)

    L = data['L']
    res = data['simdata']
    maxth = np.arctan(2/L)
    n, nx = np.shape(res)
    twopin = reduce(lambda x, y: x+1 if y[4]==1 else x, res, 0)
    onepin = reduce(lambda x, y: x+1 if y[2]<=maxth else x, res, 0)
    thlist = [y[2] for y in res]
    thlist2 = [y[2] for y in res if y[4] == 1 ]
    print "Maxtheta:",maxth
    print "Total atoms:",n 
    print "Counts:", twopin, "twopinhole;", onepin, "onepinhole"
    print "Capture fraction:",twopin/n, onepin/n

    pone = np.sin(maxth)**2
    th = np.linspace(0, maxth, 1001)
    y = np.tan(th) * L / 2
    a1 = 2*2/np.pi*(np.arccos(y) - y * np.sqrt(1 - y*y))*np.sin(th)*np.cos(th)
    ptwo = integ.trapz(a1, th)

    print "Predicted capture fraction:", ptwo, pone
    print "Reduction:", pone/ptwo
