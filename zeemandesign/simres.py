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
    thlist2 = [y[2] for y in res if y[4]==1]
    print "Maxtheta:",maxth
    print "Total atoms:",n 
    print "Counts:", twopin, "twopinhole;", onepin, "onepinhole"
    print "Capture fraction:",twopin/n, onepin/n
    # print twopin/onepin, 4/(3*np.pi), 8/(3*L)/(1+ 8/(3*L))
    print "Predicted capture fraction:", 8/(3*L*np.pi), np.sin(maxth)

    pl.figure()
    pdf, bins, patches = pl.hist(thlist, 40)
    pl.xlabel(r'\theta')
    pl.title('single hole')

    pl.figure()
    # pl.polar(bins[1:], pdf/n*40)
    # pl.polar(bins, np.cos(bins), 'r--')
    pl.plot(bins[1:], pdf, label='pdf')
    pl.plot(bins[1:], np.cos(bins[1:])*n*np.diff(bins), 'r--', label='singlehole theo')
    pl.xlabel(r'\theta')
    pl.legend(loc='best')
    pl.title('single hole')

    thlist = thlist2
    pl.figure()
    pdf, bins, patches = pl.hist(thlist, 40)
    pl.xlabel(r'\theta')
    pl.title('double hole')

    pl.figure()
    # pl.polar(bins[1:], pdf/n*40)
    # pl.polar(bins, np.cos(bins), 'r--')
    pl.plot(bins[1:], pdf, label='pdf')
    pl.plot(bins[1:], np.cos(bins[1:])*n*np.diff(bins), 'r--', label='singlehole theo')
    pl.xlabel(r'\theta')
    pl.legend(loc='best')
    pl.title('double hole')

    pl.show()
