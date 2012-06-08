from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ
import sys

if __name__ == "__main__":
    print "Analyze many"
    filename = "simloads.npz"
    try:
        data = np.load(filename)
    except (IOError):
        print "No such data"
        sys.exit(1)

    L = float(raw_input("Length of collimator"))

    res = data['simdata']
    maxth = np.arctan(2/L)
    n, nx = np.shape(res)
    onepin = reduce(lambda x, y: x+1 if y[2]<=maxth else x, res, 0)
    thlist = [y[2] for y in res]
    rsq = lambda x, y, th, phi: (x + L * np.tan(th) * np.cos(phi))**2 + (y + L * np.tan(th) * np.sin(phi))**2
    thlist2 = [y[2] for y in res if rsq(y[0], y[1], y[2], y[3]) <= 1]
    twopin = len(thlist2)
    onepin = reduce(lambda x, y: x+1 if y[2]<=maxth else x, res, 0)

    print "Maxtheta:",maxth
    print "Total atoms:",n 
    print "Counts:", twopin, "twopinhole;", onepin, "onepinhole"
    print "Capture fraction:",twopin/n, onepin/n
    # print twopin/onepin, 4/(3*np.pi), 8/(3*L)/(1+ 8/(3*L))
    print "Predicted capture fraction:", 8/(3*L*np.pi), np.sin(maxth)

    pl.figure()
    pdf, bins, patches = pl.hist(thlist, 50)
    pl.xlabel(r'\theta')
    pl.title('single hole')

    # pl.figure()
    # rings = 2*np.pi*(np.cos(bins[:-1])-np.cos(bins[1:]))
    # angles = bins[:-1]+np.diff(bins)/2
    # # pl.plot(angles, pdf*rings)
    # pl.plot(angles, np.cos(angles)*np.sin(angles))

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
    rings = 2*np.pi*(np.cos(bins[:-1])-np.cos(bins[1:]))

    pl.figure()
    # pl.polar(bins[1:], pdf/n*40)
    # pl.polar(bins, np.cos(bins), 'r--')
    pl.plot(bins[1:], pdf/n/rings, label='pdf')
    # pl.plot(bins[1:], np.cos(bins[1:])*n*np.diff(bins), 'r--', label='singlehole theo')
    pl.xlabel(r'\theta')
    pl.legend(loc='best')
    pl.title('double hole')

    pl.show()
