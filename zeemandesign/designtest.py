"""
Repeat calculation from CLHung's thesis with Cesium
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ

import flux
import zeemanslower as zs
import cesiumrepeat as cs

## Constants
kB = 1.3806488e-23
##


if __name__ == "__main__":
    #### All parameters
    atom = zs.Rb85()
    T = 80 + 273 # Kelvin
    eta = 0.5
    l1 = 0
    l3 = 200e-3
    rmot = (25.4e-3)/2
    vf = 30
    #### No parameters after this

    l2 = 0.6

    v0new = np.sqrt(2 * l2 * atom.aslow * eta + vf**2)
    print("Maximum capture velocity: %g" %(v0new))
    u = np.sqrt(2 * kB * T / atom.m)
    print "Most probable velocity: ", u
    v0 = v0new /u
    vf /= u

    collimator = (1, 2*(l1+l2+l3)/rmot)
    Dvrmax = lambda x: rmot/cs.ttime(x, vf, l1, l2, l3, atom.aslow*eta, u)/u # Divergence limit
    Gvrmax = lambda x: rmot/(l1+l2+l3)*x # Geometric limit

    print "\nNo slowing", "="*10
    spinhole = (2/collimator[1])**2/(1+(2/collimator[1])**2)
    dpinhole = integ.dblquad(cs.intfunc, 0, np.inf, lambda x: 0, Gvrmax, args=collimator)[0]
    print("Single / Double pinhole / Ratio: %g / %g / %g" %(spinhole, dpinhole, spinhole/dpinhole))

    print "\nSlowing", "="*10
    GLim = integ.dblquad(cs.intfunc, 0, v0, lambda x: 0, Gvrmax, args=collimator)[0]
    DLim = integ.dblquad(cs.intfunc, 0, v0, lambda x: 0, Dvrmax, args=collimator)[0]
    print("No divergence / Divergence: %g / %g " %(GLim, DLim))
