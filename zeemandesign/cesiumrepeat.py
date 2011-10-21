"""
Repeat calculation from CLHung's thesis with Cesium
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ

import flux
import zeemanslower as zs

## Constants
kB = 1.3806488e-23
##

def twopin(th, r, l):
    d = np.tan(th)*l
    c1 = r**2*np.arccos(d/(2*r))
    c3 = -0.5 * np.sqrt( d**2 * (4*r**2 - d**2 ))
    out = c1 + c1 + c3
    return out

def transtime(v0, vf, l1, l2, l3, a):
    t1 = l1 / v0
    t3 = l3 / vf
    L = (v0**2 - vf**2) / (2 * a)
    print L
    ta = (l2 - L) / v0
    tb = (v0 - np.sqrt(v0**2  - 2 * a * L)) / (2 * a)
    print t1, ta, tb, t3
    return t1+ta+tb+t3

def CsPv(T):
    return 10**(8.22127 - 4006.048 / T - 0.00060194 * T - 0.19623 * np.log10(T))

def doublePinFrac(r, L):
    return (L**2 + 2*r**2 - L*np.sqrt(L**2 + 4*r**2))/(2*r**2)

if __name__ == "__main__":

    #### All parameters
    atom = zs.Cs133()
    T = 60 + 273 # Kelvin
    P = 3e-5 # vapour pressure in torr
    r = 1e-3 # pinhole radius
    lcoll = 76.2e-3 # collimator length
    # eta = 0.5
    # v0 = 154
    # vf = 42
    eta = 0.82
    v0 = 205
    vf = 50
    r2 = 3.5e-3
    l1 = 200e-3
    l2 = 400e-3
    l3 = 150e-3
    #### No parameters after this
    P = CsPv(T) # fix vapour pressure from doc
    print "Calculated vapour pressure %g torr" %(P)

    l2new = zs.slowerlength(atom.aslow, eta, v0, vf)
    print "Slower correct length: %g" %(l2new)
    l2 = l2new


    baseflux = flux.flux(T, atom.m, P)
    print "Baseflux: ", baseflux

    influx = r**2 * np.pi * baseflux
    print "Flux through first pinhole: ", influx

    dPin = doublePinFrac(r, lcoll)
    print "Fractional flux out of double pinhole:", dPin
    influx2 = influx * dPin
    print "Absolute flux out of double pinhole: %g" %(influx2)

    thmax = np.arctan(2 * r / lcoll)

    u = np.sqrt(2 * kB * T / atom.m)
    print "Most probable velocity: ", u

    C = 1/u**2
    vfraq = 1 - np.exp(-C*v0**2) * (1 + C*v0**2)
    print "Velocity capture fraction", vfraq
    influx3 = vfraq * influx2
    print("Flux after slower: %g" %(influx3))

    thmaxz = np.arctan(r2/(l1+l2))
    print "Pinhole max angle", thmax
    print "Zeeman max angle", thmaxz
    if (thmaxz < thmax):
        print "System is limited by Zeeman tube diameter"
        # fraqDPrevised = np.sin(thmaxz)**2 / 4
        fraqDPrevised = np.sin(thmax)**2 / 4 * 0.17
        # fraqDPrevised = np.sin(thmax)**2 / 4 * (thmaxz/thmax)
        influx4 = vfraq * fraqDPrevised * influx
        print "Pinhole plus slower limited fraction", fraqDPrevised
        print("Revised flux --> %g" %(influx4))
    else:
        influx4 = influx3
