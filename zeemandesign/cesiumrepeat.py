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

def partialPinhole(ymax, r, L):
    X = lambda y: 2/np.pi*(np.arccos(y) - y*np.sqrt(1-y**2))*2*y*(2*r*L / (L**2+4*r**2*y**2))**2
    return integ.quad(X, 0, ymax)[0]

def intfunc(vr, vz, r, L):
    """ Double pinhole solution in the cylindrical coordinate case """
    y = L/(2*r)*(vr/vz)
    if y > 1 or y < 0:
        return 0
    X = 2/np.pi*(np.arccos(y) - y*np.sqrt(1-y**2)) # overlapping pinholes case
    intpart = vr * np.exp(-vr**2) * vz * np.exp(-vz**2) * 4 * X # velocity distribution, normalization and double pinhole
    return intpart

def ueintfunc(vr, vz, r1, r2, L):
    y = L/(r1+r2)*(vr/vz)
    if y > 1 or y < 0:
        return 0
    d = (r1+r2)*y
    X = overlap(r1, r2, d) / np.minimum(r1, r2)**2 / np.pi
    intpart = vr * np.exp(-vr**2) * vz * np.exp(-vz**2) * 4 * X # velocity distribution, normalization and double pinhole
    return intpart

def ttime(v, vf, l1, l2, l3, a, u=1):
    """ Transit time """
    if type(v) == type(1.0) or type(v) == type(1):
        v = np.array([v])
    slowindex = np.nonzero(v > vf)[0]
    nonslowindex = np.nonzero(v <= vf)[0]
    tout = np.zeros(len(v))
    
    v *= u
    vf *= u
    vs = v[slowindex]
    lB = (vs**2 - vf**2)/(2*a)
    t1 = l1/vs
    t2 = (l2 - lB)/vs
    t3 = (vs - np.sqrt(vs**2 - 2 * a * lB))/(2*a)
    t4 = l3 /vf
    tout[slowindex] = t1+t2+t3+t4
    tout[nonslowindex] = (l1+l2+l3)/v[nonslowindex]
    return tout


def overlap(r1, r2, d):
    """ 
    Overlapping area of circles

    r1, r2: the two circle's radius
    d: distance of their centre
    """
    rmax = r1 if r1 >= r2 else r2
    rmin = r1 if r1 <= r2 else r2

    if d >= (r1+r2):
        ret = 0
    elif (d+rmin) > rmax:
        c1 = r1*r1*np.arccos( (d*d + r1*r1 - r2*r2) / (2 * d * r1) )
        c2 = r2*r2*np.arccos( (d*d + r2*r2 - r1*r1) / (2 * d * r2) )
        c3 = -0.5 * np.sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2) )
        ret = c1 + c2 + c3
    else:
        ret = np.pi*rmin*rmin
    return ret


if __name__ == "__main__":

    #### All parameters
    atom = zs.Cs133()
    T = 80 + 273 # Kelvin
    P = 3e-5 # vapour pressure in torr
    r = 1e-3 # pinhole radius
    lcoll = 76.2e-3 # collimator length
    eta = 0.5
    v0 = 154
    vf = 42
    # eta = 0.82
    # v0 = 205
    # vf = 50
    r2 = 3.5e-3
    l1 = 200e-3
    l2 = 400e-3
    l3 = 150e-3
    rmot = (25.4e-3)/2
    #### No parameters after this
    P = CsPv(T) # fix vapour pressure from doc
    print "Calculated vapour pressure %.3e torr" %(P)

    # l2new = zs.slowerlength(atom.aslow, eta, v0, vf)
    # print "Slower correct length: %g" %(l2new)
    v0new = np.sqrt(2 * l2 * atom.aslow * eta + vf**2)
    print "Calculated capture velocity: %g m/s" %(v0new)

    u = np.sqrt(2 * kB * T / atom.m)
    print "Most probable velocity: ", u
    v0 = v0new /u
    vf /= u

    baseflux = flux.flux(T, atom.m, P)
    print "Baseflux: ", baseflux

    influx = r**2 * np.pi * baseflux
    print "Flux through first pinhole: ", influx

    # dPin = doublePinFrac(r, lcoll)
    # print "Fractional flux out of double pinhole:", dPin
    DPvrmax = lambda x: 2*r/lcoll*x # example of allowed limit
    dPin = integ.dblquad(intfunc, 0, np.inf, lambda x: 0, DPvrmax, args=(r, lcoll))[0]
    print "Fractional flux out of double pinhole:", dPin

    influx2 = influx * dPin
    print "Absolute flux out of double pinhole: %g" %(influx2)


    DPvrmax = lambda x: 2*r/lcoll*x # example of allowed limit
    vPin = integ.dblquad(intfunc, 0, v0, lambda x: 0, DPvrmax, args=(r, lcoll))[0]
    print "Fractional flux through pinhole and slower:", vPin
    print "Velocity capture fraction", vPin/dPin

    Zvrmax = lambda x: r2/(l1+l2)*x # example of allowed limit
    zLim = integ.dblquad(intfunc, 0, v0, lambda x: 0, Zvrmax, args=(r, lcoll))[0]
    print "Fractional flux limited by Zeeman slower geomtery:", zLim
    print "Geometry fraction", zLim / vPin
    print "=> \"Naive geometric limitation\" %e (%e)" %(zLim, zLim * influx)

    # Geometric limit
    bvrmax = lambda x: 2*r2/(lcoll+l1+l2) * x
    uePin = integ.dblquad(ueintfunc, 0, v0, lambda x: 0, bvrmax, args=(r, r2, lcoll+l1+l2))[0]
    print "=> Unequal pinholes: %e (%e)" %(uePin, uePin*influx)

    #### Experimental:
    print "\nExperiment", "="*10
    eta = 0.82
    vf = 50
    # eta = 0.5
    # vf = 42
    v0new = np.sqrt(2 * l2 * atom.aslow * eta + vf**2)
    vf /= u
    v0 = v0new/u
    print "New max velocity:", v0new

    geovrmax = lambda x: 2*r2/(lcoll+l1+l2) * x # geometry
    bvrmax = lambda x: rmot/ttime(x, vf, l1, l2, l3, atom.aslow*eta, u) # broadening
    bgvrmax = lambda x: np.minimum(r2/(l1+l2) * x, rmot/ttime(x, vf, l1, l2, l3, atom.aslow*eta, u)) # double limitation of geometry and broadening

    args2 = (r, r2, lcoll+l1+l2)
    uePin = integ.dblquad(ueintfunc, 0, v0, lambda x: 0, geovrmax, args=args2)[0]
    uePin2 = integ.dblquad(ueintfunc, 0, v0, lambda x: 0, bvrmax, args=args2)[0]
    uePin3 = integ.dblquad(ueintfunc, 0, v0, lambda x: 0, bgvrmax, args=(r, r, lcoll))[0]
    print "=> Unequal pinholes: %e (%e)" %(uePin, uePin*influx)
    print "=> Unequal pinholes + broadening: %e (%e)" %(uePin2, uePin2*influx)
    print "=> Equal pinholes + geometry + broadening: %e (%e)" %(uePin3, uePin3*influx)

    # # According to the plots, the broadening is calculated correctly
    # v = np.linspace(0, 500, 201)/u
    # tt = ttime(v, vf, l1, l2, l3, atom.aslow*eta, u)
    # pl.plot(v, rmot/(l1+l2+l3)*v, 'r-', label="Geometric limit", linewidth=3)
    # pl.plot(v, rmot/tt, 'k--', label='Broadening limit', linewidth=3)
    # pl.legend(loc='best')
    # pl.xlabel('v_z')
    # pl.ylabel('v_r')
    # pl.show()

    eta = 0.82
    vf = 50
    vf /= u

    sl = np.linspace(0.2, 1.2, 10)
    cfraqb = np.array([])
    for l2 in sl:
        print "Slower length: ", l2
        v0new = np.sqrt(2 * l2 * atom.aslow * eta + vf**2)

        v0 = v0new/u

        args2 = (r, r2, lcoll+l1+l2)
        bvrmax = lambda x: rmot/ttime(x, vf, l1, l2, l3, atom.aslow*eta, u) # broadening

        uePin2 = integ.dblquad(ueintfunc, 0, v0, lambda x: 0, bvrmax, args=args2)[0]
        cfraqb = np.append(cfraqb, uePin2)

    pl.plot(sl, cfraqb)
    pl.show()
