from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ

import zeemanslower as zs

def intfunc(vr, vz, r, L):
    """ Double pinhole solution in the cylindrical coordinate case """
    y = L/(2*r)*(vr/vz)
    if y > 1 or y < 0:
        return 0
    X = 2/np.pi*(np.arccos(y) - y*np.sqrt(1-y**2)) # overlapping pinholes case
    intpart = vr * np.exp(-vr**2) * vz * np.exp(-vz**2) * 4 * X # velocity distribution, normalization and double pinhole
    return intpart

def ttime(v, vf, l1, l2, l3, a):
    """ Transit time """
    # fixme: whenever v < vf, the transit time is just (l1+l2+l3)/v
    if type(v) == type(1.0) or type(v) == type(1):
        v = np.array([v])
    slowindex = np.nonzero(v > vf)[0]
    nonslowindex = np.nonzero(v <= vf)[0]
    tout = np.zeros(len(v))
    
    vs = v[slowindex]
    lB = (vs**2 - vf**2)/(2*a)
    t1 = l1/vs
    t2 = (l2 - lB)/vs
    t3 = (vs - np.sqrt(vs**2 - 2 * a * lB))/(2*a)
    t4 = l3 /vf
    tout[slowindex] = t1+t2+t3+t4
    tout[nonslowindex] = (l1+l2+l3)/v[nonslowindex]
    return tout

if __name__ == "__main__":

    atom = zs.Cs133()
    eta = 0.5
    v0 = 154
    vf = 42
    r2 = 3.5e-3
    l1 = 200e-3
    l2 = 400e-3
    l3 = 150e-3
    rmot = (25.4e-3)/2

    r2 = 0.1
    # Double pinhole dimensions
    r = 1
    L = 76.4

    W = 2*r/L
    X = W**2/(1+W**2)
    print "Single pinhole capture fraction", X

    ### Can adjust these
    vrmax = lambda x: np.inf # Normal limit
    vrmax = lambda x: 2*r/L*x / 2 # example of half the allowed limit
    vrmax = lambda x: 2*r/L*x / 2.8925 # example of half the allowed limit
    vrmax = lambda x: r2/ttime(x, vf, l1, l2, l3, atom.aslow*eta) ### !!!!!! Gives incorrect result at the moment
    ###

    translim = integ.dblquad(intfunc, 0, np.inf, lambda x: 0, vrmax, args=(r, L))[0]
    print "Double pinhole capture fraction", translim

    print "Ratio", X/translim

    vz = np.linspace(1, 100, 200)
    pl.plot(vz, r2/ttime(vz, vf, l1, l2, l3, atom.aslow*eta))
    pl.plot(vz, 2*r/L*vz)
    pl.show()
