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

def relflux(params):
    atom = params['atom']
    T = params['T']
    eta = params['eta']
    l1, l2, l3 = params['lengths'] # tuple
    rmot = params['eta']
    vf = params['vf']

    v0new = np.sqrt(2 * l2 * atom.aslow * eta + vf**2) # Maximum capture velocity
    u = np.sqrt(2 * kB * T / atom.m) # Most probable velocity
    # Normalize
    v0 = v0new /u
    vf /= u

    collimator = (1, 2*(l1+l2+l3)/rmot)
    Dvrmax = lambda x: rmot/cs.ttime(x, vf, l1, l2, l3, atom.aslow*eta, u)/u # Divergence limit
    Gvrmax = lambda x: rmot/(l1+l2+l3)*x # Geometric limit

    spinhole = (2/collimator[1])**2/(1+(2/collimator[1])**2)
    dpinhole = integ.dblquad(cs.intfunc, 0, np.inf, lambda x: 0, Gvrmax, args=collimator)[0]
    GLim = integ.dblquad(cs.intfunc, 0, v0, lambda x: 0, Gvrmax, args=collimator)[0]
    DLim = integ.dblquad(cs.intfunc, 0, v0, lambda x: 0, Dvrmax, args=collimator)[0]

    return (spinhole, dpinhole, GLim, DLim)

def RbPv(T):
    return 10**(15.88253 - 4529.635 / T + 0.00058663 * T - 2.991381 * np.log10(T))


if __name__ == "__main__":
    #### All parameters
    atom = zs.Rb85()
    T = 80 + 273 # Kelvin
    eta = 0.5
    l1 = 0
    l3 = 200e-3
    rmot = (25.4e-3)/2
    vf = 30
    r = 1e-3 # pinhole radius
    #### No parameters after this
    P = RbPv(T) # vapor pressure

    l2list = np.linspace(0.1, 1.2, 41)

    spin, dpin, gl, dl = [], [], [], []
    for l2 in l2list:
        params = {'atom': atom,
                  'T': T,
                  'eta': eta,
                  'rmot': rmot,
                  'vf': vf,
                  'lengths': (l1, l2, l3),
                  }

        spinhole, dpinhole, GLim, DLim = relflux(params)
        spin += [spinhole]
        dpin += [dpin]
        gl += [GLim]
        dl += [DLim]

    gla = np.array(gl)
    dla = np.array(dl)

    # # Need normalization of gla, dla before this can be used
    # baseflux = flux.flux(T, atom.m, P)
    # influx = r**2 * np.pi * baseflux
    # print baseflux, influx

    fig = pl.figure(num=1, figsize=(11.69, 8.27))
    pl.plot(l2list, gla, 'k--', label='no divergence', linewidth=3)
    pl.plot(l2list, dla, 'r-', label='divergence due to slowing', linewidth=3)
    pl.xlabel('Slower length [m]', fontsize=15)
    pl.ylabel('Flux [AU]', fontsize=15)
    pl.legend(loc='best')
    pl.savefig('designtest.png')
    pl.savefig('designtest.pdf')
    pl.show()
