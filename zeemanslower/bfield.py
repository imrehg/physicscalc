#!/usr/bin/env python
"""
Simulation and data fitting to design a Zeeman slower.
Follows in the footsteps of Bell et.al., Rev. Sci. Inst 81 (2010) 013105
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ
from mpl_toolkits.mplot3d import Axes3D

def bideal(z, params):
    Ba = 0.5
    hbar = 1
    k = -1
    mu = 1
    v0 = params['v0']
    nu = 1
    a = 1
    return Ba - hbar*k/mu*np.sqrt(v0**2 - 2 * nu * a * z)

def pos(p, params):
    """ parameteric coil """
    c = params['c']
    R = params['R']
    theta = np.polyval(c[0:7], p)
    x = R * np.cos(theta)
    y = R * np.sin(theta)
    z = np.polyval(c[7:9], p)
    return zip(x, y, z)

def deriv(p, params):
    """ dr/dp """
    c = params['c']
    # The coefficients of the polynomial after derivation
    cnew = [c[i]*(6-i) for i in range(0,6)]
    R = params['R']
    theta = np.polyval(c[0:7], p)
    cpderiv = np.polyval(cnew, p)
    dx = -R * np.sin(theta) * cpderiv
    dy = R * np.cos(theta) * cpderiv
    dz = c[8] + p*0
    return zip(dx, dy, dz)

def bzfield(p, params, z):
    cord = params['cord']
    cp = params['cp']
    cn = params['cn']
    cpnew = [cp[i]*(cord-i) for i in range(0,cord)]
    cnnew = [cn[i]*(cord-i) for i in range(0,cord)]
    R = params['R']
    upperp = np.polyval(cpnew, p)
    lowerp = (R**2 + (np.polyval(cp[cord+1:], p)-z)**2)**(3/2)
    uppern = np.polyval(cnnew, p)
    lowern = (R**2 + (np.polyval(cn[cord+1:], p)-z)**2)**(3/2)
    return upperp/lowerp + uppern/lowern


## Coil parameters
cp = [0, 1.6e-3, 0, -0.3, 8.7, 3, 0, -0.1, 0.58]
cn = [-0.0012, 0.0035, 0.02, -0.23, -1.257, -2, 3.14, 0.043, 0.62]
R = 0.04
cord = 6

params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord}
zl = np.linspace(-0.2, 1.1, 201)
res =  [integ.quad(bzfield, 0, 2*np.pi, args=(params, zz))[0] for zz in zl]
pl.plot(zl, res)
pl.xlim([zl[0], zl[-1]])
pl.xlabel('Length')
pl.ylabel('"Field"')
pl.show()
