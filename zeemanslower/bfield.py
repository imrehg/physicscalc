#!/usr/bin/env python
"""
Simulation and data fitting to design a Zeeman slower.
Follows in the footsteps of Bell et.al., Rev. Sci. Inst 81 (2010) 013105
"""

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

# z = np.linspace(0, 1, 201)
# params = {'v0': 1}
# pl.plot(z, bideal(z, params))
# pl.show()


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

c = [1e-4, 3, 8.7, -0.3, 1e-4, 1.6e-3, 1e-4][::-1] + [-0.1, 0.58]
R = 1
params = {'R': R, 'c': c}
prange = np.linspace(0, 2*np.pi, 3001)
helix = pos(prange, params)
x, y, z = zip(*helix)
dhelix = deriv(prange, params)
dx, dy, dz = zip(*dhelix)

# fig = pl.figure()
# ax = Axes3D(fig)
# ax.plot(x, y, z, label='Helix')

dx1 = np.array(np.diff(x))
dy1 = np.array(np.diff(y))

dp1 = np.array(np.diff(prange))
pl.subplot(211)
pl.plot(prange[0:-1], dx1/dp1)
pl.plot(prange, dx)
pl.subplot(212)
pl.plot(prange[0:-1], dy1/dp1)
pl.plot(prange, dy)

pl.show()
