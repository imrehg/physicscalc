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
import scipy.odr as odr

## Constants:
mu = 1.26e-6*1.5e-3
hbar = 1.05457148e-34

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
    cord = params['cord']
    cp = params['cp']
    cn = params['cn']
    R = params['R']
    thetap = np.polyval(cp[0:(cord+1)], p)
    xp = R * np.cos(thetap)
    yp = R * np.sin(thetap)
    zp = np.polyval(cp[(cord+1):(cord+3)], p)
    thetan = np.polyval(cn[0:(cord+1)], p)
    xn = R * np.cos(thetan)
    yn = R * np.sin(thetan)
    zn = np.polyval(cn[(cord+1):(cord+3)], p)
    return (xp, yp, zp, xn, yn, zn)

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
    I = params['I']
    upperp = np.polyval(cpnew, p)
    lowerp = (R**2 + (np.polyval(cp[cord+1:], p)-z)**2)**(3/2)
    uppern = np.polyval(cnnew, p)
    lowern = (R**2 + (np.polyval(cn[cord+1:], p)-z)**2)**(3/2)
    return mu*I/(4*np.pi)*(upperp/lowerp + uppern/lowern)

def plotcoils(params):
    prange = np.linspace(0, 2*np.pi, 40001)
    xp, yp, zp, xn, yn, zn = pos(prange, params)
    pl.figure()
    pl.plot(zp, yp, label="Positive coil")
    pl.plot(zn, yn, label="Negative coil")
    pl.xlabel('z')
    pl.ylabel('y')
    pl.title('Coils - side view')
    pl.legend(loc="best")

def fitfunction(c, z):
    ## Coil parameters
    cp = c[0:9]
    cn = c[9:18]
    # I = c[18]
    I = 110
    cord = 6
    ## Input parameters
    R = 0.0383
    params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
    # print c
    if type(z) == type(1):
        z = np.array([z])
    prange = np.linspace(0, 2*np.pi, 4001)
    res = np.array([integ.trapz(bzfield(prange, params, zz), prange) for zz in z])
    return res

def bideal(z):
    C1 = 0.015
    C2 = 0.03
    C3 = 1
    C4 = 1.13
    return -(C1 - C2 * np.sqrt(C3**2  - C4 * z))    

def dofit(c0, npoints):
    z = np.linspace(0, 0.85, npoints)
    goal = bideal(z)
    data = odr.Data(z, goal)
    model = odr.Model(fitfunction)
    fit = odr.ODR(data, model, c0)
    fit.set_job(fit_type=2)
    # fit.set_iprint(iter_step=9)
    out = fit.run()
    return out

## Coil parameters
cp = [0, 1.58e-3, 0, -0.302, 8.698, 1.999, 0, -0.106, 0.581]
cn = [-0.0012, 0.0035, 0.02, -0.23, -1.257, -2, 3.14, 0.043, 0.62]
cord = 6
## Input parameters
R = 0.0383
I = 110
params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}

# prange = np.linspace(0, 2*np.pi, 101)
# xp, yp, zp, xn, yn, zn = pos(prange, params)
# print xp, zn

# plotcoils(params)
# pl.show()


c = np.array(cp+cn)
c = [ -1.33845386e-04,   2.77900966e-03,   1.58008750e-02,  -6.24862682e-01,
   1.21973039e+01,   3.10753454e+00,   0.00000000e+00,  -1.33363401e-01,
   6.46170422e-01,  -1.99539130e-02,   9.22830584e-02,  -3.17120729e-01,
   5.87834745e-01,  -5.18028816e+00,   9.91083098e-02,   3.14000000e+00,
   7.37785810e-02,   6.64226717e-01]


# cp = c[0:9]
# cn = c[9:18]
# params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
# plotcoils(params)
# pl.show()

nump = 41
zl = np.linspace(0, 0.85, nump)
out = dofit(c, 21)
out.pprint()
fit = fitfunction(out.beta, zl)
# fit = fitfunction(c, zl)
pl.subplot(211)
pl.plot(zl, fit)
pl.plot(zl, bideal(zl))
pl.subplot(212)
pl.plot(zl, fit-bideal(zl))
pl.show()

# res =  [integ.quad(bzfield, 0, 2*np.pi, args=(params, zz))[0]*1e3 for zz in zl]
# pl.plot(zl, res)
# pl.xlim([zl[0], zl[-1]])
# pl.xlabel('Length')
# pl.ylabel('"Field"')
# pl.show()
