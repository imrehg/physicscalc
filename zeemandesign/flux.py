from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ

import pinholes as ph
import zeemanslower as zs

#### Dimensioned
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 * 2 * np.pi # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
# mRb = 86.909 * mU
mRb = 85 * mU  # Rb85
h = 6.626e-34
hbar = h / 2 / np.pi
uprimehbar = 1.399e10 * 2 * np.pi
uBr = 9.27e-24
kB = 1.3806488e-23

def captureFraction(r1, r2, l, sl, toplot=False):
    r_spot = (25.4e-3)/2  # one inch diameter MOT beam
    th_lim = np.arctan(r_spot/sl)
    # th_max = np.pi/2
    th_max = np.arctan((r1+r2)/l)

    th = np.linspace(0, th_max, 201)
    totalout = integ.trapz(angledist(th, r1, r2, l), th)
    th2 = np.linspace(0, th_lim, 201)
    captureout = integ.trapz(angledist(th2, r1, r2, l), th2)
    if toplot:
        pl.figure()
        pl.plot(th, angledist(th, r1, r2, l), 'k-')
        pl.fill_between(th2, 0, angledist(th2, r1, r2, l))
        pl.xlim([0, np.min((th_lim*1.5, th_max))])
    return (totalout, captureout)


def angledist(th, r1, r2, l):
    return np.array([ph.angleArea(thi, r1, r2, l)*np.cos(thi) for thi in th]) / (r1**2*np.pi)

def vapourPressure(T):
    """ Liquid phase vapour pressure for Rb-85 in torr"""
    return 10**(2.881 + 4.312 - 4040 / T)

def numberDensity(T, P=None):
    """ Number density """
    if P is None:
        P = vapourPressure(T)
    return 9.66e24*P/T # atoms / m^3

def flux(T, m, P=None):
    """ Flux from a reservoir """
    n0 = numberDensity(T, P)
    vavg = np.sqrt(8*kB*T/(np.pi*m))
    return 0.25 * n0 * vavg

def influx(r, T, m):
    """ Flux through a small aperture (atoms/s)"""
    f = flux(T, m)
    return f * r * r * np.pi

def intMaxwell(maxv, T):
    u = np.sqrt(2 * kB * T / mRb)
    vdist = lambda v: v**3 * np.exp(-v*v/u/u)
    total = integ.quad(vdist, 0, maxv)[0]
    return total

if __name__ == "__main__":

    # T = 273 + 80
    # print vapourPressure(T)
    # print numberDensity(T)
    # print flux(T, mRb)
    # print influx(1e-3, T, mRb)

    r1 = 1e-3
    r2 = r1
    l = 100e-3
    th_max = np.pi/2

    atom = zs.Rb85()
    eta, vf = 0.7, 20
    v0list = np.linspace(200, 550, 41)

    sllist = []
    for v0 in v0list:
        sllist += [zs.slowerlength(atom.aslow, eta, v0, vf)]

    pl.figure(figsize=(11.69, 8.27))
    pl.subplot(221)
    pl.plot(v0list, sllist)
    pl.xlabel('Maximum capture velocity (m/s)')
    pl.ylabel('Slower length (m)')
    pl.title(r'$v_f$=%d m/s, $\eta$=%g' %(vf, eta))

    T = 273+80
    vnorm = intMaxwell(1200, T)
    vfrac = []
    for v0 in v0list:
        vfrac += [intMaxwell(v0, T)/vnorm]
    vfrac = np.array(vfrac)

    pl.subplot(222)
    pl.plot(sllist, vfrac)
    pl.xlabel('Slower length (m)')
    pl.ylabel('Velocity capture fraction')
    pl.title('T=%d K ; r1=%g r2=%g l=%g (mm)' %(T, r1*1e3, r2*1e3, l*1e3))

    dl1 = 0.2
    dl2 = 0.2
    sl = 0.4
    l_tot = dl1 + sl + dl2
    cfrac = []
    outfrac = []
    for sl in sllist:
        l_tot = dl1 + sl + dl2
        tot, capt = captureFraction(r1, r2, l, l_tot, False)
        outfrac += [capt]
        cfrac += [capt/tot]
    outfrac = np.array(outfrac)
    cfrac = np.array(cfrac)

    pl.subplot(223)
    pl.plot(sllist, cfrac)
    pl.xlabel('Slower length (m) (+%g+%g extra)' %(dl1, dl2))
    pl.ylabel('Geometrical capture fraction')

    # pl.subplot(224)
    # pl.plot(sllist, cfrac*vfrac)
    # pl.xlabel('Slower length (m)')
    # pl.ylabel('Capture total fraction')

    ax1 = pl.subplot(224)
    captfract = outfrac*vfrac
    ax1.plot(sllist, captfract)
    pl.xlabel('Slower length (m)')
    ax1.set_ylabel('Total flux capture fraction')
    ax1.set_ylim([captfract[0], np.max(captfract)*1.1])

    fl = influx(r1, T, mRb)
    ax2 = ax1.twinx()
    totalflux = fl*outfrac*vfrac/1e12
    ax2.plot(sllist, totalflux)
    ax2.set_ylabel('Total flux (10^12 atoms/s)')
    ax2.set_ylim([totalflux[0], np.max(totalflux)*1.1])

    pl.show()

    
