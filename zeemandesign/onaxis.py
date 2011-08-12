"""
Zeeman slower design

Basic on-axis calculation.
"""
from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import scipy.integrate as integ
from scipy.interpolate import interp1d
import pylab as pl

# ##### Units
# GRb = 38.117e6 # 2pi x 6.07MHz
# kRb = 1281654.9389 * 2 * np.pi # 1/m
# mU = 1.667e-27 # Mass: Atomic Mass Unit
# tU = 1/GRb # Time: 1/linewidth
# xU = 1/kRb # Distance: 1/wavenumber
# # Planck constant rescaled
# hU = 6.626e-34 / mU / xU / xU * tU
# hbar = hU / 2 / np.pi
# # Bohr-magneton rescaled
# uB = hU * 1.399e6 * tU
# BU = 1/uB # Magnetic field unit

# # Dimensionless numbers
# mRb = 86.909


# # Some variables
# k = 1.0
# G = 1.0
# s0 = 10
# det0 = -2
# Bb = 1 / BU
# Bt = 3 / BU
# gP3 = 1.3362
# gS1 = 2.002
# me = -4
# mg = -3
# uprime = uB * (gP3*me - gS1*mg)


#### Dimensioned
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 * 2 * np.pi # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
# mRb = 86.909 * mU
mRb = 85 * mU  # Rb85
h = 6.626e-34
hbar = h / 2 / np.pi
uprimehbar = 1.399e10 * 2 * np.pi


# def bfield(z, z0, Bb, Bt):
#     """ Currently non-functional """
#     return z * 0
#     # return Bb - Bt * np.sqrt(1 - z / z0)

def bfzero(z):
    return z*0

def decelerate(y, t, bfield, s0, det0):
    """ Decelartion function """
    x, v = y[0], y[1]
    return [v, -hbar * kRb * GRb / (2*mRb) * s0 / (1 + s0 + 4*(det0 + kRb * v - bfield(x)*uprimehbar)**2/GRb**2)]
    # var = [v, -hbar * kRb / (2*mRb) * s0 / (1 + s0 + 4*(det0 + kRb * v - bfield(x)*uprimehbar)**2/GRb**2)]
    # var = [v, -hbar * kRb * GRb / (2*mRb)]
    # return var

def simulate(v0, bfield):
    s0 = 10
    det0 = -283.125e6 * 2 * np.pi
    det0 = -260e6 * 2 * np.pi
    
    print (det0 + kRb*v0 - bfield(0)*uprimehbar)/GRb
    print (det0 + kRb * v0 - bfield(0)*uprimehbar)/GRb

    nu = 0.7
    z0 = mRb * (365**2 - 20**2)/(nu * hbar *kRb * GRb)
    print "Slower length: ", z0
    
    t = sp.linspace(0, 1e-2, 400)
    y0 = [-0.5, v0]

    y, info = odeint(decelerate, y0, t, args=(bfield, s0, det0), full_output=True)
    xi = y[:, 0]
    vi = y[:, 1]
    # # pl.plot(xi, (det0 + kRb * vi - bfield(xi)*uprimehbar) )
    # # y = odeint(decelerate, y0, t, args=(bfzero, s0, det0))
    # xi = y[:, 0]
    # vi = y[:, 1]

    pl.figure(figsize=(11.69, 8.27))
    Bn = 90
    pl.suptitle('v0 = %.1f, s0 = %.2f, detuning = %.3fMHz, In = %dA' %(v0, s0, det0/2/np.pi/1e6, Bn))
    # pl.plot(xi, (det0 + kRb * vi - bfield(xi)*uprimehbar))
    pl.subplot(221)
    pl.plot(xi, vi)
    pl.xlabel('x [m]')
    pl.ylabel('v [m/s]')

    pl.subplot(222)
    pl.plot(xi, (det0 + kRb * vi - bfield(xi)*uprimehbar)/GRb)
    pl.xlabel('x')
    pl.ylabel('detuning [1/G]')


    ai = [decelerate((x, v), 0, bfield, s0, det0)[1] for x, v in zip(xi, vi)]
    pl.subplot(223)
    # pl.plot(xi, ai)
    # pl.xlabel('x')
    pl.plot(vi, ai)
    pl.xlabel('v [m/s]')
    pl.ylabel('a [m/s^2]')

    pl.subplot(224)
    xt = np.linspace(-0.3, 1.2, 100)
    # pl.plot(xt, bfield(xt), '.-')
    # pl.plot([xt[0], xt[-1]], [0, 0])

    # pl.plot(xt, bfield(xt)*uprimehbar/2/np.pi/1e6, '.')
    # pl.plot(xi, (det0 + kRb * vi)/2/np.pi/1e6, '-')
    pl.plot(xt, bfield(xt)*1e3, '.')
    pl.plot(xi, (det0 + kRb * vi)/uprimehbar*1e3, '-')
    pl.plot([xt[0], xt[-1]], [0, 0])

    # pl.plot(vi, (det0 + kRb * vi - bfield(xi)*uprimehbar)/GRb+1/2)
    # pl.xlim([-50, 50])

    pl.xlabel('x [m]')
    pl.ylabel('B(x)/Detuning as magnetic field [mT]')


def fig10(field):
    s0 = 10
    det0 = -260e6 * 2 * np.pi
    t = sp.linspace(0, 3e-2, 400)
    v0 = 365
    y0 = [-0.5, v0]

    z, Bp, Bn = field['z'], field['pos'], field['neg']
    Ip = 110
    Inlist = range(0, 95, 3)
    terminal = 1
    vval = []
    bval = []
    for In in Inlist:
        print In
        Bz = interp1d(z, Bp*Ip+Bn*In, bounds_error=False, fill_value=0)
        y, info = odeint(decelerate, y0, t, args=(Bz, s0, det0), full_output=True)
        speed = interp1d(y[:,0], y[:,1])
        vval += [speed(terminal)]
        bval += [-min(Bz(np.linspace(0.5, 0.9, 201)))*1e3]
    pl.figure(figsize=(11.69, 8.27))
    pl.plot(bval, vval, '.-')
    pl.xlabel('Closing field maximum [mT]')
    pl.ylabel('Velocity [m/s]')
    pl.xticks([0, 5, 10, 15])
    pl.yticks([0, 50, 100, 150, 200])
    pl.xlim([0, 16])
    pl.ylim([0, 210])


# detu = -260e6 * 2*np.pi * tU
# v = 365 / xU * tU
# print detu, v

# # print bfinter(0)*uprime/hbar
# z = np.linspace(-0.3, 1.2, 401)
# field = (bfinter(z) * uprime / hbar) / tU
# pl.plot(z, field)
# pl.plot(0, (detu+k*v)/tU, 'x')
# pl.ylabel('detuning [1/s]')
# pl.show()

# b = 0.015 * 1e4 * uprime / hbar
# print b



# Ls = 1 / xU
# v0 = 1 / xU * tU

# t = sp.linspace(0, 2e6, 100000)

# # # Btlist = xrange(0, 100, 1000)
# # # Btlist = [(0, 0)]
# # for B in Btlist:
# #     Bb, Bt = B
# #     y0 = [0, 100*v0]
# #     print y0
# #     y = odeint(decelerate, y0, t)

# #     # pl.subplot(211)
# #     # pl.plot(t, y[:, 0]/Ls)
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Position')
# #     # pl.ylim([0, 1])

# #     # pl.subplot(212)
# #     # pl.plot(t, y[:, 1])
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Velocity')

# #     # pl.plot(t, det0 + y[:, 1])
# #     pl.plot(t, det0 + y[:, 1] + bfield(y[:, 0], Ls, Bb, Bt)*uprime/hbar)
# #     # pl.plot(y[:, 0], y[:, 0]/Ls)
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Velocity')

# vlist = range(10, 200, 10)
# vfinal = []
# for v in vlist:
#     y0 = [0, v*v0]
#     y = odeint(decelerate, y0, t)

#     xi = y[:, 0] / Ls
#     vi = y[:, 1]

#     if xi[-1] < 1:
#         vfinal += [0]
#     else:
#         vfinal += [vi[xi <= 1][-1]]

#     # # filter negative direction
#     # xi = xi[vi > 0]
#     # t = t[vi > 0]
#     # vi = vi[vi > 0]

#     # # filter slower length
#     # t = t[xi < 1]
#     # vi = vi[xi < 1]
#     # xi = xi[xi < 1]

#     # filter negative direction
#     # xi = xi[vi > 0]
#     # t = t[vi > 0]
#     # vi = vi[vi > 0]
#     # pl.plot(t[vi>0], xi[vi>0])

#     # pl.subplot(211)
#     # # pl.plot(t*tU, xi)
#     # # pl.plot(xi)
#     # pl.xlabel('Time (s)')
#     # pl.ylabel('Position (1 = slower length)')
#     # # pl.legend(loc='best')
#     # # pl.ylim([0, 1])

#     # pl.subplot(212)
#     # pl.plot(t*tU, vi)
#     # pl.xlabel('Time (s)')
#     # pl.ylabel('Velocity')
#     pl.plot(xi[vi>0], det0 + vi[vi>0])
#     pl.xlim([0, 1])

# # pl.plot(vlist, vfinal)
# pl.show()

if __name__ == "__main__":
    """ Basic calculation """
    bfieldfile = 'coilfield.npz'
    # z, B = np.loadtxt(bfieldfile, unpack=True)
    field = np.load(bfieldfile)
    z, Bp, Bn = field['z'], field['pos'], field['neg']
    Ip = 110
    In = 90
    Bz = interp1d(z, Bp*Ip+Bn*In, bounds_error=False, fill_value=0)

    z = np.linspace(0, 1, 201)
    detu = 1.399e4 * Bz(z)
    # pl.plot(z, detu)
    # pl.show()
    v0 = 365
    # v0 = 280
    simulate(v0, Bz)


    fig10(field)
    pl.show()
