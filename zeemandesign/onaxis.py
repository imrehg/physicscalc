"""
Zeeman slower design

Basic on-axis calculation.
"""
from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import pylab as pl

##### Units
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
tU = 1/GRb # Time: 1/linewidth
xU = 1/kRb # Distance: 1/wavenumber
# Planck constant rescaled
hU = 6.626e-34 / mU / xU / xU * tU
hbar = hU / 2 / np.pi
# Bohr-magneton rescaled
uB = hU * 1.399e6 * tU
BU = 1/uB # Magnetic field unit

# Dimensionless numbers
mRb = 86.909


# Some variables
k = 1.0
G = 1.0
s0 = 10
det0 = 0
Bb = 1 / BU
Bt = 3 / BU
gP3 = 1.3362
gS1 = 2.002
me = -4
mg = -3
uprime = uB * (gP3*me - gS1*mg)

def bfield(z, z0, Bb, Bt):
    """ Currently non-functional """
    return z * 0
    # return Bb - Bt * np.sqrt(1 - z / z0)

def decelerate(y, t):
    """ Decelartion function """
    x, v = y[0], y[1]
    return [v, -hbar / (2*mRb) * s0 / (1 + s0 + 4*(det0 + v + bfield(x, Ls, Bb, Bt)*uprime/hbar)**2)]

Ls = 0.54 / xU
v0 = 1 / xU * tU

t = sp.linspace(0, 2.9e5, 1000)

# # Btlist = xrange(0, 100, 1000)
# # Btlist = [(0, 0)]
# for B in Btlist:
#     Bb, Bt = B
#     y0 = [0, 100*v0]
#     print y0
#     y = odeint(decelerate, y0, t)

#     # pl.subplot(211)
#     # pl.plot(t, y[:, 0]/Ls)
#     # pl.xlabel('Time')
#     # pl.ylabel('Position')
#     # pl.ylim([0, 1])

#     # pl.subplot(212)
#     # pl.plot(t, y[:, 1])
#     # pl.xlabel('Time')
#     # pl.ylabel('Velocity')

#     # pl.plot(t, det0 + y[:, 1])
#     pl.plot(t, det0 + y[:, 1] + bfield(y[:, 0], Ls, Bb, Bt)*uprime/hbar)
#     # pl.plot(y[:, 0], y[:, 0]/Ls)
#     # pl.xlabel('Time')
#     # pl.ylabel('Velocity')

vlist = range(10, 250, 30)
for v in vlist:
    y0 = [0, v*v0]
    print y0
    y = odeint(decelerate, y0, t)

    pl.subplot(211)
    pl.plot(t*tU, y[:, 0]/Ls, label='%d' %v)
    pl.xlabel('Time (s)')
    pl.ylabel('Position (1 = slower length)')
    # pl.legend(loc='best')
    pl.ylim([0, 1])

    pl.subplot(212)
    pl.plot(t*tU, y[:, 1])
    pl.xlabel('Time (s)')
    pl.ylabel('Velocity')

pl.show()

