from __future__ import division
import numpy as np
import pylab as pl

hbar = 1.05457148e-34
uB = hbar * 2 * np.pi * 1.399e6

B = 150

## This is for RB87
gP3 = 1.3362
gS1 = 2.002
me = -4
mg = -3
uprime = uB * (gP3*me - gS1*mg)

detu = uprime * B / hbar
print detu / 2 / np.pi / 1e6


##### Units
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 * 2 * np.pi # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
tU = 1/GRb # Time: 1/linewidth
xU = 1/kRb # Distance: 1/wavenumber

# Dimensionless detuning
# print detu * tU

v0 = 365
k = 1281655*2*np.pi

# dd = -2*np.pi*265e6 + k*v0
# print dd/ 2/ np.pi / 1e6


# print hbar*k/uprime*v0

a = -v0
t = np.linspace(0, 1, 101)
vt = v0 + a*t
xt = v0*t + a*t**2/2
dt = (-2*np.pi*260e6 + k*vt)/2/np.pi/1e6

b0 = -285
dd = (b0 + hbar*k/uprime*np.sqrt(v0**2 + 2*a*xt))*uprime/hbar/2/np.pi/1e6
print dd
dfield = dt/uprime*hbar*2*np.pi*1e6
pl.plot(xt, dfield)
# pl.plot(0, detu/2/np.pi/1e6, 'x')
# pl.plot(xt, dd, '.-')
# pl.plot(xt, v0**2 + 2*a*xt)
print dfield[-1]
print dfield[-1]/-150
# print "field:", 2*np.pi*260e6/uprime*hbar


# def gf(gi, gj, i, j, f):
#     return gj*(f*(f+1) - i*(i+1) + j*(j+1))/(2*f*(f+1)) + gi*(f*(f+1) + i*(i+1) - j*(j+1))/(2*f*(f+1))

# # P3/2
# ge = gf(-0.0003, 1.3362, 5/2, 3/2, 4)
# gg = gf(-0.0003, 2.002, 5/2, 1/2, 3)
# print ge*4 - gg*3

print "="*10

okb = 1.38e-23
om = 86.9 * 1.660e-27
oT = 273+80
ovp = np.sqrt(2 * okb * oT / om)
print "Most probable temperature:", ovp

# TempU = hbar*hbar*kRb*kRb / okb / om
TempU = (hbar*kRb)**2 / okb / om
EU = mU * xU * xU / tU / tU
# # print

xT = oT / TempU
xkb = okb / EU * TempU
xm = 86.9
xvp = np.sqrt(2 * xkb * xT / xm)
print xvp, xvp * xU / tU
print ovp / xU * tU, ovp
print xvp / (ovp / xU * tU)
print (xvp * xU / tU) / ovp
