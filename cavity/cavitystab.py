#!/usr/bin/python
from __future__ import division
import numpy as np
import pylab as pl

def mirrors(L, R):
    prop = np.mat([[1, L],[0, 1]])
    mirr = np.mat([[1, 0],[-2/R, 1]])

    total = mirr * prop * mirr * prop
    return total

def cavity(q, mirror, n):
    res = []
    qbar = q
    for i in xrange(n):
        temp = mirror * np.matrix([[qbar], [1]])
        qbar = temp[0,0] / temp[1,0]
        res += [qbar]
    return res

def gbeam(z, lam, w0):
    z += 1e-17
    zR = np.pi * w0**2 / lam
    R = z * (1 + (zR / z)**2)
    w = w0 * np.sqrt(1 + (z / zR)**2)
    q = 1 / ((1 / R) - (1j * l / (np.pi * w**2)))
    return (w, R, q)

def getwidths(L, MR, q, n):
    mirror = mirrors(L, MR)
    res = cavity(q, mirror, n)
    widths = np.array([np.sqrt(-1/(np.imag(1/r))*l/np.pi) for r in res])
    return widths

l = 852e-9
w0 = 0.5e-3
n = 20
w, R, q = gbeam(0, l, w0)

print q
pl.figure(figsize=(11.69, 8.27), dpi=100)

L = 0.5
MR = 0.5
name = "R = L = 0.5 (confocal)"
widths = getwidths(L, MR, q, n)
pl.plot(range(1,n+1), widths/w0, 'bx-', label=name)

L = 0.5
MR = 0.8
name = "R = 0.8 ,  L = 0.5"
widths = getwidths(L, MR, q, n)
pl.plot(range(1,n+1), widths/w0, 'ro-', label=name)

L = 0.8
MR = 0.5
name = "R = 0.5 ,  L = 0.8"
widths = getwidths(L, MR, q, n)
pl.plot(range(1,n+1), widths/w0, 'gd-', label=name)

# L = 0.9
# MR = 0.5
# name = "R = 0.5 ,  L = 0.8"
# widths = getwidths(L, MR, q, n)
# pl.plot(range(1,n+1), widths/w0, 'ks-', label=name)

pl.xlabel('Round trip number')
pl.ylabel('Beam waist at output (compared to input)')
pl.xlim([0, n])

pl.legend(loc='best')
pl.show()
