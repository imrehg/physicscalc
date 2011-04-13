#!/usr/bin/python2
from __future__ import division
import numpy as np
import pylab as pl

def mirrors(L, R):
    prop = np.mat([[1, L],[0, 1]])
    mirr = np.mat([[1, 0],[-2/R, 1]])

    total = prop * mirr * prop * mirr
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
    q = 1 / ((1 / R) - (1j * lam / (np.pi * w**2)))
    return (w, R, q)

def getwidths(L, MR, q, lam, n):
    mirror = mirrors(L, MR)
    res = cavity(q, mirror, n)
    widths = np.array([np.sqrt(-1/(np.imag(1/r))*lam/np.pi) for r in res])
    return widths

def basiccompare():
    l = 852e-9
    w0 = 0.5e-3
    n = 20
    w, R, q = gbeam(0, l, w0)

    print q
    pl.figure(figsize=(11.69, 8.27), dpi=100)

    L = 1
    MR = L
    name = "R = L = %g (confocal)" %(L)
    widths = getwidths(L, MR, q, l, n)
    pl.plot(range(1,n+1), widths/w0, 'bx-', label=name)

    L = 0.8
    MR = 1.0
    name = "R = %g,  L = %g" %(MR, L)
    widths = getwidths(L, MR, q, l, n)
    pl.plot(range(1,n+1), widths/w0, 'ro-', label=name)

    L = 1.0
    MR = 0.8
    name = "R = %g,  L = %g" %(MR, L)
    widths = getwidths(L, MR, q, l, n)
    pl.plot(range(1,n+1), widths/w0, 'gd-', label=name)

    pl.xlabel('Round trip number')
    pl.ylabel('Beam waist at output (compared to input)')
    pl.xlim([0, n])

    pl.legend(loc='best')
    pl.show()


def modematch():
    l = 852e-9
    w0 = 0.364e-3
    n = 20


    pl.figure(figsize=(11.69, 8.27), dpi=100)

    zlist = np.linspace(-1, 1, 51)
    minz = 0
    minval = 100
    for z in zlist:
        L = 0.8
        MR = 1.0
        w, R, q = gbeam(z, l, w0)
        name = "z = %g" %(z)
        widths = getwidths(L, MR, q, l, n)/w0
        pl.plot(range(1,n+1), widths, '.-', label=name)
        wrange = max(widths)-min(widths)
        if wrange < minval:
            minval = wrange
            minz = z
    print "Best: %g (%g)" %(minz, minval)

    pl.xlabel('Round trip number')
    pl.ylabel('Beam waist at output (compared to input)')
    pl.xlim([0, n])

    # pl.legend(loc='best')
    pl.show()

def waistsize(R1, R2, d, lam):
    return ((lam / np.pi)**2 * (d * (R1 - d) * (R2 - d) * (R1 + R2 - d)) / (R1 + R2 - 2 * d)**2)**(0.25)

def wfromq(q, lam):
    temp = 1/q
    R = 1/np.real(temp)
    w = np.sqrt(-lam / (np.pi * np.imag(temp)))
    return (w, R)

def qfromw(w0, z, lam):
    zR = np.pi * w0**2 / lam
    R = z * (1 + (zR / z)**2)
    w = w0 * np.sqrt(1 + (z / zR)**2)
    q = 1 / ((1 / R) - (1j * lam / (np.pi * w**2)))
    return q



def findlens(R, d, lam, w0):
    w0c =  waistsize(R, R, L, lam)
    zc = d/2
    zR = np.pi * w0c**2 / lam
    wR = zc * (1 + (zR / zc)**2)
    w = w0 * np.sqrt(1 + (zc / zR)**2)
    q = 1 / ((1 / wR) - (1j * lam / (np.pi * w**2)))
    q2 = (1 * q + 0) / (0.51 * q + 1)
    
    wout, Rout = wfromq(q2, lam)
    return (wout, Rout)

    # return w0c, zR, wR, w, q, q2

# def lenseffect(d, z, f, w0, lam)



if __name__ == "__main__":
    # basiccompare()
    # modematch()
    R = 1
    L = 0.8
    lam = 852e-9
    # print waistsize(R, R, L, 852e-9) * 1e6
    # print (1 - L/R)*(1 - L/R)

    print findlens(R, L, lam, 0.5e-3)
    q = qfromw(5e-4, 1, lam)
    print wfromq(q, lam)
