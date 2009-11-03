from __future__ import division
import scipy as sp
#import pylab as pl

## Input parameters
# Arm lengths, in (m)
L1 = 0.5
L2 = 1.1

R = 0.1
th2 = 14.5 / 180 * sp.pi
f = R * sp.cos(th2) / 2

d = 0.02
n = 1.76
## End:Input parameters


def matprop(d):
    ''' ABCD matrix of straigth propagation '''
    return sp.array([[1, d],[0, 1]])


def matlens(f):
    ''' ABCD matrix of a thin lens '''
    return sp.array([[1, 0],[-1/f, 1]])


def abcd0(x, z, L1, L2, f, d, n):
    ''' ABCD matrix of a single pass from M1 to M2
    Tangential plane.'''
    m = matprop(L1)
    m = sp.dot(m, matlens(f))
    m = sp.dot(m, matprop(x))
    m = sp.dot(m, matprop(d / (n ** 3)))
    m = sp.dot(m, matprop(z-x-d))
    m = sp.dot(m, matlens(f))
    m = sp.dot(m, matprop(L2))
    return m


def stab(x, z, L1, L2, f, d, n):
    ''' S = A0 * D0 + B0 * D0, as in Cerullo 1994 Eq. 3
    Required -1 < S < 1 for stability '''
    abcd = abcd0(x, z, L1, L2, f, d, n)
    return abcd[0,0] * abcd[1,1] + abcd[0,1] * abcd[1,0]


