from __future__ import division
import scipy as sp
from scipy.optimize import leastsq
import pylab as pl

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
    m = sp.dot(matlens(f), m)
    m = sp.dot(matprop(x), m)
    m = sp.dot(matprop(d / (n ** 3)), m)
    m = sp.dot(matprop(z-x-d), m)
    m = sp.dot(matlens(f), m)
    m = sp.dot(matprop(L2), m)
    return m


def abcd1(xx, f, L):
    ''' ABCD matrix of a return trip from the crystal face
    to pane mirror and back'''
    m = matprop(xx)
    m = sp.dot(matlens(f), m)
    m = sp.dot(matprop(2 * L), m)
    m = sp.dot(matlens(f), m)
    m = sp.dot(matprop(xx), m)
    return m


def stab(x, z, L1, L2, f, d, n):
    ''' S = A0 * D0 + B0 * D0, as in Cerullo 1994 Eq. 3
    Required -1 < S < 1 for stability '''
    abcd = abcd0(x, z, L1, L2, f, d, n)
    return abcd[0,0] * abcd[1,1] + abcd[0,1] * abcd[1,0]


# Find stability regions
def err(z, x, L1, L2, f, d, n, upper):
    if upper:
        lim = 1
    else:
        lim = -1
    return stab(x, z, L1, L2, f, d, n) - lim

xstart = 0.041
xfinish = 0.057
zstart = 0.113
zfinish = 0.121

#HMS region
z0 = 0.114
x = (z0 - d) / 2
zup1, success = leastsq(err, z0, args=(x, L1, L2, f, d, n, True))
zdo1, success = leastsq(err, z0, args=(x, L1, L2, f, d, n, False))
pl.plot([xstart, xfinish],[zup1, zup1],'k--')
pl.plot([xstart, xfinish],[zdo1, zdo1],'k--')

#LMS region
z0 = 0.120
zup2, success = leastsq(err, z0, args=(x, L1, L2, f, d, n, True))
zdo2, success = leastsq(err, z0, args=(x, L1, L2, f, d, n, False))
pl.plot([xstart, xfinish],[zup2, zup2],'k--')
pl.plot([xstart, xfinish],[zdo2, zdo2],'k--')
pl.xlabel('x (m)')
pl.ylabel('z (m)')
#pl.show()

def delta1(x, z, L1, L2, f, d, n):
    temp1 = abcd1(x, f, L1)
    alpha1 = temp1[0,1] / (d / n **3) + temp1[0,0]
    temp2 = abcd1(z - x - d, f, L2)
    alpha2 = temp2[0,1] / (d / n **3) + temp2[0,0]
    S = stab(x, z, L1, L2, f, d, n)
    d1 = - (alpha1 + alpha2 * S) / \
        (2 * (alpha1 + alpha2 + 2 * alpha1 * alpha2 * S))
    return d1

x = 0.045
z = 0.115
print "delta1 = %f" %(delta1(x, z, L1, L2, f, d, n))

#npts = 15
#px = sp.linspace(xstart,xfinish,npts)
#pz = sp.linspace(zstart,zfinish,npts)
#PX,PZ = pl.meshgrid(px, pz)
#Z = sp.zeros((npts, npts))
#
#for i in range(0,npts):
#    for j in range(0,npts):
#        Z[i,j] = delta1(PX[i, j], PZ[i, j], L1, L2, f, d, n)
#        if Z[i, j] < 0:
#            Z[i, j] = 0
#
##im2 = pl.imshow(Z, cmap=pl.cm.jet, alpha=.9, interpolation='bilinear')
#
#pl.contour(PX, PZ, Z)
#pl.show()
