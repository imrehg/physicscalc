from __future__ import division
import numpy as np
import pylab as pl

import layeroptimize as lo
import zeemanslower as zs
import sys

from scipy.integrate import odeint
from scipy.optimize import fmin
import time


fieldfile = "16_AWG12Coat_10.npz"
atom = zs.Rb87()

sim = np.load(fieldfile)['simulation'][()]
wire = sim['wire']
setup = sim['setup']
eta = sim['eta']
v0 = sim['v0']
vf = sim['vf']
looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']

eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']
sl = zs.slowerlength(atom.aslow, eta, v0, vf)
z = np.array([0])
bfield = zs.bideal(atom, z, eta, v0, vf, detu)
z = np.append([-10*R], np.append(np.linspace(0, sl, 61), [sl+10*R]))
ze = np.linspace(z[0], z[-1], 201)

ratio = 1
filename = 'trajectory01.npz'
simed = True
if simed:
    setup['csign'] = [s*ratio if s < 0 else s for s in setup['csign']]
    # Get the needed current to run the coils:
    bf = lo.fieldcalc(ze, setup)
    B0 = zs.bideal(atom, 0, eta, v0, vf, detu)
    I = B0[0] / lo.fieldcalc(0, setup)[0]
    print "Current:", I
    bf *= I
    field = lambda z: lo.fieldcalc(z, setup)*I
else:
    bf = zs.bideal(atom, ze, eta, v0, vf, detu)
    field = lambda z: zs.bideal(atom, z, eta, v0, vf, detu)

k = atom.k
s = 0.9
detuw = -2*np.pi*detu*1e6

# vfield = (zs.uprimehbar*field(ze)-detuw)/k

data = np.load(filename)
# v = data['arr_0']
y = data['arr_1']
n = 0
z = y[n][:, 0]
v = y[n][:, 1]

vfield = (zs.uprimehbar*field(z)-detuw)/k

## http://adorio-research.org/wordpress/?p=4493
from math import *
def polyCubicRoots(a,b, c):
    print "input=", a,b,c
    aby3 = a / 3.0
    p = b - a*aby3
    q = (2*aby3**2- b)*(aby3) + c
    X = (p/3.0)**3
    Y = (q/2.0)**2
    Q = X + Y
    D = 18*a*b*c-4*a**3*c+a**2*b**2-4*b**3-27*c**2
    print "Q=", Q, D
    if Q >= 0:
       sqQ = sqrt(Q)
       A0 = -q/2.0 + sqQ
       A = A0**(1/3.0) if A0 > 0 else -(-A0)**(1/3.0)
       B0 = -q/2.0 - sqQ
       B = B0**(1/3.0) if B0 > 0 else -(-B0)**(1/3.0)
       r1 = A + B - aby3
       re = -(A+B)/2.0-aby3
       im = sqrt(3.0)/2.0*(A-B)
       r2 = re+im*1j
       r3 = re-im*1j
    else:
       # This part has been tested.
       p3by27= sqrt(-p**3/27.0)
       costheta = -q/2.0/ p3by27
       print "@@@ costheta=", costheta
       alpha = acos(costheta)
       mag = 2 * sqrt(-p/3.0)
       alphaby3 = alpha/3.0
       r1 = mag  * cos(alphaby3) - aby3
       r2 = -mag * cos(alphaby3+ pi/3)-aby3
       r3 = -mag * cos(alphaby3- pi/3) -aby3
    return r1, r2, r3

def cubic(B, Bz):
    x = detuw - zs.uprimehbar * B
    y = zs.uprimehbar * Bz
    a = 4*k**2/atom.G**2
    b = 8*k*x/atom.G**2
    c = (1 + s + 4*x**2/atom.G**2)
    d = zs.hbar*k**2*atom.G*s/(2*atom.m*y)

    b /= a
    c /= a
    d /= a
    a = 1

    polyCubicRoots(b, c, d)

    return polyCubicRoots(b, c, d)

zerod = lambda v, B, Bz: zs.hbar * k**2 * atom.G / (v * 2*atom.m)*s/(1+s+4/atom.G**2 * (detuw + k*v - zs.uprimehbar*B)**2) + zs.uprimehbar*Bz

res = []
vals = []
nums = range(len(z))
# nums = [26]
for num in nums:
    zx = z[num]
    print ">>> ", num
    if zx < 0 or zx > 0.58:
        res += [np.nan]
        vals += [0]
        continue
    vx = v[num]
    B = field([zx])[0]
    vx0 = (zs.uprimehbar*B - detuw)/k
    vx = vx0
    zstep = 0.0001
    Bz = (np.diff(field([zx-zstep/2, zx+zstep/2]))/zstep)[0]
    # print zx, B, Bz, vx

    # zdd =  zerod(vx, B, Bz)
    # print zerod(vx*1.01, B, Bz)
    # fminfunc = lambda v, B, Bz: abs(zerod(v, B, Bz))
    # vxnew = fmin(fminfunc, vx, args=(B, Bz))

    try:
        vxpos = cubic(B, Bz)
    except (RuntimeError):
        res += [np.nan]
        continue
    vxnew = 0
    for vval in vxpos:
        # print vxpos, vval
        if np.imag(vval) != 0:
            # We have complex roots, that means we jumped the shark
            vxnew = np.nan
            break
        try:
            if vval < vx0 and vval > vxnew:
                vxnew = vval
        except (ValueError):
            vxnew = np.nan
            break
    # vals += [zerod(vxnew, B, Bz)]
    print "=>>", vxnew
    res += [vxnew]
res = np.array(res)


dres = (detuw + k*res - zs.uprimehbar*field(z))/(atom.G/2)
print res

# print zerod(vxnew, B, Bz)
# print cubic(vxnew, B, Bz)
# pl.figure()
# vt = np.linspace(0, vxnew*2, 101)
# pl.plot(vt, cubic(vt, B, Bz))

# pl.figure()
# pl.plot(z, v)
# pl.plot(z, vfield)
# pl.plot(z, res, '--')
# pl.xlim([0, 0.6])

pl.figure()
pl.plot(z, vfield)
pl.plot(z, v)
pl.plot(z, res, 'k-', linewidth=2)
pl.xlim([0, 0.6])


pl.figure()
pl.plot(z, vfield)
pl.plot(z, v)
pl.plot(z, res, 'k-', linewidth=2)
pl.xlim([0, 0.6])


pl.subplot(2,2,1)
pl.plot(z, vfield)
pl.plot(z, v)
pl.plot(z, res, 'k-', linewidth=2)
pl.xlim([0, 0.6])

pl.subplot(2,2,2)
delta = detuw + k * v - zs.uprimehbar * field(z)
pl.plot(z, delta/(atom.G/2))
pl.plot(z, dres)
pl.xlim([0, 0.6])
pl.ylim([-5, 1])

pl.subplot(2,2,3)
dz = np.diff(delta)/np.diff(z)
Bz = np.diff(field(z))/np.diff(z)
vz = np.diff(v)/np.diff(z)

pl.plot(z[1:], dz, label='dz')
pl.plot(z[1:], -zs.uprimehbar*Bz, label='Bz')
pl.plot(z[1:], k*vz, label='vz')
pl.plot(z[1:], k*vz-zs.uprimehbar*Bz)
# pl.plot(z[1:], zerod(res[1:], field(z[1:]), Bz), 'k--')
# pl.plot(z, vals, 'b--')
pl.legend(loc='best')

pl.xlim([0, 0.6])

pl.subplot(2,2,4)
# dv = k*np.diff(delta)/np.diff(vz)
# dB = -zs.uprimehbar * np.diff(delta)/np.diff(field(z))
dv = np.diff(delta)/np.diff(v)
dB = np.diff(delta)/np.diff(field(z))
vz = np.diff(v)/np.diff(z)
Bz = np.diff(field(z))/np.diff(z)
pl.plot(z[1:], dv*vz)
pl.plot(z[1:], dB*Bz)
pl.xlim([0, 0.6])
# lim = 1
# pl.ylim([-lim, lim])

pl.show()
