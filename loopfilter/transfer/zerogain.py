#!/usr/bin/env python2
from __future__ import division
import numpy as np
import pylab as pl
import scipy.optimize as op
import scipy.interpolate as ip

def angle(v):
    """ Complex number to angle conversion """
    try:
        deg = np.arctan(v.imag / v.real) / np.pi * 180
    except ZeroDivisionError:
        if v.imag >= 0:
            deg = 90
        else:
            deg = 270
    deg = (deg - 180) if deg > 0 else deg
    return deg

def angles(v):
    """ Make a numpy array into angles """
    return np.array([angle(vi) for vi in v])

def tfsplit(v):
    if type(v) != type(np.array([1])):
        v = np.array([v])
    return np.log10(np.abs(v)), angles(v)

def tfsplitdb(v):
    return np.abs(v), angles(v)

### Overall settings
flim0 = 10
flim1 = 5000
###

# Phase detector
Kf = 5e-3 / 2*np.pi  # 5mA / 2pi rad

# Loop filter
def loopf(p):
    C1, C2, R2 = p
    a = 1 / (C2 * R2)
    b = (C1 + C2) / (C1 * C2 * R2)
    c = 1/C1
    return lambda s: c / s * (s + a) / (s + b)

p0 = (10e-6, 47e-6, 68)
loop = loopf(p0)

# def loop2(p):
#     C1, C2, R2 = p
#     T1 = R2 * C1*C2 / (C1 + C2)
#     T2 = R2 * C2
    

# # Experimental loop filter
# dloop = np.loadtxt("transfer_120413_151455.log", delimiter=",", comments="#")
# fl = dloop[:, 0]
# al = dloop[:, 1]
# phl = dloop[:, 3]
# al = al/al[0]
# alinterp = ip.interp1d(fl, al)
# plinterp = ip.interp1d(fl, phl)
# pllinterp = lambda f : alinterp(f)*(np.cos(plinterp(f)) + 1j*np.sin(plinterp(f)))


# PID loop (from experiment
dpid = np.loadtxt("sweep_120409_103951.log", delimiter=",", comments="#")
fp = dpid[:, 0]
ap = dpid[:, 1]
ppd = dpid[:, 2]
# pl.semilogx(fp, ppd)
ppdi = ip.interp1d(fp, ppd)
pp = ppd/180*np.pi
ainterp = ip.interp1d(fp, ap)
pinterp = ip.interp1d(fp, pp)
ap = ap[fp <= flim1]
pp = pp[fp <= flim1]
fp = fp[fp <= flim1]
ap = ap[fp >= flim0]
pp = pp[fp >= flim0]
fp = fp[fp >= flim0]
# ainterp = ip.interp1d(fp, ap)
# pinterp = ip.interp1d(fp, pp)
# pidinterp = lambda f : 10**(ainterp(f)/20)*(np.cos(pinterp(f)) + 1j*np.sin(pinterp(f)))
pidinterp = lambda f : 10**(ainterp(f)/10)*(np.cos(pinterp(f)) + 1j*np.sin(pinterp(f)))

# PZT from experiment
dpzt = np.loadtxt("pzt_120409_160557_2_fixed.log", delimiter=",", comments="#")
fz = dpzt[:, 0]
az = dpzt[:, 1]
pz = dpzt[:, 3]
az = az/az[0]
azinterp = ip.interp1d(fz, az)
pzinterp = ip.interp1d(fz, pz)
pzzinterp = lambda f : azinterp(f)*(np.cos(pzinterp(f)) + 1j*np.sin(pzinterp(f)))

Klaser = 1.6e7  # Hz/Volt for the laser



# x1, x2 = tfsplit(pidinterp(fp))
# # x1 = -np.abs(pidinterp(fp))
# # pl.semilogx(fp, x1, 'x')
# pl.semilogx(fp, x2, 'x')

f = np.logspace(np.log10(flim0), np.log10(flim1), 101)
imags = lambda f : 2*np.pi*1j*f
s = imags(f)


# t = Kf * loop(s) * pidinterp(f) * pzzinterp(f) * Klaser / s
t = Kf * loop(s) * pidinterp(f) * 15 *  Klaser / s * pzzinterp(f) / 32

ftest = 101
print "Testing at %g Hz" %(ftest)
tl = loop(imags(ftest))
print "Loop", tfsplit(tl)
tl = pidinterp(ftest)
print "PID", tfsplit(tl)
tl = pzzinterp(ftest)
print "PZT", tfsplit(tl)
tl = Kf * loop(imags(ftest)) * pidinterp(ftest) * 15 *  Klaser / imags(ftest) * pzzinterp(ftest) / 32
print "total", tfsplit(tl)



# t1 = Kf * loop(s) * pidinterp(f)

# t1 = pllinterp(f)

# t = pzzinterp(f)
# t = pidinterp(f)
# t = Kf * loop(s)
amp, ph = tfsplit(t)
# amp1, ph1 = tfsplit(t1)

pl.figure(figsize=(8.27, 11.69))
pl.figtext(0.5, 0.95,  '', ha='center', size='x-large')

pl.subplot(211)
pl.semilogx(f, amp, 'k-')
pl.semilogx([f[0], f[-1]], [0, 0], 'k--')
# pl.semilogx(f, amp1)
# pl.semilogx(f, amp+ainterp(f))
pl.xlabel("frequency (Hz)")
pl.ylabel("LogAmplitude")

pl.subplot(212)
pl.semilogx(f, ph)
# pl.semilogx(f, ph1)
# pl.semilogx(f, ppdi(f))
# pl.semilogx(f, ph + ppdi(f))
pl.xlabel("frequency (Hz)")
pl.ylabel("Phase")



# # t = Kf * loop(s)
# # amp, ph = tfsplit(t)
# pl.subplot(211)
# pl.semilogx(f, amp)
# pl.xlabel("frequency (Hz)")
# pl.ylabel("Amplitude")

# pl.subplot(212)
# pl.semilogx(f, ph)
# # pl.semilogx(f, ppdi(f))
# pl.semilogx(f, ph + ppdi(f))
# pl.xlabel("frequency (Hz)")
# pl.ylabel("Phase")


pl.show()
