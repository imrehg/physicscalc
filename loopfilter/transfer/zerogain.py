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
    return deg

def angles(v):
    """ Make a numpy array into angles """
    return np.array([angle(vi) for vi in v])

def tfsplit(v):
    return 20*np.log10(np.abs(v)), angles(v)

def tfsplitdb(v):
    return np.abs(v), angles(v)

### Overall settings
flim0 = 10
flim1 = 5000
###

# Phase detector
Kf = 5e-3  # 5mA / 2pi rad

# Loop filter
def loopf(p):
    C1, C2, R2 = p
    a = 1 / (C2 * R2)
    b = (C1 + C2) / (C1 * C2 * R2)
    c = 1/C1
    return lambda s: c / s * (s + a) / (s + b)

p0 = (10e-6, 47e-6, 68)
loop = loopf(p0)

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
pidinterp = lambda f : 10**(ainterp(f)/20)*(np.cos(pinterp(f)) + 1j*np.sin(pinterp(f)))

# x1, x2 = tfsplit(pidinterp(fp))
# # x1 = -np.abs(pidinterp(fp))
# # pl.semilogx(fp, x1, 'x')
# pl.semilogx(fp, x2, 'x')

f = np.logspace(np.log10(flim0), np.log10(flim1), 101)
s = 2*np.pi*1j*f


# t = Kf * loop(s) * pidinterp(f)
t = Kf * loop(s)
amp, ph = tfsplit(t)

pl.figure(figsize=(8.27, 11.69))
pl.subplot(211)
pl.semilogx(f, amp)
pl.xlabel("frequency (Hz)")
pl.ylabel("Amplitude")

pl.subplot(212)
pl.semilogx(f, ph)
# pl.semilogx(f, ppdi(f))
pl.semilogx(f, ph + ppdi(f))
pl.xlabel("frequency (Hz)")
pl.ylabel("Phase")



pl.show()
