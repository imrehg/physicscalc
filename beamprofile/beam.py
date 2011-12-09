"""
Beam diameter calculation using the D4s, or second moment width
from http://en.wikipedia.org/wiki/Laser_beam_profiler#Measurements
"""
from __future__ import division
import numpy as np
import pylab as pl
import time

if __name__ == "__main__":
    X, Y = np.mgrid[0:150.,0:100.]
    xc =50
    yc = 50
    sx = 15
    sy = sx/2

    xg = np.exp(-(X - xc)**2 / (2 * sx**2))
    yg = np.exp(-(Y - yc)**2 / (2 * sy**2))
    gg = (xg * yg).T

    P = np.sum(gg)
    xx = np.sum(gg * X.T) / P
    yy = np.sum(gg * Y.T) / P
    xx2 = np.sum(gg * (X.T - xx)**2)/P
    yy2 = np.sum(gg * (Y.T - yy)**2)/P
    xy = np.sum(gg * (X.T - xx) * (Y.T - yy)) / P
    gamm = np.sign(xx2 - yy2)
    angle = 0.5 * np.arctan(2*xy / (xx2 - yy2))

    dx = 2 * np.sqrt(2) * (xx2 + yy2 + gamm * ( (xx2 - yy2)**2 + 4*xy**2)**0.5)**(0.5)
    dy = 2 * np.sqrt(2) * (xx2 + yy2 - gamm * ( (xx2 - yy2)**2 + 4*xy**2)**0.5)**(0.5)

    print dx, dy
    print angle

    pl.imshow(gg)
    a = np.linspace(0, np.pi*2, 101)
    xr = np.cos(a)*dx/2 + xx
    yr = np.sin(a)*dy/2 + yy
    pl.plot(xr, yr, 'k-')
    pl.xlim([0, 150])
    pl.ylim([0, 100])
    pl.show()
