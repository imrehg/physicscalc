from __future__ import division
import numpy as np
import pylab as pl
from scipy import optimize

def gauss2d(cx, cy, sx, sy, h, X, Y):
    res = h * np.exp(-(X - cx)**2 / (2*sx**2)) * np.exp(-(Y - cy)**2 / (2*sy**2))
    return res

def waves(th, period, phase, amp, X, Y):
    pos = np.cos(th)*X + np.sin(th)*Y
    ret = np.sin(2*np.pi*pos/period + phase) * amp
    return ret

def wobble(cx, cy, sx, sy, h, th, period, phase, amp, offset):
    return lambda X, Y: gauss2d(cx, cy, sx, sy, h, X, Y) * (1 + waves(th, period, phase, amp, X, Y)) + offset

def fitwobble(data, X, Y, p0):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    errorfunction = lambda p: np.ravel(wobble(*p)(X, Y) - data)
    p, success = optimize.leastsq(errorfunction, p0)
    return p

if __name__ == "__main__":
    filename = "img_418733.npy"
    data = np.load(filename)
    dimy, dimx = np.shape(data)
    X, Y = np.mgrid[0:dimx,0:dimy]
    X = X.T
    Y = Y.T

    th = -56/180*np.pi
    period = 36.5
    phase = 0.2
    amp = 0.2
    offset = data[0, 0]
    cx = 314
    cy = 260
    sx = 100
    sy = 100
    h = 175

    p0 = (cx, cy, sx, sy, h, th, period, phase, amp, offset)
    # # fitfun = np.ravel(wobble(*p)(X, Y) - data)
    # # print fitfun
    p = fitwobble(data, X, Y, p0)
    # print p
  #   p = [ 3.22891320e+02,   2.59297580e+02,   1.40109060e+02,   1.40215415e+02,
  #  1.42835103e+02,  -1.05651955e+00,   3.65765445e+01,   1.36331418e+00,
  # -1.22638586e-01,   8.96411783e-01,]
    # p = p0
    vals = wobble(*p)(X, Y)
    cx, cy, sx, sy, h, th, period, phase, amp, offset, = p
    print th, period, phase, amp

    fh, fw = 8.27, 11.6
    # Plotting
    pl.figure(num=1, figsize=(fw, fh))
    # Original data
    pl.subplot(221)
    pl.imshow(data)
    pl.title('Original data')

    pl.subplot(223)
    pl.imshow(vals)
    pl.title('Fitted function')

    # Residuals
    pl.subplot(224)
    pl.imshow(data - vals)
    pl.title('residuals')
    pl.xlabel('x')
    pl.xlabel('y')

    thdeg = (th/np.pi*180)
    print "Angle: %g deg " %thdeg
    pl.subplot(222)
    pl.imshow(data)
    pl.title('angle %.1f deg, period: %.1f pixel, sx/y: (%.1f/%.1f pixel)' %(thdeg, period, sx, sy))
    centre, = pl.plot(cx, cy, '+', markersize=20)
    thlist = np.linspace(0, 2*np.pi, 101)
    xe = 2*sx*np.cos(thlist) + cx
    ye = 2*sy*np.sin(thlist) + cy
    pl.plot(xe, ye, 'k-')

    ax = period*np.cos(th) + cx
    ay = period*np.sin(th) + cy
    pl.plot([cx, ax], [cy, ay], 'k-', linewidth=2)

    centre.axes.set_xlim([0, dimx])
    centre.axes.set_ylim([dimy-1, 0])

    print np.sqrt(np.sum((data - vals)**2))

    # pl.figure()
    # pl.subplot(211)
    # xstep = np.array(range(0, 640))
    # cuty = cy-200
    # xcut = data[cuty, :]
    # # x2 = np.exp( - (xstep - cx)**2 / (2 * sx**2))*h + offset
    # x3 = np.exp( - (xstep - cx)**2 / (2 * sx**2)) * h * (1 + 0.2*np.sin(2*np.pi*xstep*np.cos(th)/36 + 2.5))
    # ww = wobble(cx, cy, sx, sy, h, th, period, phase, amp, offset)
    # x2 = ww(xstep, cuty)

    # pl.plot(xstep, xcut, 'k-')
    # pl.plot(xstep, x2, 'r-')
    # pl.title(cy)

    # pl.subplot(212)
    # pl.plot(xstep, (xcut-x2) / xcut)

    pl.show()
