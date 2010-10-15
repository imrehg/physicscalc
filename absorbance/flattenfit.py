#!/usr/bin/env python
#################
# Fitting absorbance data, all together now
from __future__ import division
from numpy import *
import numpy as np
import pylab as pl
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.integrate import odeint
import sys

### Load data
filename = "grab_results.txt"
grab = loadtxt(filename, 
              dtype={'names': ('series', 'power', 'sens'),
                     'formats': ('i','f','S3')},
                 comments=";", delimiter=",")

sy = [170, 240]
sx = [170, 440]
y = array(range(sy[0], sy[1]))
x = array(range(sx[0], sx[1]))
X, Y = pl.meshgrid(x, y)

def loadseries(ser):
    pics = []
    i = 0
    while True:
        try:
            pics += [pl.imread('grab_%04d_%02d.png' %(ser, i))]
            i += 1
        except:
            break
    a = pics[0]
    for k in range(1, i):
        a += pics[k]
    a /= i
    a = a[sy[0]:sy[1], sx[0]:sx[1]]
    print "Series %04d : ...  -> %f" %(ser, np.max(a))
    return a

pics = [loadseries(ser) for ser in grab['series']]
ipics = [sum(pic, axis=0) for pic in pics]
gains = [2**int(sens, 2) for sens in grab['sens']]
#############


# Generate interpolation for z->I and I->z mapping
def beamprop(I, z, a):
    """ Differential equation for the beam propagation with saturation """
    return [-a/sqrt(1+I[0])*I[0]]
a0 = 1
I0 = 1000
z = linspace(0, 120, 501)
Iz = odeint(beamprop, I0, z, args=(a0,)).T[0]
zfromi = interp1d(-Iz, z, bounds_error=False, fill_value=0)
ifromz = interp1d(z, Iz, bounds_error=False, fill_value=0)

ssq = lambda x : sqrt(sum(x**2))

def fitx(p, x, zfromi, ifromz):
    """ Fitting function """
    I0, x0, a0, bg, A = p
    zI0 = zfromi(-I0)
    x0 = 150
    out = A*ifromz((x-x0)*a0+zI0) + bg
    return out

def errf(p, x, data, zfromi, ifromz):
    """ Fit a single picture """
    return fitx(p, x, zfromi, ifromz) - data

def errfmany(p, x, ipics, power, zfromi, ifromz):
    """ Fit multiple pictures together """
    x0, a0 = p[0:2]
    out = array([])
    for i, ipic in enumerate(ipics):
        I0 = p[2+i]
        A = p[2+2*len(ipics)+i]
        # bg = p[2+len(ipics)+i]
        # if bg < 0:
        #     bg = 0
        # fixed background
        bg = 0.45
        # pin = [power[i]/I0, x0, a0, bg, A]
        pin = [I0, x0, a0, bg, A]
        funcdiff = fitx(pin, x, zfromi, ifromz)-ipic
        out = append(out, funcdiff)
    return out

##### Starting fitting parameters
I0 = 3
x0 = 75
a0 = 1/80
bg = 0.5
A = 0.5
p = [x0, a0]
p += [I0]*len(ipics)
p += [bg]*len(ipics)
p += [A]*len(ipics)

### Single fit
# ni = 0
# pl.plot(x, ipics[ni], '.')
# pout, success = leastsq(errf, p, args=(x, ipics[ni], zfromi, ifromz))
# print pout
# pl.plot(x, fitx(pout, x, zfromi, ifromz))


pout, success = leastsq(errfmany, p, args=(x, ipics, grab['power'], zfromi, ifromz))
print success
# pout = p
x0, a0 = pout[0:2]
print "Rss             : %f" %(ssq(errfmany(pout, x, ipics, grab['power'], zfromi, ifromz)))
print "Absorbance 1/a  : %f" %(1/a0)
print "Beam start x0   : %f" %(x0)
print "Powers:"
for I0 in pout[2:2+len(ipics)]:
    print "%f" %I0
# print "Background:"
# for bg in pout[2+len(ipics):2+2*len(ipics)]:
#     print "%f" %bg
print "Gain:"
for A in pout[2+2*len(ipics):2+3*len(ipics)]:
    print "%f" %A
filenames = []
for i, ipic in enumerate(ipics):
    # bg = pout[4+i]
    I0 = pout[2+i]
    bg = pout[2+len(ipics)+i]
    A = pout[2+2*len(ipics)+i]
    bg = 0.45
    # pin = [grab['power'][i]/I0, x0, a0, bg, A]
    pin = [I0, x0, a0, bg, A]
    yhat = fitx(pin, x, zfromi, ifromz)
    pl.figure()
    pl.plot(x, ipic, 'kx')
    pl.plot(x, yhat, 'r-')
    pl.title("Series %04d: %.1f uW / %.1f Isat; 1/abs: %.1f px" %(grab['series'][i], grab['power'][i], I0, 1/a0))
    pl.xlim([x[0], x[-1]])
    pl.xlabel("x, Position (pixel)")
    pl.ylabel("Signal (integrated in y)")
    filename = "grab_%04d_fit.pdf" %grab['series'][i]
    filenames += [filename]
    pl.savefig(filename)

### Save one merged pdf
# import os
# mergepdf = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sPAPERSIZE=a4 -sOutputFile=grab_fitted.pdf"
# for fn in filenames:
#     mergepdf += " %s" %fn
# os.system(mergepdf)
pl.show()
