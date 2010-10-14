#!/usr/bin/env python
#################
# Fitting absorbance data
# usage: python analyze_avg.py <series number>
# where <series number> is e.g. 0004 in the filename grab_0004_00.png 
from __future__ import division
from numpy import *
import numpy as np
import pylab as pl
from scipy.optimize import leastsq
import sys

# a = pl.imread('grab_0000_0.png')
# a = pl.imread('0004.png')

### Which series?
if len(sys.argv) >= 2:
    ser = int(sys.argv[1])
    print "Series: %04d" %ser
else:
    ser = int(raw_input("Series number:"))
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
# Cropping data
sy = [175, 240]
sx = [170, 440]
y = range(sy[0], sy[1])
x = range(sx[0], sx[1])
X, Y = pl.meshgrid(x, y)
a = a[sy[0]:sy[1], sx[0]:sx[1]]
### Output: averages; check range values
amin = np.min(a)
amax = np.max(a)
print "Min/Max: ", amin, amax
saturate = ""
if amax > 0.96:
    saturate = ">>> SATURATION?"
    print saturate

### Print info file if exists
try:
    info = open('grab_%04d_info.txt' %ser, 'r')
    print "## Info:"
    for line in info.readlines():
        print line
    info.close()
    print "##"
except:
    pass

def beam(p, x, y):
    y0, th, A, xt, sy, off, offc = p
    # Coordinate transformation
    y = y - y0
    x = cos(th)*x + sin(th)*y
    y = -sin(th)*x + cos(th)*y
    return A/(sqrt(2*pi)*sy)*exp(-x/xt)*exp(-y**2/2/sy**2)+off*x + offc

gauss = lambda p, x: p[2]/(sqrt(2*pi)*p[1])*exp(-(x-p[0])**2/2/p[1]**2)+p[3]
           
pin = [200, 0.01, 1, 100, 3, -4e-5, 2e-2]
errfunc = lambda p, X, Y, a: ravel(beam(p, X, Y) - a)
p, success = leastsq(errfunc, pin, args=(X, Y, a))
fit = beam(p, X, Y)
print "RSS       : %f" %sqrt(sum(ravel(fit - a)**2))
print "Angle     : %f rad" %p[1]
print "Absorbance: %f" %p[3]
print "Width     : %f" %p[4]
print "Amplitude : %f" %p[2]
print "BG        : slope = %f ; offset = %f" %(p[5], p[6])

# Plotting results
pl.subplot(2,2,1)
pl.imshow(a)
pl.gray()
pl.xlabel('x (pixel)')
pl.ylabel('y (pixel)')
pl.title('Original %s' %saturate)

pl.subplot(2,2,2)
pl.imshow(fit)
pl.gray()
pl.xlabel('x (pixel)')
pl.ylabel('y (pixel)')
pl.title('Fit: xt=%.2f, sy=%.2f' %(p[3], p[4]))

pl.subplot(2,2,3)
pl.imshow(fit-a)
pl.gray()
pl.xlabel('x (pixel)')
pl.ylabel('y (pixel)')
pl.title('Residuals')


pl.subplot(2,2,4)
n = 0
yhat = linspace(y[0], y[-1], 401)
pg = [p[0], p[4], p[2]*exp(-x[n]/p[3]), p[5]*x[n]+p[6]]
pl.plot(y,a[:, n], '.')
pl.plot(yhat,gauss(pg, yhat))
pl.xlabel('y (pixel)')
pl.ylabel('Signal')
pl.title('Cross-section fit at x=%d' %n)
pl.xlim([y[0], y[-1]])
pl.show()
