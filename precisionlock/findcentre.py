#!/usr/bin/python
import numpy as np
import pylab as pl

lorentz = lambda p, x: p[2]*p[1]/((x-p[0])**2 + p[1]**2)



x = np.linspace(-5, 5, 501)
p = [0, 1, 1]

y = lorentz(p, x) + (x+1)**2 * 0.05

pl.subplot(211)
pl.plot(x, y)
pl.subplot(212)
dy = np.diff(y)
dx = x[0:-1]+(x[1]-x[0])

pl.plot(dx, dy)
pl.show()
