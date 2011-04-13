#!/usr/bin/python
from __future__ import division
import numpy as np
import pylab as pl

def fbres(R, d):
    T = 1 - R
    

    maxn = 30
    # Iout2 = np.abs(np.sum([np.exp(1j*i*d)*np.sqrt(T)*R**i for i in xrange(maxn)]))
    E = sum(np.array([np.exp(1j*i*d)*np.sqrt(T)*R**i for i in xrange(maxn)]))
    Iout2 = E * np.conj(E)
    Iout = (1-R)**2 / ((1-R)**2 + 4*R*np.sin(d/2)**2)
    return (Iout, Iout2*(1-R))


# Iout = fbres(0.96, 0.5)
# print Iout
R = 0.5
d = np.linspace(0, 4*np.pi, 101)

pl.plot(d, [fbres(R, di)[0] for di in d], 'b-')
pl.plot(d, [fbres(R, di)[1] for di in d], 'r.')
pl.show()
