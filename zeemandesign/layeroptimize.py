#!/usr/bin/env python2.7
"""
Simulate multi-layer and optimize layer boundaries
"""
import numpy as np
import numpy.random as random
import pylab as pl
import scipy.odr as odr

mu0 = 1.257e-6 # Tm/A
R0 = 0.0383
d0 = 0.004

def loopfield(x, n=1, R=R0, d=d0):
    """ Magnetic field by a multilayer loop """
    Ri = R
    if n < 1:
        return x*0
    field = mu0 * Ri*Ri / (2 * (Ri*Ri+x*x) ** (1.5))
    for i in xrange(1, n):
        Ri += d
        field += mu0 * Ri*Ri / (2 * (Ri*Ri+x*x) ** (1.5))
    return field

def fieldcalc(z, setup, curr=1):
    if type(z) == type(1):
        z = np.array([z])
    elif type(z) == type([1]):
        z = np.array(z)
    out = np.zeros(z.size)
    looppos = setup['looppos']
    layernum, csign = getLayerNumber(setup)
    for i, val in enumerate(zip(looppos, layernum)):
        czi, di = val
        multi = csign[i] * curr
        out += multi*loopfield(z-czi, int(di), R=R0, d=d0)
    return out

def bideal(z):
    C1 = 0.015
    C2 = 0.03
    C3 = 1
    C4 = 1.13
    return -(C1 - C2 * np.sqrt(C3**2  - C4 * z))    

def getLayerNumber(setup):
    """ Calculate the number of layers at the different loop positions """
    out = []
    outlayer = []
    for si, seglim in enumerate(setup['segments']):
        for i in range(seglim[0], seglim[1]):
            out += [setup['layer'][si]]
            outlayer += [setup['csign'][si]]
    return out, outlayer

if __name__ == "__main__":

    # The field that we want to match
    nump = 41
    zl = np.linspace(0, 0.88, nump)

    # Our coil
    coilpos = (-0.2, 1.1)
    nloops = int((coilpos[1] - coilpos[0]) / d0)
    loops = [coilpos[0]+i*d0 for i in xrange(nloops)]
    # The number of layers
    layer = [0,10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 0]
    # The direction of current
    csign = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1]

    defaultLayerLength = int(nloops / len(layer))
    start = 0
    segments = [[defaultLayerLength*i, defaultLayerLength*(i+1)] for i in range(len(layer))]
    segments[-1][-1] = nloops

    setup = {'looppos': loops,
             'layer': layer,
             'csign': csign,
             'segments': segments,
             }
    
    # # How our cross cut looks 
    # pl.plot(setup['loops'], getLayerNumber(setup))
    # pl.show()

    # pl.plot(zl, bideal(zl))
    # pl.show()

    # Looks like good field shape
    pl.plot(zl, fieldcalc(zl, setup))
    pl.show()
