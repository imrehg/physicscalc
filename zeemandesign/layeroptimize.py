#!/usr/bin/env python2.7
"""
Simulate multi-layer and optimize layer boundaries
"""
import numpy as np
import numpy.random as random
import pylab as pl
import scipy.odr as odr
from time import strftime

mu0 = 1.257e-6 # Tm/A
R0 = 0.0383
# R0 = 0.025
d0 = 0.008

class memoized(object):
   """Decorator that caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned, and
   not re-evaluated.

   TODO: add keyword arguments:
   http://pko.ch/2008/08/22/memoization-in-python-easier-than-what-it-should-be/
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      try:
         return self.cache[args]
      except KeyError:
         value = self.func(*args)
         self.cache[args] = value
         return value
      except TypeError:
         # uncachable -- for instance, passing a list as an argument.
         # Better to not cache than to blow up entirely.
         return self.func(*args)
   def __repr__(self):
      """Return the function's docstring."""
      return self.func.__doc__
   def __get__(self, obj, objtype):
      """Support instance methods."""
      return functools.partial(self.__call__, obj)

@memoized
def loopfield(pos, n=1, R=R0, d=d0):
    """ Magnetic field by a multilayer loop """
    x = np.array(pos)
    Ri = R
    if n < 1:
        return x*0
    field = mu0 * Ri*Ri / (2 * (Ri*Ri+x*x) ** (1.5))
    for i in xrange(1, n):
        Ri += d
        field += mu0 * Ri*Ri / (2 * (Ri*Ri+x*x) ** (1.5))
    return field

def coillength(n, R, d):
    """ Approximate length of wires in a single multiple layer setting """
    Ri = R
    if n < 0:
        return 0
    length = 2*np.pi* Ri
    for i in xrange(1, n):
        Ri += d
        length += 2*np.pi* Ri
    return length

def totalcoillength(setup, R=R0, d=d0):
    layernum, csign = getLayerNumber(setup)
    cl = sum([coillength(layer, R, d) for layer in layernum])
    return cl

def fieldcalc(z, setup, curr=1, R=R0, d=d0):
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
        out += multi*loopfield(tuple(z-czi), int(di), R, d)
    return out

def bideal(z):
    C1 = 0.015
    C2 = 0.03
    C3 = 1
    C4 = 1.13
    return -(C1 - C2 * np.sqrt(C3**2  - C4 * z))    

def bideal2(z):
    extra = 5
    zextra = np.append(np.append(-extra*R0+z[0], z), z[-1]+extra*R0)
    field = normalize(np.append(np.append(0, bideal(z)), 0), 1)
    return (zextra, field)

def getLayerNumber(setup):
    """ Calculate the number of layers at the different loop positions """
    out = []
    outlayer = []
    for si, seglim in enumerate(setup['segments']):
        for i in range(seglim[0], seglim[1]):
            out += [setup['layer'][si]]
            outlayer += [setup['csign'][si]]
    return out, outlayer

def normalize(values, i=0):
    """ Normalize by values """
    return values / values[i]

def ssq(val):
    """ squareroot of sum of squares """
    return np.sqrt(np.sum(val**2))

def seglen(segment):
    return segment[1] - segment[0]

def optimize(zl, setup, maxtry=100):
    zideal, fideal = bideal2(zl)
    foriginal = normalize(fieldcalc(zideal, setup), 1)

    oldval = ssq(fideal - foriginal)
    
    nchoice = len(setup['layer']) - 1
    n = 0
    Temp = 0.1
    while n < maxtry:
        if (n+1) % 100 == 0:
            Temp /= 1.05
            print n+1, oldval
        segment = random.randint(0, nchoice)
        lowseg, highseg = setup['segments'][segment], setup['segments'][segment+1]
        decrease = random.randint(0, 2)
        if (((decrease == 0) and seglen(lowseg) < 1) or
            ((decrease == 1) and seglen(highseg) < 1)):
            continue

        # Change segments
        if decrease == 0:
            lowseg[1] -= 1
            highseg[0] -= 1
        else:
            lowseg[1] += 1
            highseg[0] += 1

        field = normalize(fieldcalc(zideal, setup), 1)
        newval = ssq(fideal - field)
        # print np.exp((-newval + oldval)/Temp)
        if np.exp((-newval + oldval)/Temp) > random.rand():
            # print "Kept: ", newval, oldval
            oldval = newval
        else:
            # print "Reversed", newval, oldval
            if decrease == 0:
                lowseg[1] += 1
                highseg[0] += 1
            else:
                lowseg[1] -= 1
                highseg[0] -= 1
        n += 1
    return setup

if __name__ == "__main__":


    # The field that we want to match
    nump = 61
    zl = np.linspace(0, 0.88495, nump)

    ze, fe = bideal2(zl)

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
    
    # # # How our cross cut looks 
    # # pl.plot(setup['looppos'], getLayerNumber(setup))
    # # pl.show()

    

    # # ideal = bideal(zl)
    # # ideal = ideal / ideal[0]
    # # calcfield = fieldcalc(zl, setup)
    # # calcfield = calcfield / calcfield[0]
    # # pl.plot(zl, ideal, 'x')
    # # pl.plot(zl, calcfield)
    # # pl.plot([0, 0.9], [0, 0])
    # # pl.show()

    # # # import time
    # # # start = time.time()
    # # # for i in xrange(100):
    # # #     fieldcalc(zl, setup)
    # # # print time.time()-start

    # pl.figure()
    # pl.plot(ze, fe, 'x')
    # pl.plot(ze, normalize(fieldcalc(ze, setup), 1), 'g--')
    newsetup = optimize(zl, setup, maxtry=25000)
    zz = np.append(np.linspace(-6*R0, 0.1, 10),  np.linspace(0, 1.1, 101))
    pl.plot(ze, fe, 'x')
    pl.plot(ze, normalize(fieldcalc(ze, setup), 1), 'g--')
    pl.plot(zz, normalize(fieldcalc(zz, newsetup), 10), 'r-')
    # print newsetup['segments']

    # from tempfile import NamedTemporaryFile
    # outfile = NamedTemporaryFile(prefix='layer', suffix='.npz', dir='.', delete=False)
    outfile = "layers_%s.npz" %(strftime("%y%m%d_%H%M%S"))
    print outfile
    np.savez(outfile, setup=setup, ideal=bideal2(zl), d0=d0)

    pl.show()
