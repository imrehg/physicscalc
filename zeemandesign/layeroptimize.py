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
# R0 = 0.0383
# # R0 = 0.025
# d0 = 0.008

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
def loopfield(pos, n, R, d):
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

def totalcoillength(setup, R, d):
    layernum, csign = getLayerNumber(setup)
    cl = sum([coillength(layer, R, d) for layer in layernum])
    return cl

def fieldcalc(z, setup, curr=1):
    R = setup['R']
    d = setup['d']
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

# def bideal(z):
#     C1 = 0.015
#     C2 = 0.03
#     C3 = 1
#     C4 = 1.13
#     return -(C1 - C2 * np.sqrt(C3**2  - C4 * z))    

# def bideal2(z):
#     extra = 5
#     zextra = np.append(np.append(-extra*R0+z[0], z), z[-1]+extra*R0)
#     field = normalize(np.append(np.append(0, bideal(z)), 0), 1)
#     return (zextra, field)

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

def optimize(zideal, fideal, setup, maxtry=100, printprogress=True):
    # zideal, fideal = bideal2(zl)
    fideal = fideal/fideal[1]
    foriginal = fieldcalc(zideal, setup)/fieldcalc(0, setup)

    oldval = ssq(fideal - foriginal)
    # pl.figure()
    # pl.plot(zideal, fideal, 'o', label='target')
    # pl.plot(zideal, foriginal, 'x', label='current')
    # pl.legend(loc='best')

    nchoice = len(setup['layer']) - 1
    n = 0
    Temp = 0.1
    while n < maxtry:
        if (n+1) % 100 == 0:
            Temp /= 1.035
            if printprogress:
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

        field = fieldcalc(zideal, setup)/fieldcalc(0, setup)
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

def createStruct(z, d0, n):
    coilpos = (z[0], z[1])
    nloops = int((coilpos[1] - coilpos[0]) / d0)
    loops = [coilpos[0]+i*d0 for i in xrange(nloops)]

    np, nn = n
    layer = [0] + range(np, 0, -1) + range(0, nn+1) + [0]
    csign = [1]*(np+2) + [-1]*(nn+1)

    defaultLayerLength = int(nloops / len(layer))
    start = 0
    segments = [[defaultLayerLength*i, defaultLayerLength*(i+1)] for i in range(len(layer))]
    segments[-1][-1] = nloops

    return nloops, loops, layer, csign, segments

if __name__ == "__main__":
    import wires
    import zeemanslower as zs

    # Parameters
    R = 0.30 / 2 # larger diameter slower tube
    eta = 0.5 # efficiency
    Ls = 0.6 # set slower length
    vf = 30
    # The field that we want to match
    nz = 61
    atom = zs.Rb85()
    v0 = np.sqrt(2 * l2 * atom.aslow * eta + vf**2)
    print v0
    sl = zs.slowerlength(atom.aslow, eta, v0, vf)
    z = np.append([-5.5*R], np.append(np.linspace(0, sl, nz-2), sl+5.5*R))

    # bfield = zs.bideal(atom, z, eta, v0, vf, detu)


    # wire = wires.AWG9
    # # # Our coil
    # nloops, loops, layer, csign, segments = createStruct((z[0],z[-1]), wire[0], (10, 10))

    # setup = {'looppos': loops,
    #          'layer': layer,
    #          'csign': csign,
    #          'segments': segments,
    #          'R': R,
    #          'd': wire[0],
    #          }
    
    # newsetup = optimize(z, bfield, setup, maxtry=10000)
    # # newsetup = setup
    # # pl.plot(z, bfield, 'x')
    # ze = np.linspace(z[0], z[-1], 201)
    # pl.figure(figsize=(11.69, 8.27))
    # pl.plot(ze, fieldcalc(ze, newsetup)/fieldcalc(0, newsetup), 'r-', linewidth=2, label='coil field')
    # pl.plot(z, bfield/bfield[1], 'ko', markersize=5, label='target field')
    # pl.title("%s, R: %g mm, v: %d-%d m/s, %d MHz" %(wire[2], R*1e3, v0, vf, detu))
    # pl.xlabel('position (m)')
    # pl.ylabel('normalized magnetic field')
    # pl.legend(loc='best')
    # # pl.plot(zz, normalize(fieldcalc(zz, newsetup), 10), 'r-')
    # # print newsetup['segments']

    # simulation = {'R': R,
    #               'eta': eta,
    #               'v0': v0,
    #               'vf': vf,
    #               'detu': detu,
    #               'wire': wire,
    #               'setup': newsetup
    #               }
    # series = 0
    # savefile = "%d_%s" %(series, wire[2])
    # np.savez("%s", simulation=simulation)
    # pl.savefig('%s.png' %savefile)
    # pl.show()
    
