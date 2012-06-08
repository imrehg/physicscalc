from __future__ import division
import numpy as np
import pylab as pl
import csv

import layeroptimize as lo
import zeemanslower as zs

simfile = "19_AWG12Coat_10.npz"
sim = np.load(simfile)['simulation'][()]
atom = zs.Rb87()
wire = sim['wire']
setup = sim['setup']
eta = sim['eta']
v0 = sim['v0']
vf = sim['vf']
sl = zs.slowerlength(atom.aslow, eta, v0, vf)

looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']
eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']

bideal = zs.bideal(atom, 0, eta, v0, vf, detu)

z = np.append([-5.5*R], np.append(np.linspace(0, sl, 61), [sl+5.5*R]))
ze = np.linspace(z[0], z[-1], 201)
ze = np.linspace(0.4, z[-1], 201)

# for i, s in enumerate(segments):
#     print i, s, layer[i]

f1 = lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup)
pl.plot(ze, lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup), 'k-', linewidth=1)

print setup['segments']
segments[18] = [302, 302]
segments[19] = [302, 304]
segments[20] = [304, 304]

print setup['segments']

f2 = lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup)
pl.plot(ze, lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup), 'r-', linewidth=1)


simulation = {'R': R,
              'eta': eta,
              'v0': v0,
              'vf': vf,
              'detu': detu,
              'wire': wire,
              'setup': setup
              }
series = 19
nlayer = 10
savefile = "NEW%d_%s_%d" %(series, wire[2], nlayer)
np.savez('%s' %(savefile), simulation=simulation)


# pl.plot(ze, (f2-f1)/abs(f1))

pl.show()
