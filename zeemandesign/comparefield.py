import numpy as np
import pylab as pl

import zeemanslower as zs
import layeroptimize as lo


current = 2
outname = "theory1"
delta = 0
title = 'Zeeman slower, I=%.1fA' %(current)

current = 2
outname = "theory2"
delta = -0.01
title = 'Zeeman slower, I=%.1fA, Negative coil shift=%.3fm' %(current, delta)

measurement = "measuredata.csv"
mesdata = np.loadtxt(measurement, skiprows=1, delimiter=",")

simfile = "19_AWG12Coat_10.npz"
atom = zs.Rb87()
sim = np.load(simfile)['simulation'][()]
wire = sim['wire']
setup = sim['setup']
eta = sim['eta']
v0 = sim['v0']
vf = sim['vf']
looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']

# delta = -0.01
# delta = 0
for i in range(len(looppos)):
    if i > 239:
        looppos[i] += delta
cstart = looppos[segments[1][0]]

eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']
sl = zs.slowerlength(atom.aslow, eta, v0, vf)
z = np.array([0])
bfield = zs.bideal(atom, z, eta, v0, vf, detu)
z = np.append([-5.5*R], np.append(np.linspace(0, sl, 61), [sl+5.5*R]))
ze = np.linspace(z[0], z[-1], 201)

# delta2 = 0.008
delta2 = 0
zadjust = -0.019-0.058+cstart+delta2
zpos = mesdata[:, 0]/100 + zadjust
mesfield = mesdata[:, 1]


pl.figure(figsize=(11.69, 8.27))
pl.plot(ze, lo.fieldcalc(ze, setup)*current * 1e4, 'r-', linewidth=3, label='Simulation')
pl.plot(zpos, mesfield, 'ko', label='Measured')
pl.xlabel("Position (m)", fontsize=16)
pl.ylabel("Field (G)", fontsize=16)
pl.legend(loc='upper right')
pl.title(title)

pl.savefig("%s.pdf" %(outname))
pl.savefig("%s.png" %(outname))

pl.show()


