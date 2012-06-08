import numpy as np
import pylab as pl

import zeemanslower as zs
import layeroptimize as lo


## Measured positions
## Format: Start, finish, layer, turns, current
coils = [(0, 0.435, 1, 207, 1),
         (0, 0.395, 2, 186, 1),
         (0, 0.345, 3, 163, 1),
         (0, 0.290, 4, 137, 1),
         (0, 0.239, 5, 110, 1),
         (0, 0.173, 6, 80, 1),
         (0, 0.108, 7, 49, 1),
         (0.005, 0.045, 8, 18, 1),
         (0.005, 0.045, 9, 18, 1),
         (0.473, 0.635, 1, 75, -1),
         (0.508, 0.635, 2, 58, -1),
         (0.537, 0.635, 3, 43, -1),
         (0.558, 0.635, 4, 31, -1),
         (0.595, 0.629, 5, 13, -1),
         (0.595, 0.629, 6, 13, -1),
         (0.595, 0.629, 7, 13, -1),
         (0.595, 0.629, 8, 13, -1),
         (0.600, 0.627, 9, 11, -1),
         (0.600, 0.627, 10, 11, -1),
         ]

wthick = 0.0021
R0 = 0.019/2
current = 2
title = 'Field comparison, current=%.1fA' %(current)
outname = 'slowerfield3'

# for c in coils:
#     l = (c[1]-c[0])
#     ul = l / c[3]
#     print l, ul, ul/wthick


zex = np.linspace(-0.03225, 0.66225, 501)
xfield = zex*0
for c in coils:
    R = R0 + (c[2]-0.5)*wthick
    pos = np.linspace(c[0], c[1], c[3])
    print c
    for p in pos:
        xfield += lo.loopfield(zex-p, 1, R, wthick)*c[4]*1e4

# current = 2
# outname = "theory1"
# delta = 0
# title = 'Zeeman slower, I=%.1fA' %(current)

# current = 2
# outname = "theory2"
# delta = -0.01
# title = 'Zeeman slower, I=%.1fA, Negative coil shift=%.3fm' %(current, delta)

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
delta = 0
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
print ze[0], ze[-1]

# delta2 = 0.008
delta2 = 0
zstart = -0.019
zadjust = -0.019-0.058+cstart+delta2
zpos = mesdata[:, 0]/100 + zadjust
mesfield = mesdata[:, 1]



pl.figure(figsize=(11.69, 8.27))
pl.plot(ze, lo.fieldcalc(ze, setup)*current * 1e4, 'r-', linewidth=3, label='Simulation with calculated positions')
pl.plot(zex+zstart, xfield*current, 'b-', linewidth=3, label='Simulation with measured positions')
pl.plot(zpos, mesfield, 'ko', label='Measured')
pl.xlabel("Position (m)", fontsize=16)
pl.ylabel("Field (G)", fontsize=16)
pl.legend(loc='upper right')
pl.title(title)

pl.savefig("%s.pdf" %(outname))
pl.savefig("%s.png" %(outname))

pl.show()
