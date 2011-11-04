"""
Repeat the calculation of the Cs slower from CLH's thesis
Short, low power consumption
"""
import numpy as np
import pylab as pl

import wires
import zeemanslower as zs
from layeroptimize import *


filename = "0_squaremagnet.npz"
atom = zs.Cs133()

sim = np.load(filename)['simulation'][()]
eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']
print eta, v0, vf, detu

sl = zs.slowerlength(atom.aslow, eta, v0, vf)
z = np.array([0])
bfield = zs.bideal(atom, z, eta, v0, vf, detu)

def getdata(filename):
    # Change from 0-d array into dict
    sim = np.load(filename)['simulation'][()]
    setup = sim['setup']
    coilfield = fieldcalc(z, setup, curr=1)
    current = bfield / coilfield
    totallen = totalcoillength(sim['setup'], sim['R'], sim['wire'][0])
    resistance = sim['wire'][1] * totallen
    power = resistance * current * current

    layernum, layerlayer = getLayerNumber(setup)
    nlayer = len(layernum)
    optsides = 0
    truesides = 0
    for i in xrange(nlayer):
        if  layernum[i] > 0:
            optsides += layernum[i]*4
            # The top is always open
            truesides += 1
            if (i > 0) and (layernum[i-1] < layernum[i]):
                truesides += layernum[i]-layernum[i-1]
            if (i < (nlayer-1)) and (layernum[i+1] < layernum[i]):
                truesides += layernum[i]-layernum[i+1]
                
    reducedarea = 1.0*truesides/optsides
    return current, totallen, power, resistance, reducedarea

newsetup = sim['setup']
layers = getLayerNumber(sim['setup'])
print layers

current, totallen, power, resistance, reducedarea = getdata(filename)
print "Current: %g A" %(current)
print "Totallen: %g m" %(totallen)
print "Power: %g W" %(power)
print "Resistance: %g w" %(resistance) 

z = np.append([-5.5*R], np.append(np.linspace(0, sl, 61), [sl+5.5*R]))
ze = np.linspace(z[0], z[-1], 201)

bf = zs.bideal(atom, z, eta, v0, vf, detu)
wire = sim['wire']

fig = pl.figure(figsize=(11.69, 8.27))
ax1 = fig.add_subplot(111)
# pl.plot(ze, fieldcalc(ze, newsetup)/fieldcalc(0, newsetup), 'r-', linewidth=2, label='coil field')
# pl.plot(z, bf/bfield[0], 'ko', markersize=5, label='target field')
p1 = ax1.plot(ze, fieldcalc(ze, newsetup)*current*1e4, 'r-', linewidth=2, label='coil field')
p2 = ax1.plot(z, bf*1e4, 'ko', markersize=5, label='target field')
ax2 = ax1.twinx()
p3 = ax2.plot(newsetup['looppos'], layers[0], linewidth=2, label='Layers (right scale)')
pl.title("%s, R: %g mm, v: %d-%d m/s, %d MHz, current: %g A" %(wire[2], R*1e3, v0, vf, detu, current))
ax1.set_xlabel('position (m)', fontsize=14)
ax1.set_ylabel('Magnetic field (G)', fontsize=14)
ax2.set_ylabel('Layer number', fontsize=14)
ax2.set_ylim([0, 10])
# leg = pl.legend(loc='best')
# leg.text[0].set_color(p1.get_color())
# # pl.legend(loc='best')
ax2.legend(loc='best')

pl.show()

# print getdata("102.npz")

# colours = ['b', 'g', 'r', 'y', 'k', 'm']
# for k, maxwind in enumerate(range(5, 11)):
#     powers = []
#     currents = []
#     lengths = []
#     resistances = []
#     reducedareas = []
#     csize = range(1, 11)
#     for i in csize:
#         # maxwind = 5
#         filename = "%d_AWG%d.npz" %(maxwind, i)
#         current, totallen, power, resistance, reducedarea = getdata(filename)
#         powers += [power]
#         currents += [current]
#         lengths += [totallen]
#         resistances += [resistance*1e3]
#         reducedareas += [reducedarea]

#     filename = "%d_hollow.npz" %(maxwind)
#     hcurrent, hlength, hpower, hresistance, hreducedarea = getdata(filename)

#     pl.figure(num=1, figsize=(11.69, 8.27))
#     pl.subplot(221)
#     pl.plot(csize, powers, colours[k]+'o-', label='W=%d' %maxwind)
#     pl.plot(11, hpower, colours[k]+'o')
#     pl.xlabel('AWG number')
#     pl.ylabel('Total power (W)')
#     pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
#     pl.xlim([0, 12])

#     pl.subplot(222)
#     pl.plot(csize, currents, colours[k]+'o-', label='W=%d' %maxwind)
#     pl.plot(11, hcurrent, colours[k]+'o')
#     pl.xlabel('AWG number')
#     pl.ylabel('Total current (I)')
#     pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
#     pl.xlim([0, 12])

#     pl.subplot(223)
#     pl.plot(csize, lengths, colours[k]+'o-', label='W=%d' %maxwind)
#     pl.plot(11, hlength, colours[k]+'o')
#     pl.xlabel('AWG number')
#     pl.ylabel('Total wire length (m)')
#     # pl.legend(loc='best')
#     pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
#     pl.xlim([0, 12])

#     pl.subplot(224)
#     pl.plot(csize, resistances, colours[k]+'o-', label='W=%d' %maxwind)
#     pl.plot(11, hresistance*1e3, colours[k]+'o')
#     pl.xlabel('AWG number')
#     pl.ylabel('Total resistance (mOhm)')
#     pl.legend(loc='best')
#     pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
#     pl.xlim([0, 12])

#     # pl.plot(csize, reducedareas, 'o-', label='W=%d' %maxwind)
#     # pl.xlabel('AWG number')
#     # pl.ylabel('Area reduction')

# pl.savefig('plotstudy.png')
# pl.show()


# # pl.plot(out)
# # pl.show()
