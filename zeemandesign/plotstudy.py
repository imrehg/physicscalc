import numpy as np
import pylab as pl

import wires
import zeemanslower as zs
from layeroptimize import *



filename = "5_AWG1.npz"
sim = np.load(filename)['simulation'][()]
eta, v0, vf, detu = sim['eta'], sim['v0'], sim['vf'], sim['detu']

atom = zs.Rb85()
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

# getdata("10_AWG7.npz")

separate = False

colours = ['b', 'g', 'r', 'y', 'k', 'm']
for k, maxwind in enumerate(range(5, 11)):
    powers = []
    currents = []
    lengths = []
    resistances = []
    reducedareas = []
    csize = range(1, 11)
    for i in csize:
        # maxwind = 5
        filename = "%d_AWG%d.npz" %(maxwind, i)
        current, totallen, power, resistance, reducedarea = getdata(filename)
        powers += [power]
        currents += [current]
        lengths += [totallen]
        resistances += [resistance*1e3]
        reducedareas += [reducedarea]

    filename = "%d_hollow.npz" %(maxwind)
    hcurrent, hlength, hpower, hresistance, hreducedarea = getdata(filename)

    if separate:
        pl.figure(num=1, figsize=(11.69, 8.27))
    else:
        pl.figure(num=1, figsize=(11.69, 8.27))
        pl.subplot(221)
    pl.plot(csize, powers, colours[k]+'o-', label='W=%d' %maxwind)
    pl.plot(11, hpower, colours[k]+'o')
    pl.xlabel('AWG number')
    pl.ylabel('Total power (W)')
    pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
    pl.xlim([0, 12])

    if separate:
        pl.figure(num=2, figsize=(11.69, 8.27))
    else:
        pl.subplot(222)
    pl.plot(csize, currents, colours[k]+'o-', label='W=%d' %maxwind)
    pl.plot(11, hcurrent, colours[k]+'o')
    pl.xlabel('AWG number')
    pl.ylabel('Total current (I)')
    pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
    pl.xlim([0, 12])

    pl.subplot(223)
    pl.plot(csize, lengths, colours[k]+'o-', label='W=%d' %maxwind)
    pl.plot(11, hlength, colours[k]+'o')
    pl.xlabel('AWG number')
    pl.ylabel('Total wire length (m)')
    # pl.legend(loc='best')
    pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
    pl.xlim([0, 12])

    pl.subplot(224)
    pl.plot(csize, resistances, colours[k]+'o-', label='W=%d' %maxwind)
    pl.plot(11, hresistance*1e3, colours[k]+'o')
    pl.xlabel('AWG number')
    pl.ylabel('Total resistance (mOhm)')
    pl.legend(loc='best')
    pl.xticks( csize[0::]+[11],  (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'hollow') )
    pl.xlim([0, 12])

    # pl.plot(csize, reducedareas, 'o-', label='W=%d' %maxwind)
    # pl.xlabel('AWG number')
    # pl.ylabel('Area reduction')

pl.savefig('plotstudy.png')
pl.show()


# # pl.plot(out)
# # pl.show()
