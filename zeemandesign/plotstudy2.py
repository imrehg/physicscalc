import numpy as np
import pylab as pl

import wires
import zeemanslower as zs
from layeroptimize import *

# Whether or not separate the images
separate = True

# This the testing data to get the B-field scaling
series = 3
filename = "%d_AWG1_5.npz" %(series)
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

    volume = totallen * ((sim['wire'][0]/2)**2 * np.pi)
    radius = sim['wire'][0]/2
    return current, totallen, power, resistance, reducedarea, volume, radius

weightUnit = 8881.5  # Copper, kg/m^3

colours = ['b', 'g', 'r', 'y', 'k', 'm']
for k, maxwind in enumerate(range(5, 11)):
    powers = []
    currents = []
    lengths = []
    resistances = []
    reducedareas = []
    volumes = []
    weights = []
    currentdense = []
    csize = range(1, 16)
    for i in csize:
        filename = "%d_AWG%d_%d.npz" %(series, i, maxwind)
        current, totallen, power, resistance, reducedarea, volume, radius = getdata(filename)
        powers += [power]
        currents += [current]
        lengths += [totallen]
        resistances += [resistance]
        reducedareas += [reducedarea]
        volumes += [volume]
        weights += [volume * weightUnit]
        currentdense += [current / (radius*radius*np.pi)]

    if separate:
        fig = pl.figure(num=1, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        fig = pl.figure(num=1, figsize=(8.27, 11.69))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
        pl.subplot(321)
    pl.plot(csize, powers, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('AWG number')
    pl.ylabel('Total power (W)')
    pl.xlim([0, 16])

    if separate:
        fig = pl.figure(num=2, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        pl.subplot(322)
    pl.plot(csize, currents, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('AWG number')
    pl.ylabel('Total current (I)')
    pl.xlim([0, 16])

    if separate:
        fig = pl.figure(num=3, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        pl.subplot(323)
    pl.plot(csize, lengths, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('AWG number')
    pl.ylabel('Total wire length (m)')
    pl.xlim([0, 16])

    if separate:
        fig = pl.figure(num=4, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        pl.subplot(324)
    pl.plot(csize, resistances, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('AWG number')
    pl.ylabel('Total resistance (Ohm)')
    pl.legend(loc='best')
    pl.xlim([0, 16])

    if separate:
        fig = pl.figure(num=5, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        pl.subplot(325)
    pl.plot(csize, weights, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('AWG number')
    pl.ylabel('Weight (kg)')
    pl.xlim([0, 16])

    if separate:
        fig = pl.figure(num=6, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
    else:
        pl.subplot(326)
    pl.plot(np.array(currentdense)/1e6, powers, colours[k]+'o-', label='W=%d' %maxwind)
    pl.xlabel('Current density (A/mm^2)')
    pl.ylabel('Power (W)')

    if separate:
        pl.figure(num=7, figsize=(11.69, 8.27))
        fig.text(0.5, 0.95, 'Slower length = %g m, Max B field = %g G' %(sl, bfield*1e4), horizontalalignment='center')
        pl.plot(powers, weights, colours[k]+'o-')
        pl.xlabel('Power (W)')
        pl.ylabel('Total weight (kg)')

# Save figures
if separate:
    nfig = 7
    for i in range(1, nfig+1):
        pl.figure(i)
        pl.savefig('plotstudy_%d_%d.png' %(series, i))
        pl.savefig('plotstudy_%d_%d.pdf' %(series, i))
else:
    pl.savefig('plotstudy_%d.png' %(series))
    pl.savefig('plotstudy_%d.pdf' %(series))
pl.show()
