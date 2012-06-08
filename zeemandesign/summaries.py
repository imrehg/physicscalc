from __future__ import division
import numpy as np
import pylab as pl
import csv

import layeroptimize as lo
import zeemanslower as zs

fw, fh = 11.69, 8.26  # A4 size in inches

## Settings for this run
sumname = "summary01"
series = [15,
          16,
          17,
          18,
          19
          ]
lnums = range(5, 16)

nseries = len(series)
nnums = len(lnums)

def savecsv(filename, savearray, headers, colhead):
    """ Save data structure in CSV file with column and row headers"""
    resout = csv.writer(open(filename, "wb"),
                        delimiter=' ',
                        quotechar='|',
                        quoting=csv.QUOTE_MINIMAL
                        )
    resout.writerow(headers)
    row, col = np.shape(savearray)
    for i in range(row):
        resout.writerow([colhead[i]]+list(savearray[i, :]))

weightUnit = 8881.5  # Copper, kg/m^3
atom = zs.Rb87()

sizes = np.zeros((nnums, nseries))
currents = np.zeros((nnums, nseries))
powers = np.zeros((nnums, nseries))
weights = np.zeros((nnums, nseries))
deco = ['ko-', 'ro--', 'bo:', 'yo-', 'go--']
pl.figure(figsize=(fw, fh))
for i, s in enumerate(series):
    for j, l in enumerate(lnums):
        simfile = "%d/%d_AWG12Coat_%d.npz" %(s, s, l)
        sim = np.load(simfile)['simulation'][()]
        wire = sim['wire']
        setup = sim['setup']
        eta = sim['eta']
        v0 = sim['v0']
        vf = sim['vf']
        looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']
        eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']

        ## Slower length
        d = wire[0]/2  # half wire diameter
        l1 = segments[1][0]
        l2 = segments[-2][1]-1
        slength = (looppos[l2]+d) - (looppos[l1]-d)
        sizes[j, i] = slength * 100  # turn into cm

        ## Current
        bideal = zs.bideal(atom, 0, eta, v0, vf, detu)
        bthis = lo.fieldcalc(0, setup)
        I0 = bideal / bthis
        currents[j, i] = I0

        ## Power
        totallen = lo.totalcoillength(sim['setup'], sim['R'], sim['wire'][0])
        resistance = sim['wire'][1] * totallen
        power = resistance * I0**2
        powers[j, i] = power

        ## Weight
        volume = totallen * ((sim['wire'][0]/2)**2 * np.pi)
        weight = volume * weightUnit
        weights[j, i] = weight

    ## Plot the Power vs. Weight curve
    pl.plot(weights[:, i], powers[:, i], deco[i], label='Series %d' %s, linewidth=2)
    pl.legend(loc='best')
    pl.xlabel('Slower weight [kg]', fontsize=15)
    pl.ylabel('Power [W]', fontsize=15)
pl.savefig("%s.png" %sumname)
pl.savefig("%s.pdf" %sumname)

# Save data
savecsv("%s_length.csv" %(sumname), sizes, ["Layers\Length(cm)"]+["Series %s" %s for s in series], lnums)
savecsv("%s_current.csv" %(sumname), currents, ["Layers\Current(A)"]+["Series %s" %s for s in series], lnums)
savecsv("%s_power.csv" %(sumname), powers, ["Layers\Power(W)"]+["Series %s" %s for s in series], lnums)

# Finally
pl.show()
