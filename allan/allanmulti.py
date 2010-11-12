from __future__ import division
import numpy as np
import pylab as pl
import scipy as sp
import re
import os
import scipy.odr as odr

def allan(t, freq, tau, base):
    """
    allan(t, y, tau, base)
    Allan variance calculation

    Input variables:
    ----------------
    t    : time of measurement
    freq : measured frequency
    tau  : averaging time
    base : base frequency

    Output variables:
    -----------------
    s : Squared Allan variance
    """
    # Divide time up to 'tau' length units for averaging
    times = np.arange(min(t), max(t), tau)
    # Create temporary variable for fractional frequencies
    vari = np.zeros(len(times))
    for tstep in range(0, len(times)):
        # Get the data within the time interval
        data = freq[(t >= times[tstep]) & (t < (times[tstep] + tau))]
        # Fractional frequency calculation
        vari[tstep] = (sp.mean(data) - base) / base
    # Squared Allan variance
    s = sp.mean((vari[0:-1] - vari[1:]) ** 2) / 2
    return s


def noisefit(h, t):
    for i in xrange(len(h)):
        if h[i] < 0:
            h[i] = 0
    return h[0]*t**-3 + h[1]*t**-2 + h[2]*t**-1 + h[3] + h[4]*t + h[5]*t**2

def dofit(x, y, h0):
    data = odr.Data(x, y)
    model = odr.Model(noisefit)
    fit = odr.ODR(data, model, h0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output

path = '.'
dirList=os.listdir(path)

avgtime = []
allanv = []
for fname in dirList:
    m = re.match(r"(?P<header>\w+?)_(?P<date>\d{6})_(?P<time>\d{6})_(?P<parameter>[0-9.]+)\.log", fname)
    if m:
        out = m.groupdict()
        out['parameter'] = float(out['parameter'])
        out['filename'] = fname
        
        data = np.loadtxt(out['filename'])
        logtime = data[:, 0]
        freq = data[:, 1]
        # allanv += [np.mean(abs(np.diff(freq))/np.mean(freq))]
        allanv += [np.sqrt(np.mean(0.5*(np.diff(freq/np.mean(freq)))**2))]

        avgtime += [out['parameter']]



pl.loglog(avgtime, allanv, '.', markersize=10)
t = np.logspace(-1, 2, 1001)
h = [1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]

out = dofit(avgtime, allanv, h)
fitbeta = out.beta
for i in xrange(len(fitbeta)):
    fitbeta[i] = max(0, fitbeta[i])
names = ['White PM', 'Flicker PM', 'White FM', 'Flicker FM', 'Random Walk', 'Drift']


pl.loglog(t, noisefit(fitbeta, t))
total = np.zeros(len(t))
for i in xrange(len(h)):
    if fitbeta[i] > 0:
        hx = [0]*6
        hx[i] = out.beta[i]
        nf = noisefit(hx, t)
        total += nf
        pl.loglog(t, nf, 'k--')
        print "%s : %g" %(names[i], hx[i])
    else:
        # print "%d is wrong range" %(i)
        print "%s don't seem to play a role" %(names[i])

pl.loglog(t, total, 'r:')
pl.ylim([1e-10, 1e-8])
pl.xlabel('averaging time (s)')
pl.ylabel('Allan deviation $\sigma$')
pl.show()
