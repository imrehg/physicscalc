"""
Quick Allan variance calculation for the "allan_cont_long" type of data
(big stream of frequency records with constant time spacing)
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy as sp
import re
import os
import scipy.odr as odr
import ourgui

# Todo: rewrite it so does not have to use time data at all
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

def main(filename, gatetime=0.1, marker='k.',title=None):
    print filename

    data = np.loadtxt(filename, comments="#")
    n = len(data)
    time = np.array(range(0, n))*gatetime
    out = []
    while gatetime <= (time[-1]/3):
        print("Averaging time: %g" %(gatetime))
        out += [[gatetime, allan(time, data, gatetime, np.mean(data))]]
        gatetime *= 2
    out = np.array(out)

    pl.figure(1)
    pl.title(filename)
    pl.subplot(211)
    pl.plot(time, data-np.mean(data))
    pl.xlabel("Time (s)")
    pl.ylabel("Frequency deviation (Hz)")
    pl.xlim([0, time[-1]])

    pl.subplot(212)
    pl.loglog(out[:,0], out[:,1], marker)
    pl.xlabel("Averaging time (s)")
    pl.ylabel("Squared Allan variance")
    pl.xlim([out[0, 0], out[-1, 0]])

    pl.show()

if __name__ == "__main__":
    filename  = ourgui.openFile(type="log")
    main(filename, marker='k.-')
