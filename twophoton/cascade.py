#!/usr/bin/env python
# Cascade configuration two photon absorption calc
# Based on: DOI: 10.1103/PhysRevLett.100.203001

from numpy import * 
from pylab import *
from time import time

### Adjustable variables
# Wavelength of levels
# [6P, 8S]
wavellevels = array([852e-9, 794e-9])
# decay rates / 2
decayrate = array([5e6, 1e6])*2*pi / 2
# Repetition frequency
#freqrep = 80e6
# Modes to include
modes = array([-100, 100])
########################

### Constants and derived variables
# Speed of light
c = 299792458
# Frequency of levels (angular untits)
wlevels = c/wavellevels*2*pi
# Center frequency:
fcentre = sum(wlevels)/2 / (2*pi)
# 2-level freqency in angular units: centre locked on it
wtwo = 2 * (2 * pi * fcentre) 
########################

twopifcentre = 2 * pi * fcentre
fourpifcentre = 2 * twopifcentre
offs = wlevels[0] - twopifcentre
offs2 = wtwo - fourpifcentre
decay = decayrate[0]

def cmode(params, x):
    (wtwo, fcentre, freqrep, decayrate, (m, n)) = params
    return 1 / (1j * (offs - x * 2 * pi * freqrep) + decay)

def cgf(params):
    (wtwo, fcentre, freqrep, decayrate, (m, n)) = params
    return 1 / (1j * (offs2 - (m + n) * 2 * pi * freqrep) +
        decayrate[1]) * (cmode(params, n) + cmode(params, m))

def gettrans(params):
    (wtwo, fcentre, freqrep, decayrate, (m, n)) = params
    modesn  = floor((wlevels/(2*pi) - fcentre) / freqrep)
    goodmode = abs(int(modesn[0]))
    ctot = 0j
    for m in range(goodmode-11,goodmode+11):
        for n in range(goodmode-11, goodmode+11):
            params = (wtwo, fcentre, freqrep, decayrate, (-m, n))
            ctot += cgf(params)
    return ctot

def scan(freqrepscan):
    trans = []
    for f in freqrepscan:
        params = (wtwo, fcentre, f, decayrate, (0, 0))
        trans.append(abs(gettrans(params)))
    return trans

# Do a scan of frequency
freqscan = linspace(90e6-1e3,90e6+1e3,401)
start = time()
res = scan(freqscan)
elapsed = time() - start
print "Elapsed time:", elapsed
semilogy((freqscan-90e6)/1e3, res, '.-')
xlabel('Repetition frequency - 90 Mhz (kHz)')
ylabel('|Ctot|')
show()
