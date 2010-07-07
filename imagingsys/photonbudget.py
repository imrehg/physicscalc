#!/usr/bin/env python

# Photon budget calculation for the MOT imaging system

from __future__ import division
import numpy as np

def area(r):
    """ Area of a circle """
    return (np.pi * r**2)

def surface(r):
    """ Surface of sphere """
    return (4 * np.pi * r**2)

def radiate(I, d):
    """ Fluorescence rate of single Cs atoms
    Units are: Hz, mW/cm**2
    """
    Gamma = 5e6
    I0 = 3
    return (I / I0 * np.pi * Gamma) / (1 + I / I0 + 4*(d / Gamma)**2)
    

# Energy of a single photon in the 852nm range
c = 299792458
lam = 852e-9
Ep = 6.6e-34 * c / lam
# Approximate fluo rate of a single atom in our beam
R1 = radiate(I = 15, d = 10e6)
# Expected number of atoms
N = 1e5

# Detector area
A0 = area(1.25)
# Full 4pi area at the detector
A = surface(10)

# Conversion efficiency: A/W
C = 0.6 

## Power collected by area a0
p_rate = N * R1 * (A0 / A)
print "Photon rate: %f MHz" %(p_rate / 1e6)
print "Power rate : %f nW" %(p_rate * Ep * 1e9)
P = p_rate * Ep * C
print "Signal photocurrent: %f nA" %(P * 1e9)
# Gain of the detector
G = 1e10

## generated voltage
V = P * G
print "Generated voltage: %f mV" %(V*1000)
