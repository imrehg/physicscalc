from __future__ import division
import numpy as np

hbar = 1.05457148e-34
uB = hbar * 2 * np.pi * 1.399e6

B = 150

## This is for RB87
gP3 = 1.3362
gS1 = 2.002
me = -4
mg = -3
uprime = uB * (gP3*me - gS1*mg)

detu = uprime * B / hbar
print detu / 2 / np.pi / 1e6


##### Units
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
tU = 1/GRb # Time: 1/linewidth
xU = 1/kRb # Distance: 1/wavenumber

# Dimensionless detuning
print detu * tU
