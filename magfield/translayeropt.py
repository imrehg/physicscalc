#!/usr/bin/env python
# Some code based on Sumner, et. a.l, JPhysD, 20 (1987), 1095-1101

from __future__ import division
from scipy.optimize import fmin
from translayer import *

def findlayer(dr, *args):
    """ Function to optimize: three-layer shielding, how to place outer layers """
    u, t, R = args
    R[1] = dr[0]
    R[2] = dr[1]
    Sti = [u*t/(2*r) for r in R]
    return -getscreening(Sti, R)

## Settings for the simulation
# relative permittivity
u = 1e3
# thickness
t = 0.001
# Radiuses : but not all will be used
R = [0.1, 0.12, 0.2]

p0 = [0.1, 0.12]
res = fmin(findlayer, p0, args=(u, t, R), full_output=1)
print "Layer position: ", res[0]
print "Total shielding factor: ", res[1]
