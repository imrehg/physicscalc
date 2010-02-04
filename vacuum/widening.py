#!/usr/bin/env python
from __future__ import division
from scipy import linspace, diff

# Vacuum system calculation for non uniform tube thickness
# Based on Building Scientific Apparatus, p96-97

# Average velocity of molecules (cm/s)
v = 2e4

# Length of tube (cm)
L = 5

# Diameter (cm), start, finish, dD/dL
D0 = 0.5
Dt = 4
k = (Dt - D0) / L

# number of tube subdivisions
n = 10000


def conduct(D, l, v):
    """ Vacuum conductance

    Input:
        D = tube diameter, cm
        l = tube segment length, cm
        v = average moleculat velocity, cm/s

    Output:
        Conductance
    """
    return 2.6e-4 * v * D ** 3 / l


# Subdivide the tube into smaller segments
l = linspace(0, L, n+1)
dl = diff(l)

# Connect these segments into series
# Cout = 1 / (Sum(1 / Ci))
Cout = 0
for i in range(0, n):
    li = l[i] + dl[i]/2
    D = k * li + D0
    Cout += 1.0 / conduct(D, dl[i], v)
Cout = 1.0/Cout

print "Total conductance with %d tube segments: %.3f L/s" %(n, Cout)
