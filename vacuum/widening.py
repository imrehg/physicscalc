#!/usr/bin/env python
from __future__ import division
from scipy import *
from pylab import *

# Vacuum system calculation for non uniform tube thickness
# Based on Building Scientific Apparatus, p96-97

# Average velocity of molecules (cm/s)
v = 2e4

# Length of tube (cm)
L = 2

# Diameter (cm), start, finish, dD/dL
D0 = 1
k = 2

# number of tube subdivisions
n = 400

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


def diameter(l, params):
    """ Tube diameter
    Standard is linear change, but can add arbitrary functions
    """
    # Linear
    (k, D0) = params
    D = k * l + D0

    # Sine
    #(k, D0) = params
    #D = D0 + 1/k*sin(l)
    
    return D

def sumconduct(L, n, D0, k):
    """ Total conductance of the vacuum tube

    Input:
     L : total length
     n : number of segments
     D0, k : tube geometry parameters, depending on diameter(),
             for the linear these are the starting diameter and 
             change/unit length, respectively

    Output:
     Cout: total conductance
    """
    # Subdivide the tube into smaller segments
    l = linspace(0, L, n+1)
    dl = diff(l)

    # Connect these segments into series
    # Cout = 1 / (Sum(1 / Ci))
    Cout = 0
    for i in range(0, n):
        li = l[i] + dl[i]/2
        params = (k, D0)
        D = diameter(li, params)
        Cout += 1.0 / conduct(D, dl[i], v)
    Cout = 1.0/Cout
    return Cout

####################################

# Series calculation along a whole tube as a function of its length
mn = 30
Llist = linspace(1,10,mn)
Con = zeros(mn)
DD = zeros(mn)
for m in range(0,mn):
    Con[m] = sumconduct(Llist[m], n, D0, k)
    #params = (k, D0)
    #DD[m] = diameter(Llist[m], params)


# Plot results
figure()
title('Sums')
plot(Llist, Con)
xlabel('Distance (cm)')
ylabel('Vacuum conductance (L/s)')

#figure()
#plot(Llist, DD)
#ylabel('Diameter')

show()
