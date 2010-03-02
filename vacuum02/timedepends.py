#!/usr/bin/env python

# Vacuum system pressure during dispenser operation
# Simplistic modeling of our vacuum system

from numpy import *
import pylab as pl
from scipy import integrate

# Definition of parameters
# pumping speed
ionpump = 0.25
# dispenser loading rate
dispenser_on = 0.6
# dispenser off time
tdisp = 20
# diffusion rate between the two chambers
diffuse = 0.01

def dispenser(t):
    """ Time dependent dispenser.
    If t < tdisp, on, otherwise of """
    if t < tdisp:
        return dispenser_on
    else:
        return 0

def dn_dt(X, t=0):
    """ Differential equation for the number of atoms in the two chambers 
    Chamber 1: +dispenser, -ionpump, +- diffusion
    Chamber 2: -+ diffusion 
    """
    return array([ dispenser(t) - ionpump * X[0] - diffuse * (X[0] - X[1]),
                   diffuse * (X[0] - X[1])])

t = linspace(0, 100, 5000)
n0 = array([1, 1.2])
n, infodict = integrate.odeint(dn_dt, n0, t, full_output=True)
infodict['message']   

ch1, ch2 = n.T

pl.plot(t, ch1, label = "Chamber 1")
pl.plot(t, ch2, label = "Chamber 2")
pl.xlabel("Time")
pl.ylabel("Pressure")
pl.legend(loc = "best")
pl.show()
