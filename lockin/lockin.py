
from __future__ import division
from numpy import *
from scipy.integrate import *
import pylab as pl

def reference(t, *args):
    f, phi = args
    return sin(2*pi*f*t + phi)
    
def signal(t):
    # fs = fl*0.9 
    # return sin(2*pi*fs*t)
    # ee = 1
    # return (1-mod(floor(t),2))*(1-exp(-ee*mod(t,1))) + mod(floor(t),2)*(exp(-ee*mod(t, 1))-exp(-ee))
    return 1 + mod(floor(t),2)*0.1

def lockin(t, *args):
    return reference(t, *args)*signal(t)


phl = linspace(0,90,6)/180*pi
for phi in phl:
    T = 10
    fl = 0.5
    # phi = 20/180*pi
    t = linspace(0, 10, 45)
    out = array([quadrature(lockin, tt-T, tt, args=(fl, phi), maxiter=60) for tt in t])/T
    pl.plot(t, out[:,0], label='%d' %(phi/pi*180))

pl.legend(loc='best')
pl.show()
