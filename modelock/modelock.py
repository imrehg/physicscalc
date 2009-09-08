#!/usr/bin/env python
# Calculation of the effect of summing up a lot of
# sine functions. Models what happens in the case
# of mode-locking.
# Also, see the effect of phase-difference between
# the different modes.

from scipy import *
from pylab import *

# number of time points in time-domanin sim
np = 10000
# base frequency (omega)
w = 1
# base phase
ph = pi/2
# max number of modes - highest number in summation
maxmode = 10

def moder(maxmode,t,w,ph,phasemode,phr):
    ''' Mode-lock squence calculator '''
    s1 = sin(w*t+ph)
    for n in range(2,maxmode+1):
        if phasemode == "rand":
            ph1 = ph + phr * (0.5 - rand())
        elif phasemode == "diff":
            ph1 = ph + n * phr
        else:
            ph1 = ph
        s1 += sin(n*w*t+ph1)
    return s1


# Perfect mode locking
t = linspace(0,2*pi/w*2,np)
s1 = moder(maxmode,t,w,ph,"",0)
figure()
plot(t/(2*pi/w),sqrt(s1*s1))
# plot(t/(2*pi/w),s1)
xlim([0,t[-1]/(2*pi/w)])
title("Perfect mode locking, modes: %d" %(maxmode))
xlabel("Time (periods of base)")
ylabel("Intensity")

# Modes with random phase
t = linspace(0,2*pi/w*2,np)
s1 = moder(maxmode,t,w,ph,"rand",2*pi)
figure()
plot(t/(2*pi/w),sqrt(s1*s1))
xlim([0,t[-1]/(2*pi/w)])
title("Random phase, modes: %d" %(maxmode))
xlabel("Time (periods of base)")
ylabel("Intensity")


#
# Next two only makes sense when phase = pi/2 (cosine)
#

# Small phase difference between modes
t = 2*pi/w
phrlist = linspace(0,2*pi,1000)
maxphr = zeros(phrlist.size)
for phr in range(0,phrlist.size):
    s1 = moder(maxmode,t,w,ph,"diff",phrlist[phr])
    maxphr[phr] = s1
figure()
plot(phrlist,maxphr)
title("Small constant phase difference, modes: %d" %(maxmode))
xlabel("Phase difference (radian)")
ylabel("Mod-locked maximum value")
xlim([0,phrlist[-1]])

# Random phase of modes
t = 2*pi/w
phrlist = linspace(0,2*pi,1000)
maxphr = zeros(phrlist.size)
for phr in range(0,phrlist.size):
    s1 = moder(maxmode,t,w,ph,"rand",phrlist[phr])
    maxphr[phr] = s1
figure()
plot(phrlist,maxphr)
title("Random phase, modes: %d" %(maxmode))
xlabel("Max phase range (radian)")
ylabel("Mod-locked maximum value")
xlim([0,phrlist[-1]])

# show plots
show()
