"""
Laser linewidth calculation, based on Shevy (1993).
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integrate

def noise_white(h, t):
    return h/(2*t)

def integ(t, v, v0, h):
    s2 = noise_white(h, t) * v0 * v0 / 2
    return np.real(np.exp(-2 * np.pi**2 * t**2 * s2) * np.exp(2*np.pi*1j*(v-v0)*t))

# White noise parameter
h = 1e-26
# Carrier frequency
v0 = 300 * 1e12

# Depending on the value of "h", these two parameters might need adjustment
# The higher the h, the higher range dv has to be, and ilim (especially the upper limit)
# might need to be smaller
dv = np.linspace(-1e4, 1e4, 301)
ilim = [1e-7, 1e-1]

out = []
for d in dv:
    out += [2*integrate.quad(integ, ilim[0], ilim[1], args=(v0+d, v0, h))[0]]
out = np.array(out)
# Normalize by centre
out /= out[int(len(out)/2)+1]
pl.plot(dv, out)
pl.xlabel('Frequency detuning from carrier (Hz)')
pl.ylabel('Intensity (A.U.)')

pl.show()

