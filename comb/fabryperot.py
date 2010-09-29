from __future__ import division
from numpy import *
import pylab as pl
from scipy.fftpack import *


def pulse(t, t0, a, w, f):
    return a*exp(-(t-t0)**2/(2*w**2))*sin(2*pi*t*f)

def filter(tdom, dt, fw, phi):
    fdom = fft(tdom)
    n = len(fdom)
    F1 = array(range(0,int(n/2)))/n/dt
    F = concatenate((F1, -1*F1[::-1]))
    filt = (sin(2*pi*fw*F + phi)+1)/2
    fdom2 = fdom*filt
    return F, fdom, fdom2, filt


t = linspace(-100, 100, 20000)
a = 1
w = 0.5
f = 3
fr = 10
reps = range(-95, 96, fr)
tdom = None
for i in reps:
    if tdom is None:
        tdom = pulse(t, i, a, w, f)
    else:
        tdom += pulse(t, i, a, w, f)

pl.figure(1)
pl.plot(t, tdom+0.5)

freqs = [fr/1.333, fr/1.5, fr/2, fr/3, fr/4]
for i, fi in enumerate(freqs):
    F, fdom, fdom2, filt = filter(tdom, t[1]-t[0], fi, pi/2)
    pl.figure(1)
    pl.plot(t, ifft(fdom2)-i-0.5)
    pl.figure(2)
    if (i == 0):
        pl.plot(F, abs(real(fdom))/max(abs(real(fdom)))+0.5)
    pl.plot(F, abs(real(fdom2))/max(abs(real(fdom)))-i-0.5)

pl.figure(1)
pl.xlim([-10, 10])
pl.xlabel('time')
pl.ylabel('electric field (+offset)')
pl.figure(2)
pl.xlim([2, 4])
pl.xlabel('frequency')
pl.ylabel('fourier mode power (+offset)')

# pl.figure(1)
# pl.savefig('fp_freqdomain.png')
# pl.figure(2)
# pl.savefig('fp_timedomain.png')

pl.show()
