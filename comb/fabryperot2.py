from __future__ import division
from numpy import *
import pylab as pl
from scipy.fftpack import *


def pulse(t, t0, a, w, f):
    return a*exp(-(t-t0)**2/(2*w**2))*sin(2*pi*t*f)

def filter2(tdom, dt, fp):
    fdom = fft(tdom)
    n = len(fdom)
    F1 = array(range(0,int(n/2)))/n/dt
    F = concatenate((F1, -1*F1[::-1]))
    filt = None
    for pars in fp:
        print pars
        fc = pars[0]
        fw = pars[1]
        if filt is None:
            filt = exp(-(F-fc)**2/(2*fw**2))
            filt += exp(-(F+fc)**2/(2*fw**2))
        else:
            filt += exp(-(F-fc)**2/(2*fw**2))
            filt += exp(-(F+fc)**2/(2*fw**2))

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
F, fdom, fdom2, filt = filter2(tdom, t[1]-t[0], [(3.3, 0.01), (2.7, 0.01)])
pl.figure(1)
pl.plot(t, ifft(fdom2)-0.5)
pl.figure(2)
pl.plot(F, abs(real(fdom))/max(abs(real(fdom)))+0.5)
pl.plot(F, abs(real(fdom2))/max(abs(real(fdom2)))-0.5)

pl.figure(1)
pl.xlim([-10, 10])
pl.xlabel('time')
pl.ylabel('electric field (+offset)')
pl.figure(2)
pl.xlim([2, 4])
pl.xlabel('frequency')
pl.ylabel('fourier mode power (+offset)')

pl.figure(1)
pl.savefig('fp2_freqdomain.png')
pl.figure(2)
pl.savefig('fp2ls_timedomain.png')

pl.show()
