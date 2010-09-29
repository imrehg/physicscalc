from __future__ import division
from numpy import *
import pylab as pl
from scipy.fftpack import *


def pulse(t, t0, a, w, f):
    return a*exp(-(t-t0)**2/(2*w**2))*sin(2*pi*t*f)

def filtersin(tdom, dt, fw, phi):
    """ 
    Sinusoid filter in frequency domain
    filtersin(tdom, dt, fw, phi)

    Input
    =====
    tdom : time domain signal
    dt : time step size
    fw : filter sine function frequency
    phi: filter sine function phase offset 
         (probably only n*pi/2 n=0,1.. makes sense in this implementation)

    Output
    ======
    F : frequency components
    fdom : fft(tdom)
    fdomfilt : fdom*filt
    filt : filter function
    """
    fdom = fft(tdom)
    n = len(fdom)
    F1 = array(range(0,int(n/2)))/n/dt
    F = concatenate((F1, -1*F1[::-1]))
    filt = (sin(2*pi*fw*F + phi)+1)/2
    fdomfilt = fdom*filt
    return F, fdom, fdomfilt, filt

def filterexp(tdom, dt, fp):
    """ 
    Exponential filter in frequency domain
    filterexp(tdom, dt, fp)

    Input
    =====
    tdom : time domain signal
    dt : time step size
    fp : filter parameters in the form of
         [(fc1, fw1), (fc2, fw2), ...]
         where fcX and fwX is the Xth exponential's centre and width respectively

    Output
    ======
    F : frequency components
    fdom : fft(tdom)
    fdomfilt : fdom*filt
    filt : filter function
    """
    fdom = fft(tdom)
    n = len(fdom)
    F1 = array(range(0,int(n/2)))/n/dt
    F = concatenate((F1, -1*F1[::-1]))
    filt = None
    for pars in fp:
        fc = pars[0]
        fw = pars[1]
        if filt is None:
            filt = exp(-(F-fc)**2/(2*fw**2))
            filt += exp(-(F+fc)**2/(2*fw**2))
        else:
            filt += exp(-(F-fc)**2/(2*fw**2))
            filt += exp(-(F+fc)**2/(2*fw**2))

    fdomfilt = fdom*filt
    return F, fdom, fdomfilt, filt

def testsin(savefigs = False):
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
        F, fdom, fdom2, filt = filtersin(tdom, t[1]-t[0], fi, pi/2)
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

    if savefigs:
        pl.figure(1)
        pl.savefig('fpsin_freqdomain.png')
        pl.figure(2)
        pl.savefig('fpsin_timedomain.png')

def testexp(savefigs = False):
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

    pl.figure(11)
    pl.plot(t, tdom+0.5)

    filters = [[(3.1, 0.01), (2.9, 0.01)], [(3.2, 0.01), (2.8, 0.01)], [(3.3, 0.01), (2.7, 0.01)]]
    for i, fi in enumerate(filters):
        F, fdom, fdom2, filt = filterexp(tdom, t[1]-t[0], fi)
        pl.figure(11)
        pl.plot(t, ifft(fdom2)-i-0.5)
        pl.figure(12)
        if (i == 0):
            pl.plot(F, abs(real(fdom))/max(abs(real(fdom)))+0.5)
        pl.plot(F, abs(real(fdom2))/max(abs(real(fdom)))-i-0.5)

    pl.figure(11)
    pl.xlim([-10, 10])
    pl.xlabel('time')
    pl.ylabel('electric field (+offset)')
    pl.figure(12)
    pl.xlim([2, 4])
    pl.xlabel('frequency')
    pl.ylabel('fourier mode power (+offset)')

    if savefigs:
        pl.figure(11)
        pl.savefig('fpexp_freqdomain.png')
        pl.figure(12)
        pl.savefig('fpexp_timedomain.png')


if __name__ == '__main__':
    savefigs = False
    testsin(savefigs)
    testexp(savefigs)
    pl.show()
