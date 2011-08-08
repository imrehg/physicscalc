#!/usr/bin/env python

from __future__ import division
from numpy import *
# from scipy import *
from scipy.odr import *
import pylab as pl
from scipy.signal import convolve
from numpy import convolve
from scipy.special import wofz
from scipy.interpolate import interp1d


########### Functional forms to fit

laplace = lambda p, x: p[2]*exp(- abs(x-p[0])/p[1])
lorentz = lambda p, x: p[2]*p[1]/((x-p[0])**2 + p[1]**2)
gauss = lambda p, x: p[2]/p[1]*exp(-(x-p[0])**2/(2*p[1]**2))
voigt = lambda p, x: p[2]*real(wofz(((x-p[0])+1j*p[1])/(p[3]*sqrt(2))))/(p[3]*sqrt(2*pi))

def multilorentz(p, x):
    '''
    multiple Lorentz distributions

    p : parametes in formatt [x01, g1, i1, x02, g2, ...., in, c]
        with centre, width, height for n Lorentzians + c constans offset
    '''
    n = int((len(p)-1)/3)
    res = zeros(len(x))
    for i in xrange(n):
        params = p[(i*3):((i+1)*3)]
        res += lorentz(params, x)
    res += p[-1]
    return res

def multivoigt(p, x):
    '''
    multiple Lorentz distributions

    p : parametes in formatt [x01, g1, i1, x02, g2, ...., in, s, c]
        with centre, width, height for n Lorentzians 
        + s gaussian width, c constans offset
    '''
    n = int((len(p)-2)/3)
    res = zeros(len(x))
    for i in xrange(n):
        params = append(p[(i*3):((i+1)*3)], p[-2])
        res += voigt(params, x)
    res += p[-1]
    return res

def multivoigt_g(p, x):
    '''
    multiple Lorentz distributions

    p : parametes in formatt [x01, i1, x02, ...., in, g, s, c]
        with centre, height for n Lorentzians 
        + g lorentz width, s gaussian width, c constans offset
    '''
    n = int((len(p)-2)/3)
    res = zeros(len(x))
    for i in xrange(n):
        params = append(p[(i*3):((i+1)*3)], p[-2])
        res += voigt(params, x)
    res += p[-1]
    return res


def getu(T):
    """
    Root-mean-square velocity: u = sqrt(2*kB*T/m) [m/s]
    input: T (C)
    """
    kB = 1.38e-23
    M = 133*1.66e-27
    T = 273+65
    u = sqrt(2*kB*T/M)
    return u

def transitfunc(x, w0, Ge):
    """ Transit broadening profile - convolutuion """
    detulist = linspace(-20e6, 20e6, 601)*2*pi
    nh = 301
    u = getu(70)
    # Beam radius [m]
    w0 = w0*1e-6
    # Upper level linewidth
    Ge = Ge*1e6*2*pi
    laplaced = lambda x, u, b: exp(- abs(x-u)/b)
    lorentzd = lambda x, u, g: g/2/((x-u)**2 + g*g/4)

    a = laplaced(detulist, 0, u/w0)/laplaced(0,0,u/w0)
    b = lorentzd(detulist, 0, Ge)/lorentzd(0, 0, Ge)
    o1 = scipy.signal.convolve(a, b, mode=1)
    o1 = o1/o1[nh]
    f = interp1d(detulist, o1)
    return f(x*1e6*2*pi)

def double(p, xlim, np):
    """ Double convolution profile """
    G, u = p
    xlims = linspace(-xlim, xlim, np)
    pa = [0, u, 1]
    a = laplace(pa, xlims)/laplace(pa, 0)
    pb = [0, G, 1]
    b = lorentz(pb, xlims)/lorentz(pb, 0)
    o1 = convolve(a, b, mode=1)
    o1 = o1/o1[(np-1)/2]
    return o1

def tripple(p, xlim, np):
    """ Tripple convolution profile """
    G, u, s = p
    xlims = linspace(-xlim, xlim, np)

    pa = [0, u, 1]
    a = laplace(pa, xlims)/laplace(pa, 0)
    pb = [0, G, 1]
    b = lorentz(pb, xlims)/lorentz(pb, 0)
    pc = [0, s, 1]
    c = gauss(pc, xlims)/gauss(pc, 0)
    o1 = convolve(a, b, mode=1)
    o2 = convolve(o1, c, mode=1)
    o2 = o2/o2[(np-1)/2]
    return o2

def multipledouble(p, x):
    ''' Transit broadedning '''
    n = int((len(p)-3)/2)
    xlim = x[-1] - x[0]
    xc = linspace(-xlim, xlim, xlim*3+1)
    np = xlim*3+1
    res = zeros(len(x))
    for i in xrange(n):
        x0N = p[i*2]
        iN = p[i*2+1]
        G = p[-3]
        u = p[-2]
        params = (G, u)
        t = double(params, xlim, np)
        o2 = interp1d(xc, t)
        res += o2(x-x0N)*iN
    res += p[-1]
    return res

def multipledouble_plusg(p, x):
    n = int((len(p)-4)/3)
    xlim = x[-1] - x[0]
    xc = linspace(-xlim, xlim, xlim*3+1)
    np = xlim*3+1
    res = zeros(len(x))
    s = p[-3]
    A = p[-4]
    for i in xrange(n):
        x0N = p[i*3]
        iN = p[i*3+1]
        G = p[i*3+2]
        u = p[-2]
        params = (G, u)
        t = double(params, xlim, np)
        o2 = interp1d(xc, t)
        res += o2(x-x0N)*iN
        res += gauss([x0N, s, iN*A], x)
    res += p[-1]
    return res

# def multipledouble_gb(p, x):
#     ''' Transit broadedning with gaussian background'''
#     n = int((len(p)-5)/2)
#     xlim = x[-1] - x[0]
#     xc = linspace(-xlim, xlim, xlim*3+1)
#     np = xlim*3+1
#     res = zeros(len(x))
#     for i in xrange(n):
#         x0N = p[i*2]
#         iN = p[i*2+1]
#         G = p[-5]
#         u = p[-4]
#         params = (G, u)
#         t = double(params, xlim, np)
#         o2 = interp1d(xc, t)
#         res += o2(x-x0N)*iN
#     res += gauss([p[-3], p[-2], p[-1]], x)
#     return res


def multipledouble_freeg(p, x):
    ''' Transit broadening with G as free parameter '''
    n = int((len(p)-2)/3)
    xlim = x[-1] - x[0]
    # xc = linspace(-xlim, xlim, xlim*3+1)
    xc = linspace(-xlim, xlim, 501)
    np = len(xc)
    res = zeros(len(x))
    for i in xrange(n):
        x0N = p[i*3]
        iN = p[i*3+1]
        G = p[i*3+2]
        u = p[-2]
        params = (G, u)
        t = double(params, xlim, np)
        o2 = interp1d(xc, t)
        res += o2(x-x0N)*iN
    res += p[-1]
    return res


def multipletripple(p, x):
    ''' Tripple convolution with single width for all '''
    n = int((len(p)-4)/2)
    xlim = x[-1] - x[0]
    xc = linspace(-xlim, xlim, xlim*10+1)
    np = len(xc)
    res = zeros(len(x))
    for i in xrange(n):
        x0N = p[i*2]
        iN = p[i*2+1]
        s = p[-4]
        G = p[-3]
        u = p[-2]
        params = (G, u, s)
        t = tripple(params, xlim, np)
        o2 = interp1d(xc, t)
        res += o2(x-x0N)*iN
    res += p[-1]
    return res

def multipletripple_freeg(p, x):
    ''' Tripple convolution with flexible Lorentzian width '''
    n = int((len(p)-3)/3)
    xlim = x[-1] - x[0]
    xlim = xlim*2
    xc = linspace(-xlim, xlim, 1201)
    np = len(xc)
    # xc = linspace(-xlim, xlim, xlim*3+1)
    # np = xlim*3+1
    res = zeros(len(x))
    for i in xrange(n):
        x0N = p[i*3]
        iN = p[i*3+1]
        G = p[i*3+2]
        s = p[-3]
        u = p[-2]
        params = (G, u, s)
        t = tripple(params, xlim, np)
        o2 = interp1d(xc, t, 'linear')
        res += o2(x-x0N)*iN
    res += p[-1]
    return res
    

####################

def dofit(f, x, y, p0):
    ''' Do actual fitting of the data '''
    data = Data(x, y)
    model = Model(f)
    fit = ODR(data, model, p0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output

def completefit(f, x, y, p0, filename, name, toplot=False, tosave=False):
    print "##### %s:" %(name)
    out = dofit(f, x, y, p0)
    for i in range(len(out.beta)):
        print "Param %d: %g (%g)" %(i, out.beta[i], out.sd_beta[i])

    yh = f(out.beta, x)
    if tosave:
        fitpars = zip(out.beta, out.sd_beta)
        savetxt("%s_%s.txt"%(filename, name), fitpars, delimiter=',')
        fitdata = zip(x, y, yh)
        savetxt("%s_%s.predict.txt"%(filename, name), fitdata, delimiter=',')
    if toplot: 
        pl.figure()
        pl.plot(x, y, 'k.')
        pl.plot(x, yh, 'r-', linewidth=2)
        pl.title('%s'%name)
        pl.xlabel('x')
        pl.ylabel('Data & fit')
        pl.xlim([min(x), max(x)])
        if tosave:
            pl.savefig("%s_%s.png"%(filename, name))
            pl.savefig("%s_%s.eps"%(filename, name))
        pl.figure()
        pl.plot(x, y-yh, 'k.')
        pl.title('%s'%name)
        pl.xlabel('x')
        pl.ylabel('Residuals')
        pl.xlim([min(x), max(x)])
        if tosave:
            pl.savefig("%s_%s-residuals.png"%(filename, name))
            pl.savefig("%s_%s-residuals.eps"%(filename, name))
    return yh, out

