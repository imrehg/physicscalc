from __future__ import division
from numpy import *
from scipy.integrate import odeint, simps
from scipy.optimize import leastsq
from pylab import plot, show, legend, figure, title, xlabel, ylabel, savefig

## Input parameters
# filter frequency, in modulation frequency units
fc = 1
# Background levels: [input-independent, input-dependent]
bglevel = [0.1, 1]
# lock-in time constant
Ts = 6
# Signal model:
# 0: original low-pass
# 1: slow build up, quick decay, (similar to CPT)
model = 1
# Save figures
savefigures = True
savefigbase = "cpt02_"
savefigext = "pdf"


# Helper functions
ph2deg = lambda ph: ph/pi*180
chop = lambda t, a, b : a if mod(ceil(t*2),2) > 0 else b

# Filter response calculation
def filterresponse(t, fc):
    """ Calculate phase shift for low pass filter response to square wave
    fc: filter frequency, in chopped frequency units
    toplot: show calculation on plot? (True/False)

    returns phase shift in degrees
    """

    force = lambda t,f: mod(ceil(t*f*2),2)
    filterode = lambda u,t,b,f: -u[0]/b + force(t,f)

    u0 = array([0,1])
    b = 1/(2*pi*fc)
    f = 1.0
    u = odeint(filterode,u0,t,args=(b,f)) #b is in tuple, needs comma
    return u[:,0]

fign = 0
dt = 0.025
t = array(arange(0,20,dt))
# linspace(0, 10, 400)
ft = filterresponse(t, fc)
bg = array([chop(tt, bglevel[1]+bglevel[0], bglevel[0]) for tt in t])

if model == 1:
    for n,ti in enumerate(t):
        ft[n] = chop(ti, -ft[n], 0)

signal1 = bg + ft
signal2 = bg
figure()
plot(t, signal1, label="Signal+bg")
plot(t, signal2, label="Pure bg")
xlabel("Time")
ylabel("Measurement")
legend(loc="best")
if savefigures:
    savefig("%s%d.%s"%(savefigbase,fign,savefigext))
    fign +=1

## Lockin calculation
def lockin(t, signal, ph, T):
    """ Generate lock-in amplifier output """
    def reference(t, *args):
        f, phi = args
        return sin(2*pi*f*t + phi)
    
    lockinsignal = reference(t, *(1,ph))*signal
    dt = t[1] - t[0]
    nT = int(T / dt)
    lockout = array([simps(lockinsignal[it-nT:it],t[it-nT:it]) for it in range(nT, len(t))])/T
    return (t[nT:], lockout)

# Example lockin signal
figure()
phex = 0/180*pi
lt, ls = lockin(t, signal1, phex, Ts)
plot(lt, ls, label="Signal+bg")
lt, ls = lockin(t, signal2, phex, Ts)
plot(lt, ls, label="Pure bg")
xlabel("Time")
ylabel("Lock-in signals")
legend(loc="best")
if savefigures:
    savefig("%s%d.%s"%(savefigbase,fign,savefigext))
    fign +=1


## Lockin phase dependence
phl = linspace(0, pi, 181)
ms = []
mb = []
for n,ph in enumerate(phl):
    # figure(n)
    lt, ls = lockin(t, signal1, ph, Ts)
    # plot(lt, ls, label='Signal' %(ph/pi*180))
    ms.append(mean(ls))
    lt, ls = lockin(t, signal2, ph, Ts)
    # plot(lt, ls, label='Pure background')
    mb.append(mean(ls))
    # title("Phi = %.1f" %(ph/pi*180))
    # legend(loc="best")
ms = array(ms)
mb = array(mb)

figure()
plot(ph2deg(phl), ms, label="Signal + bg")
plot(ph2deg(phl), mb, label="Pure bg")
xlabel("Lock-in phase (deg)")
ylabel("Lockin signal")
legend(loc="best")
if savefigures:
    savefig("%s%d.%s"%(savefigbase,fign,savefigext))
    fign +=1

ratio = abs(ms/mb)
pmax = argmax(ratio)
pmin = argmin(ratio)
ptrick = ph2deg((phl[pmax] + phl[pmin])/2)
figure()
plot(ph2deg(phl), abs(ms/mb),'.')
plot([ptrick, ptrick], [min(ratio), max(ratio)])
xlabel("Lock-in phase (deg)")
ylabel("(S+B) / B")
title("Resonance at %.1f deg" %(ptrick))
if savefigures:
    savefig("%s%d.%s"%(savefigbase,fign,savefigext))
    fign +=1


show()
