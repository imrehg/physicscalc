from __future__ import division
import numpy as np
import pylab as pl
import scipy.odr as odr
import scipy.integrate as integrate

def gauss(p, x):
    return p[2]*np.exp(-(x - p[0])**2/(2*p[1]**2))


def backg(p, x):
    return p[1]*x+p[0]

def noisy(p, x):
    return gauss(p[2:], x)+backg(p[0:2], x)

def runfit(x, y, func, p0):
    data = odr.Data(x, y)
    model = odr.Model(func)
    fit = odr.ODR(data, model, p0)
    fit.set_job(fit_type=2)
    return fit.run()

def integc(t, v, v0, allan):
    s2 = allan(t)**2
    return np.real(np.exp(-2 * np.pi**2 * t**2 * s2) * np.exp(2*np.pi*1j*(v-v0)*t))

def lorentz(p, x):
    return p[2]*(p[1]/2)/np.pi/((x-p[0])**2 + (p[1]/2)**2)


filename = 'allan_110720_191644.log'
# filename = 'allan_110720_185103.log'
comment="Greg's laser"
gate, allan, avgf, minf, maxf = np.loadtxt(filename,
                                           delimiter=',',
                                           comments='#',
                                           unpack=True,
                                           )
gate = np.log10(gate)
allan = np.log10(allan)
p0 = [3, -1/4, -3, 1, 2]
fit = noisy(p0, gate)

fitfunc = noisy
out = runfit(gate, allan, fitfunc, p0)
xlim = [2, 2]
xhat = np.linspace(gate[0]-xlim[0], gate[-1]+xlim[1], 301)
# xhat = np.linspace(gate[0], gate[-1], 301)

out.pprint()
fitparam = out.beta
fitparam[4] = 0
# fitparam[1] = -0.25
# fitparam[0] -= 0.5

pl.figure(1, figsize=(11.69, 8.27))
pl.suptitle("%s : %s" %(filename, comment))
pl.subplot(221)
pl.loglog(10**gate, 10**allan,'.')
pl.loglog(10**xhat, 10**fitfunc(fitparam, xhat), 'r-')
pl.xlabel('Gating time (s)')
pl.ylabel('Allan deviation (Hz)')

# pl.plot(gate, allan-backg(out.beta, gate),'.')
# pl.plot(xhat, fitfunc(out.beta, xhat)-backg(out.beta, xhat), 'r-')


allanf = lambda x : 10**noisy(fitparam, np.log10(x))
t = 10**xhat
dv = np.linspace(-2e5, 2e5, 501)
ilim = [t[0], 5e-4]
v0 = np.mean(avgf)

pl.figure(1)

pl.subplot(222)
pl.semilogx(t, integc(t, v0, v0, allanf))
pl.xlabel("Gating time (s)")
pl.ylabel("Integrand")

pl.subplot(224)



out = []
for d in dv:
    out += [2*integrate.trapz(integc(t, v0+d, v0, allanf), t)]
out = np.array(out)
# Normalize by centre
out /= out[int(len(out)/2)+1]
pl.plot(dv/1e3, out, 'r.')
pl.xlabel('Frequency detuning from carrier (kHz)')
pl.ylabel('Intensity (A.U.)')

# lp0 = [0, 1e4, 0]
# lfit = runfit(dv, out, lorentz, lp0)
# lfit.pprint()
# pl.plot(dv, lorentz(lfit.beta, dv), 'b-')

### Find halfwidth
halfd = dv[dv>0]
halfv = np.abs(out[dv>0] - 0.5)
halfwidth = halfd[np.argmin(halfv)]/1e3
pl.plot([-halfwidth, halfwidth], [0.5, 0.5], 'k-')
pl.text(dv[0]/1e3/2, 1.05, "Halfwidth ~ %.3f kHz" %(halfwidth))

pl.show()
