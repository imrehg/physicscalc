import numpy as np
import scipy.odr as odr

# Built in pulse shapes
def gaussian(p, t):
    return p[1]/(p[0]*np.sqrt(2*np.pi))*np.exp(-(t)**2/(2*p[0]**2))

def delay(p, pulseshape, tau, w, chirp2, chirp3, ph0):
    t = np.linspace(-100, 100, 401)
    tprime = t - tau
    dt = t[1]-t[0]
    p1 = pulseshape(p, t)*np.exp(1j*(w*t + chirp2*t*t + chirp3*t**3))
    # p1b = pulseshape(p, tprime)*np.exp(1j*(w*tprime + chirp2*tprime*tprime + ph0))
    p1b = pulseshape(p, tprime)*np.exp(1j*(w*tprime + chirp2*tprime*tprime + chirp3*tprime**3 + ph0))
    return sum((abs((p1 + p1b)**2))**2*dt)

def envelope(p, pulseshape, tau):
    t = np.linspace(-100, 100, 401)
    tprime = t - tau
    dt = t[1]-t[0]
    env1 = pulseshape(p, t)
    env1p = pulseshape(p, tprime)
    return (sum((env1 + env1p)**4*dt), sum((env1 - env1p)**4*dt))

def envelopecurve(t, p, pulseshape):
    upper = np.array([])
    lower = np.array([])
    for ti in t:
        res = envelope(p, pulseshape, ti)
        upper = np.append(upper, res[0])
        lower = np.append(lower, res[1])
    return upper, lower

def autocorr(pulseparams, pulseshape, delays, carrier, chirp2, chirp3, ph0):
    tscale = 1e-15   # fs time scale
    c = 3e8
    wcarrier = 2*np.pi*c/(carrier*1e-9)*tscale

    out = np.array([])
    for tau in delays:
        out = np.append(out, delay(pulseparams, pulseshape, tau, wcarrier, chirp2, chirp3, ph0))
    return out

def fitfunc(beta, x, pulseshape):
    print beta
    pulsep = beta[0:2]
    carrier = beta[2]
    offset = beta[3]
    chirp2  = beta[4]
    chirp3 = beta[5]
    ph0 = beta[6]
    return autocorr(pulsep, pulseshape, x-offset, carrier, chirp2, chirp3, ph0)

def fitmodel(x, y, pin, pulseshape):
    data = odr.Data(x, y)
    model = odr.Model(fitfunc, extra_args=(pulseshape,))
    fit = odr.ODR(data, model, pin)
    fit.set_job(fit_type=2)
    out = fit.run()
    return out

# t = np.linspace(-100, 100, 10001)
# def delay(p, tau, w):
#     t = np.linspace(-100, 100, 10001)
#     return sum((abs((gaussian(p, t)*np.exp(1j*w*t) + gaussian(p, t-tau)*np.exp(1j*w*(t-tau)*(1+5e-5*tau)))**2))**2)
#     # return sum((abs(gaussian(p, t)*np.exp(1j*w*t) + gaussian(p, t-tau)*np.exp(1j*w*(t-tau))**2))**2)

# p = [0, 10, 1]
# print delay(p, 0, w)
# tlim = 40
# taulist = pl.linspace(-tlim, tlim, 401)
# auto = [delay(p, tau, w) for tau in taulist]
# # pl.plot(t, gaussian([0, 20, 1], t)*np.sin(w*t))
# pl.plot(taulist, auto)

# pl.show()
