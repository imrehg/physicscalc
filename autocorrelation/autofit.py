import numpy as np
import pylab as pl
# import autocorr as ac

def gaussian(p, t):
    return p[2]/(p[1]*np.sqrt(2*np.pi))*np.exp(-(t-p[0])**2/(2*p[1]**2))

def middle(p, t):
    return gaussian(p[0:3], t) + p[-1]

def inter(p, t):
    return (np.sin(p[3]*time + p[4] + p[5]*time**2 + p[6]*time**3)**4 - 0.5)*gaussian(p[0:3], t)*3

def totalfit(p, t):
    return middle(p, t)+inter(p, t)


if __name__ == "__main__":
    import scipy.odr as odr

    filename = 'F0004CH1.CSV'
    data = np.loadtxt(filename, delimiter=',', skiprows=18, usecols=(3,4))
    delayscale = 4.73   # fs/ms, from notes
    time = data[:, 0] * 1e3 * delayscale
    sig = data[:, 1]

    # Guessed parameters
    p2 = [0.5, 20, 60, 2*np.pi*0.2, 1, 0.0002, 0.00002, 0.5]

    data = odr.Data(time, sig)
    model = odr.Model(totalfit)
    fit = odr.ODR(data, model, p2)
    fit.set_job(fit_type=2)
    out = fit.run()
    out.pprint()
    pout = out.beta

    # Plots
    pl.subplot(2, 1, 1)
    pl.plot(time, sig, '.')
    pl.plot(time, totalfit(pout, time))
    pl.xlim([time[0], time[-1]])
    pl.xlabel('Delay time (fs)')
    pl.ylabel('Autocorrelatior signal')
    pl.title('Fit with quartic chirp')

    pl.subplot(2, 1, 2)
    pl.plot(time, sig - totalfit(pout, time))
    pl.xlim([time[0], time[-1]])
    pl.xlabel('Delay time (fs)')
    pl.ylabel('Residuals')
    pl.show()
