from __future__ import division
from numpy import *
from scipy.odr import *
from scipy.integrate import quad
import pylab as pl
from scipy.interpolate import interp1d

# Gaussian beam
gaussbeam = lambda p, x: exp(-2*(x-p[1])**2/p[0]**2)

# ODR fit wrapper
def totalfit(x, y, f, p0):
    d = Data(x, y)
    model = Model(f)
    fit = ODR(d, model, p0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output

# The actual clever
def chopsignal(p, t):
    """ 
    Simulated chopped beam signal, approximation
    chopsignal(p, t)

    The signal is calculated as windowed integral of a gaussian beam, with time replacing
    the space.

        0  P/2  P  3P/2
        -----   -----
        |   |   |   |
        |   |   |   |
    -----   -----   ----
    Eg. at 0 delay, on the positive side we have an open window with width equal to half the
    chopping period.
    The signal calculated for shifting the window with a delay 0..P, and the output signal is
    interpolated by mapping the time onto the first shifted period: (time % period)
    It is not the ultimate accuracy but more than an order of magnitude quicker than direct simulation
    of the signal and seems to produce good enough results

    Input:
    p : simulation parameters
        [beam size, signal height, signal offset, chopping period, time offset]
    t : time steps
    Eg: chopsignal([330, 5e-4, 0, 1/500*1e6, 700], range(4000))

    Output:
    array of simulated values
    """
    # start = clock()
    gbp = [p[0], 0]
    A = p[1]
    offset = p[2]
    gapt = p[3]
    phase = p[4]
    # Simulate just one period of signal that we'll use for interpolation later
    delay = linspace(0, gapt, 151)
    beaminteg = lambda x, p: gaussbeam(gbp, x)
    chop = zeros(len(delay))
    # Integrate over a number of periods. Should be just 2 bright ones, it's bad sign if the beam is bigger than the chopping gap
    patt = range(-2, 2, 2)
    for i in xrange(len(delay)):
        for patti in patt:
            chop[i] += A*quad(beaminteg, patti*gapt/2+delay[i], (patti+1)*gapt/2+delay[i], args=(p))[0]
    chop += offset
    chopint = interp1d(delay, chop, kind='linear')
    return chopint((t-phase)%(gapt))

def chopfit(t, v, doplot=False, saveplot=False, filename="chopfit", fileext="png"):
    """
    Chopping signal fit, takes data that we normally have from oscilloscope
    """
    # Zero and rescale time
    t -= t[0]
    t *= 1e6

    vmax = max(v)
    vmin = min(v)
    vrange = vmax - vmin
    # Starting input parameters
    p0 = [330, 5e-4, 0, 1/502*1e6, 700]
    out = totalfit(t, v, chopsignal, p0)
    pout = out.beta
    period = out.beta[3]/1e6
    period_sd = out.sd_beta[3]/1e6
    beamsize = "w = %.2f (\pm %.2f) \mu s" %(out.beta[0], out.sd_beta[0])
    repetition = "\mathrm{rep} = %.2f (\pm %.2f) \mathrm{Hz} / %.0f \mu s" %(1/period, (1/period**2 *period_sd), period*1e6)
    print "beam size: %s\nrepetition rate: %s" %(beamsize, repetition)
    if doplot:
        pl.plot(t, v, 'k.', label='Data')
        pl.plot(t, chopsignal(pout, t), 'r-', linewidth=1.5, label="Fit")
        pl.xlabel("Time ($\mu s$)")
        pl.ylabel("Signal ($V$)")
        pl.title("$\mathrm{Fit:} %s; %s$" %(beamsize, repetition))
        pl.ylim([vmin - 0.05*vrange, vmax + 0.05*vrange])
        pl.xlim([t[0], t[-1]])
        pl.legend(loc='best')
        if saveplot:
            pl.savefig("%s.%s" %(filename, fileext))
        pl.show()


if __name__  == "__main__":
    filename = "ALL0001/F0001CH1.CSV"
    # filename = "ALL0002/F0002CH1.CSV"
    data = loadtxt(filename, delimiter=',', usecols=[3,4])
    t = data[:, 0]
    v = data[:, 1]
    plotname = (filename.split('/')[1]).split('.')[0]
    chopfit(t, v, doplot=True, saveplot=True, filename=plotname, fileext="pdf")
