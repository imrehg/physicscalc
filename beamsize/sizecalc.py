from numpy import *
from scipy import special
from scipy.odr import *
from scipy.integrate import quad
import pylab as pl

fitfunc = lambda p, x: sqrt(2*pi)*(special.erf((x-p[2])*sqrt(2)/p[0]) + 1)/(4/p[0]) * p[1] + p[3]
gauss = lambda p, x: p[1]/sqrt(2*pi*p[0]**2)*exp(-(x-p[2])**2/(p[0]**2))+p[3]
gaussbeam = lambda p, x: exp(-2*(x-p[2])**2/p[0]**2)
def totalfit(x, y, f, p0):
    d = Data(x, y)
    model = Model(f)
    fit = ODR(d, model, p0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output

filename = "ALL0001/F0001CH1.CSV"
# The signals, hand-set
l = [[0.0678, 0.0688]]
l = [[0.06795, 0.0687],
     [0.06895, 0.0697],
     [0.0699, 0.0707],
     [0.0709, 0.0717]]

data = loadtxt(filename, delimiter=",", usecols=[3,4])
t = data[:, 0]
v = data[:, 1]

w = []
ws = []
ct = []
n = 0
pl.figure()
pl.subplot(2,1,1)
for li in l:
    t1 = t[(li[0] <= t) & (t <= li[1])]
    v1 = v[(li[0] <= t) & (t <= li[1])]
    p0 = [0.0007, sign(t1[-1]-t1[0])*300, mean(t1), 0.2]
    out = totalfit(t1, v1, fitfunc, p0)    
    if n==0:
        pl.plot((t1-mean(t1))*1e6, v1, 'k.')
        wg = out.beta[0]
        cg = out.beta[2]
        # pl.plot(t1*1e6, gauss(pg, t1))
        pl.plot((t1-mean(t1))*1e6, fitfunc(out.beta, t1), 'r-', linewidth=1.5)
        pl.xlabel('Time ($\mu s$)')
        pl.ylabel('Signal (V)')

        # pl.plot(t1*1e6, gaussbeam(out.beta, t1)/6)
        # pl.ylim([0, 0.25])
    n += 1
    w.append(abs(out.beta[0]))
    ws.append(out.sd_beta[0])
    ct.append(out.beta[2])


w = array(w)
ws = array(ws)
wavg = sum(ws*w)/sum(ws)
wsavg = sqrt(sum(ws**2))
# Set the title of the above plot
pl.title('%s: w0 = %d us (@ %.2f Hz)' %(filename, wavg*1e6, 1/mean(diff(ct))/2))


#### Beam analysis results
print "Beam size: %.0f (+- %.0f) us (weighted average)" %(wavg*1e6, wsavg*1e6)
print "Gap frequency: %.2f" %(1/mean(diff(ct))/2)

### Values set by guesstimation, measure later
# Wheel radius (m)
R = 0.05
# Number of gaps
N = 30
# Distance between gaps (m)
d = 2 * pi * R / N
# Time between gaps (s)
td = mean(diff(ct))*2
# Speed of rotation (m / s)
vr = d / td
print "With parameters: wheel radius = %g m; number of gaps = %d" %(R, N)
print "=> Real beam size: %.2f (+- %.2f)mm" %(wavg * vr * 1e3, wsavg * vr * 1e3)

#### Create a simulated chopping signal with the above beam
pl.subplot(2,1,2)
p = [wavg, 0, 0, 0]
gapt = td/2
beaminteg = lambda x, p: gaussbeam(p, x)
delay = linspace(-td/2, 2.2*td, 201)
chop = zeros(len(delay))
patt = range(-8, 8, 2)
for i in xrange(len(delay)):
    for patti in patt:
        chop[i] += quad(beaminteg, patti*gapt+delay[i], (patti+1)*gapt+delay[i], args=(p))[0]
pl.plot((t-ct[0])*1e6-15, v/max(v), 'k.')
pl.plot(delay*1e6, chop/max(chop), 'r-', linewidth=1.5)
pl.xlabel("Time ($\mu s$)")
pl.ylabel("Chopping signal (a.u.)")
pl.ylim([-0.05, 1.05])
pl.savefig('fitting.pdf')
print "Fractional power: %.4f " %(quad(beaminteg, -gapt/2, gapt/2, args=(p))[0]/quad(beaminteg, -5*td, 5*td, args=(p))[0])
pl.show()
