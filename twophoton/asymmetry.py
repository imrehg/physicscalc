from __future__ import division
from numpy import *
import pylab as pl
from scipy.odr import *
import cProfile



n = 1001
s = 0.5
p = [0, 1, 1000]
lorentz = lambda p, x : p[2]*p[1]/((x - p[0])**2 + p[1]**2)+p[3]

def gendata(f, p, x, s):
    out = array([random.randn()*v*s+v for v in f(p, x)])
    # return f(p, x)+e
    return out

def dofit(f, x, y, p0):
    data = Data(x, y)
    model = Model(f)
    fit = ODR(data, model, p0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output

def twosides(f, x, y, p0):
    fit1 = dofit(f, x, y, p0)
    centre = fit1.beta[0]
    w = fit1.beta[1]
    xl = x[x <= centre]
    yl = y[x <= centre]
    xr = x[x >= centre]
    yr = y[x >= centre]
    fitl = dofit(f, xl, yl, fit1.beta)
    fitr = dofit(f, xr, yr, fit1.beta)
    wl = fitl.beta[1]
    wr = fitr.beta[1]
    sym = abs(wl-wr)/w
    # if sym > 0.5:
    #     print w, wl, wr
    #     pl.plot(x, f(fit1.beta, x))
    #     pl.plot(xl, yl, 'x')
    #     pl.plot(xl, f(fitl.beta, xl))
    #     pl.plot(xr, yr, '.')
    #     pl.plot(xr, f(fitr.beta, xr))
    return sym


def doanalyse(filename):
    # filename = "885_3991.csv"
    data = loadtxt(filename)
    ns = 1670
    ne = 1770
    y = data[ns:ne]

    x = array(range(len(y)))
    p = array([ 43.00373901,  11.79285268,  28.03482431, 1])


# xc = linspace(-3, 3, n)
# out = gendata(lorentz, p, xc, s)
# # pl.plot(xc, out-lorentz(p, xc), '.')
# # fitted = dofit(lorentz, xc, out, p)
# # pf = fitted.beta
# # pl.plot(xc, out, '.')
# # pl.plot(xc, lorentz(pf, xc))
# ww = twosides(lorentz, xc, out, p)
# wcrit = 0.0006

# print ww
# if (ww > wcrit):
#     print "Bigger"
# else:
#     print "Smaller"



# pl.figure(1)
# pl.plot(x, y, '.')
    out = dofit(lorentz, x, y, p)
    pp = out.beta
# print out.beta
# # pl.plot(x, lorentz(out.beta, x))
# # pl.plot(x, (y-lorentz(out.beta, x))/(y-out.beta[3]), '.')
# ey = (y-lorentz(out.beta, x))/(y-out.beta[3])
# print std(ey)
# pl.plot(abs(ey), '.')

# pl.plot(x, y, '.')

    s = 0.16
    ww = twosides(lorentz, x, y, pp)
    rep = 10000
    # wlist = array([])
    # for r in xrange(rep):
    #     wlist = append(wlist, twosides(lorentz, x, gendata(lorentz, pp, x, s), pp))
    #     # if wlist[-1] > 0.5:
    #     #     pl.plot(x, y, 'o')
    #     #     break

    # print ww
    # print len(wlist[wlist > ww])/rep
    # pl.hist(wlist, 30)

### Compare residuals
    new = gendata(lorentz, pp, x, s)
    out2 = dofit(lorentz, x, y, p)
    pp2 = out2.beta
    pl.plot(x, new-lorentz(pp2, x), 'x')
    pl.plot(x, y-lorentz(pp, x),'o')

# pl.show()

filename = "885_3991.csv"
doanalyse(filename)
# cProfile.run('doanalyse(filename)')
pl.show()
