from __future__ import division
from numpy import *
import pylab as pl
from scipy.odr import *
import cProfile
import twophotonfit
from time import time
from scipy.interpolate import UnivariateSpline

n = 1001
s = 0.5
p = [0, 1, 1000]
lorentz = lambda p, x : p[2]*p[1]/((x - p[0])**2 + p[1]**2)+p[3]
gauss = lambda p, x: p[2]/p[1]*exp(-(x-p[0])**2/(2*p[1]**2))+p[3
]
def gendata(f, p, x, s):
    out = array([random.randn()*v*s+v for v in f(p, x)])
    # return f(p, x)+e
    return out

def gendata2(f, p, x, sf):
    out = array([random.randn()*sf(v)+v for v in f(p, x)])
    return out


def twosides(f, x, y, p0):
    fit1 = dofit(f, x, y, p0)
    centre = fit1.beta[0]
    xl = x[x <= centre]
    yl = y[x <= centre]
    xr = x[x >= centre]
    yr = y[x >= centre]
    fitl = dofit(f, xl, yl, fit1.beta)
    fitr = dofit(f, xr, yr, fit1.beta)
    # if sym > 0.5:
    #     print w, wl, wr
    xc = linspace(x[0], x[-1], 10001)
    def getwidth(x, fit):
        yc = abs(f(fit.beta, x)-fit.beta[1]/2)
        mm = argmin(yc)
        return abs(x[mm]-fit.beta[0])
    w = getwidth(xc, fit1)
    wl = getwidth(xc, fitl)
    wr = getwidth(xc, fitr)
    sym = abs(wl-wr)/w
    # print w, wl, wr
    # pl.plot(x, f(fit1.beta, x))
    # pl.plot(xl, yl, 'x')
    # pl.plot(xl, f(fitl.beta, xl))
    # pl.plot(xr, yr, '.')
    # pl.plot(xr, f(fitr.beta, xr))
    # pl.show()
    return sym



def dofit(f, x, y, p0):
    data = Data(x, y)
    model = Model(f)
    fit = ODR(data, model, p0)
    fit.set_job(fit_type=2)
    output = fit.run()
    return output
    #     pl.plot(xr, yr, '.')
    #     pl.plot(xr, f(fitr.beta, xr))
    return sym


def doanalyse(filename):
    # filename = "885_3991.csv"
    data = loadtxt(filename, delimiter=',')

    x = (data[:,0]-1.7e8)/1e6
    y = -data[:,1]

    xc = linspace(x[0], x[-1], 10001)
    ###############  First fit
    # p = [-0.36400664,  1.2765653,  0.4023463,  0.0770604]
    # fitted = dofit(lorentz, x, y, p)
    # # fitted.pprint()
    # pf = fitted.beta
    # # print twosides(lorentz, x, y, p)
    # pl.figure()
    # pl.plot(x, y, '.')
    # pl.plot(xc, lorentz(pf, xc), '-')
    # pl.title('lorentzian')
    # pl.figure()
    # pl.plot(x, y-lorentz(pf, x), '.')
    # pl.title('lorentzian')

    # p = [-0.36400664,  1.2765653,  0.4023463,  0.0770604]
    # fitted = dofit(gauss, x, y, p)
    # # fitted.pprint()
    # pf = fitted.beta
    # # print twosides(lorentz, x, y, p)
    # pl.figure()
    # pl.plot(x, y, '.')
    # pl.plot(xc, gauss(pf, xc), '-')
    # pl.title('gauss')
    # pl.figure()
    # pl.plot(x, y-gauss(pf, x), '.')
    # pl.title('gauss')


    # p = [-0.36400664,  0.22765653,  0.04023463,  0.04, 0.00770604]
    # fit = twophotonfit.dofit(twophotonfit.multipledouble_freeg, x, y, p)
    # fit.pprint()
    # yy = twophotonfit.multipledouble_freeg(fit.beta, x)

    p = [-3.6400664,  0.22765653,  0.4023463,  0.1, 0.4, 0.00770604]
    fit = twophotonfit.dofit(twophotonfit.multipletripple_freeg, x, y, p)
    # fit.pprint() 
   # yy2 = twophotonfit.multipletripple_freeg(fit.beta, x)
    # pl.plot(x, y, '.')

    # for i in range(len(fit.beta)):
    #     print "%d : %f (+- %f)" %(i, fit.beta[i], fit.sd_beta[i])

    # # pl.plot(x, y-yy2, '.')
    # pl.plot(x, yy2, '-')
    # pl.show()

    #######################

    # pf = [-3.64759824,  0.29662501,  0.60663077, -0.86707082,  0.27488604,  0.01065825]
    pf = fit.beta
    print pf
    f = twophotonfit.multipletripple_freeg
    fity = f(pf, x)
    u = len(x) - len(pf) - 1

    ee = y - fity
    sddd = std(ee)
    n = 61
    a = linspace(min(fity), max(fity), n+1)
    ey = array([])
    em = array([])
    # for i in xrange(n):
    #      eey = ee[(a[i] <= fity) & (fity <= a[i+1])]
    #      print len(eey)
    #      ey = append(ey, std(eey))
    #      em = append(em, mean(eey))
    # h = a[1]-a[0]
    # pl.plot(fity, abs(ee), '.')
    # pl.plot(a[:-1]+h/2, ey, '-')
    # # pl.plot(a[:-1]+h/2, em, '-')
    # # # aa = a[:-1]+h/2
    # aa = a[:-1]+h/2
    # z = polyfit(aa, ey, 1)
    # pz = poly1d(z)
    # pl.plot(aa, pz(aa))
    # print ey
    # # # pl.plot(x, y, '.')
    # # # pl.plot(x, fity)
    # # pl.plot(y, abs(ee), '.')
    # pl.show()


    ordery = zip(fity, ee)
    sy, se = zip(*(sorted(ordery)))
    sy = array(sy[-7000:])
    se = array(se[-7000:])
    ey = array([])
    xy = array([])
    nn = 70
    for i in range(int(len(sy)/nn)):
        ey = append(ey, std(se[(i*nn):((i+1)*nn)]))
        xy = append(xy, mean(sy[(i*nn):((i+1)*nn)]))
    pl.plot(sy, abs(se), '.')
    pl.plot(xy, ey)
    iny = linspace(sy[0], sy[-1], 1001)
    extrap = UnivariateSpline( xy, ey, k=5 )
    pl.plot(iny, extrap(iny), '-')
    pl.show()
    # chired = sum((y - fity)**2/pz(fity)**2)/u
    # print chired

    # pl.figure()
    # pl.plot(x, y, '.')
    # pl.plot(xc, f(pf, xc), '-', linewidth=2)
    # pl.xlim([xc[0], xc[-1]])
    # pl.xlabel('detuning (MHz)')
    # pl.ylabel('Signal (a.u.)')
    # pl.savefig('fit.png')

    # pl.figure()
    # pl.plot(x, y-fity, '.')
    # pl.xlim([xc[0], xc[-1]])
    # pl.xlabel('detuning (MHz)')
    # pl.ylabel('residuals')
    # pl.savefig('residuals.png')

    

    # pl.show()

    # ##### Scan systematics
    # pl.plot(fity, ee, '-')
    # pl.show()

    ### Error variance stuff
    # print x
    # yy = gendata2(f, pf, x, pz)
    # print yy
    # pl.plot(y, y-fity, '.')
    # # pl.plot(x, yy-fity,'x')
    # pl.show()

    # ssq = ee**2/pz(fity)**2
    # ssq = ee**2/sddd**2
    # pl.plot(x, ee**2/pz(fity)**2)
    # pl.plot(x, ee**2)
    # pl.show()

    # print len(x)
    # print sum(ssq)

    # ##########################3
    # ww = twosides(f, x, y, pf)
    # print ww

    ### Artificial skew
    # yy = gendata2(f, pf, x, pz)
    # sx = 0.06
    # yy = (yy - pf[-1])*(1 + sx*(x-pf[0]))+pf[-1]
    # ww2 = twosides(f, x, yy, pf)
    # # pl.plot(x, y, '.')
    # # pl.plot(x, yy, 'x')
    # pl.plot(x, yy-f(pf, x), 'x')
    # print ww, ww2
    # pl.show()



    ## Monte carlo repetition
    # # rep = 1000
    # # wlist = array([])
    # # start = time()
    # # for r in xrange(rep):
    # #     # yyy = gendata2(f, pf, x, pz)
    # #     # pl.plot(x, yyy-fity, 'x')
    # #     # pl.plot(x, y-fity, '.')
    # #     # pl.show()
    # #     # break
    # #     wlist = append(wlist, twosides(f, x, gendata2(f, pf, x, pz), pf))
    # #     print r, wlist[-1]
    # # #     # if wlist[-1] > 0.5:
    # # #     #     pl.plot(x, y, 'o')
    # # #     #     break
    # # elapsed = time() - start
    # # print "Original: ", ww
    # # print len(wlist[wlist > ww])/rep
    # # print "Time: %f (%f / rep)" %(elapsed, elapsed / rep)
    # # savetxt("montecarlo.txt", wlist)
    # # pl.hist(wlist)
    # # pl.show()
    


filename = "822_lineshape.csv"
doanalyse(filename)

# hist = loadtxt("montecarlo.txt")
# pl.hist(hist)
# pl.show()
