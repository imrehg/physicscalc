from numpy import *
from scipy.integrate import odeint
from scipy.optimize import leastsq
from pylab import plot, show, semilogx

def getphase(fc, toplot):
    """ Calculate phase shift for low pass filter response to square wave
    fc: filter frequency, in chopped frequency units
    toplot: show calculation on plot? (True/False)

    returns phase shift in degrees
    """

    force = lambda t,f: mod(ceil(t*f*2),2)
    filterode = lambda u,t,b,f: -u[0]/b + force(t,f)

    t = arange(0,10,0.05)
    u0 = array([0,1])
    # fc = 0.5
    b = 1/(2*pi*fc)
    # f = 0.235
    f = 1.0
    u = odeint(filterode,u0,t,args=(b,f)) #b is in tuple, needs comma
    fit = lambda p, t: p[2]*sin(2*pi*p[0]*t+p[1])+p[3]
    errf= lambda p, t, x: fit(p, t) - x
    p0 = [f, 0.1, 2, 0]
    p1, success = leastsq(errf, p0, args=(t, u[:,0]/b))
    if (toplot):
        plot(t,u[:,0]/b)
        plot(t,force(t, f),label='original force')
        plot(t, fit(p1, t))
        show()
    return p1[1]/pi*180

print getphase(1,False)

# fcl = linspace(0.1,5,10)
# ps = [getphase(fc, False) for fc in fcl]
# semilogx(fcl, ps)
# show()
