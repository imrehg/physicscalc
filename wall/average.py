from numpy import *
from scipy.integrate import quad

def xfunc(th, x0):
    """
    Distance of a point on the unit circle (defined by th angle) from
    a point on the axis, [x0, 0]
    """
    return sqrt(sin(th)**2 + (cos(th) - x0)**2)

def avgint(f, x0):
    """
    Average function f(x) over the unit circle, at the point
    shifted to [x0, 0]
    """
    integrand = lambda th, x0 : f(xfunc(th, x0))
    return quad(integrand, 0, 2*pi, x0)[0]/2/pi


if __name__ == "__main__":
    ## Example, for van der Waals type of potential
    import pylab as pl

    xtest = lambda x: x**(-3)
    x0l = linspace(0, 0.9)

    res = array([])
    for x0 in x0l:
        res = append(res, avgint(xtest, x0))
    pl.plot(x0l, res)

    pl.xlabel("x0, test position distance from centre")
    pl.ylabel("Averaged function value over the circle")
    pl.show()

