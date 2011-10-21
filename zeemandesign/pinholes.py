"""
Calculating and simulating the angular distribution after a collimator
in the Zeeman slower
"""
from __future__ import division
import numpy as np
import pylab as pl

def overlap(r1, r2, d):
    """ 
    Overlapping area of circles

    r1, r2: the two circle's radius
    d: distance of their centre
    """
    rmax = r1 if r1 >= r2 else r2
    rmin = r1 if r1 <= r2 else r2

    if d >= (r1+r2):
        ret = 0
    elif (d+rmin) > rmax:
        c1 = r1*r1*np.arccos( (d*d + r1*r1 - r2*r2) / (2 * d * r1) )
        c2 = r2*r2*np.arccos( (d*d + r2*r2 - r1*r1) / (2 * d * r2) )
        c3 = -0.5 * np.sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2) )
        ret = c1 + c2 + c3
    else:
        ret = np.pi*rmin*rmin
    return ret

def angleArea(th, r1, r2, l):
    ##### # d = np.sin(th)*l # Wrong
    # d = np.sin(th)*l
    d = np.tan(th)*l
    return overlap(r1, r2, d)

if __name__ == "__main__":


    # filename = "pinhole_1_20.npz"
    # filename = "pinhole_1_1.npz"
    filename = "pinhole_0.5_20.npz"
    if filename:
        data = np.load(filename)
        sim = data['sim'][()]
        r1 = sim['r1']
        r2 = sim['r2']
        l = sim['d']
        result = data['result'][()]
    else:
        r1 = 1
        r2 = 1
        l = 20
        result = None

    # th = np.linspace(0, np.pi/2, 1001)
    # area = []
    # for thval in th:
    #     area += [ angleArea(thval, r1, r2, l) ]
    # area = np.array(area)

    plotdouble = False
    if result is not None:
        toplot = [val[0] for val in result if val[2]]
        ntotal = len(result)
        ngood = len(toplot)
        print ngood/ntotal
        pl.figure(1, figsize=(11.69/1.5, 8.27/1.5))

        if plotdouble:
            total = [val[0] for val in result]
            ns, bins, patcheses = pl.hist([toplot, total], 30, alpha=0.75, label=['Simulation', 'Total atoms in angle bucket'])
            n = ns[0]
        else:
            n, bins, patches = pl.hist(toplot, 30, facecolor='green', alpha=0.75, label='Simulation')


        pos = bins[1:]-(bins[1]-bins[0])/2
        th = pos
        area = []
        for thval in th:
            area += [ angleArea(thval, r1, r2, l) ]
        area = np.array(area)

        scale = (np.sin(bins[1:]) - np.sin(bins[0:-1]))
        dist = area/(np.pi*r1*r1) * scale
        # dist2 = area/(np.pi*r1*r1) * scale

        # Still figure 1
        pl.plot(th, dist*ntotal, 'k-', linewidth=3, label='Theoretical prediction')
        pl.xlabel(r'Escape angle $\theta$')
        pl.ylabel('Number of atoms')
        pl.title('Total number: %d, r: (%g,%g), L: %g' %(ntotal, r1, r2, l))
        pl.legend(loc='best')
        # pl.savefig('pinhole_1_20.pdf')

        pl.figure(2)
        pl.plot(pos, n/ntotal,'--')
        pl.plot(th, dist, 'r-', linewidth=2)

        if plotdouble:
            pl.figure(3)
            pl.plot(pos, ns[0]/ns[1], '--')
            dist2 = area/(np.pi*r1*r1)
            pl.plot(th, dist2, 'k-', linewidth=3)
            pl.title('Escape probablity by angle')

        pl.show()
