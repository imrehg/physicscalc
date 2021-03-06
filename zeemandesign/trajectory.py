"""
Trajectory calculation for the zeeman slower
"""
import numpy as np
import pylab as pl

import layeroptimize as lo
import zeemanslower as zs
import sys

from scipy.integrate import odeint
import time

def eqmotion(y, t, b):
    atom, field, s, detu = b
    B = field([y[0]])[0]
    dy = y[1]
    ddy = -zs.hbar*atom.k*atom.G/(2*atom.m) * s / (1 + s + (4 / atom.G**2)*(detu + k * y[1] - zs.uprimehbar*B)**2)
    return [dy, ddy]


def fly(args):
    print "Starting:", y[1]
    y0, t, b  = args
    y = odeint(eqmotion, y0, t, args=(b,))
    return y

fw, fh = (11.69, 8.27)

if __name__ == "__main__":
    filename = sys.argv[1]
    atom = zs.Rb87()

    sim = np.load(filename)['simulation'][()]
    wire = sim['wire']
    setup = sim['setup']
    eta = sim['eta']
    v0 = sim['v0']
    vf = sim['vf']

    looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']

    print filename
    eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']
    sl = zs.slowerlength(atom.aslow, eta, v0, vf)
    z = np.array([0])
    bfield = zs.bideal(atom, z, eta, v0, vf, detu)
    z = np.append([-10*R], np.append(np.linspace(0, sl, 61), [sl+10*R]))
    ze = np.linspace(z[0], z[-1], 201)

    # print setup['csign']
    # bf = lo.fieldcalc(ze, setup)
    # # pl.plot(ze, bf)
    # pl.plot(ze[1:], np.diff(bf))

    series = 4
    ratio = 1

    setup['csign'] = [s*ratio if s < 0 else s for s in setup['csign']]
    # setup['csign'] = [s*ratio if s < 0 else s*ratio for s in setup['csign']]

    # Get the needed current to run the coils:
    bf = lo.fieldcalc(ze, setup)
    B0 = zs.bideal(atom, 0, eta, v0, vf, detu)
    I = B0[0] / lo.fieldcalc(0, setup)[0]
    print "Current:", I
    bf *= I
    

    k = atom.k
    # # bf = bf*0
    # # detu = 0
    # fulldetu = k*v - 2*np.pi*zs.uprimehbar*bf - detu*1e6*2*np.pi
    # pl.plot(ze, fulldetu/2/np.pi/1e6)
    # pl.show()
    
    # x = (k*v - 2*np.pi*zs.uprimehbar*bf - detu*1e6*2*np.pi)/2/np.pi/1e6
    # x = bf*1e4


    ub = zs.bohrmag/zs.h / 1e4
    # # x = k*v/2/np.pi/1e6 - bf*1e4 * ub / 1e6
    # # x = bf*1e4 *ub/1e6 - k*v/2/np.pi/1e6
    # x = 2*(zs.hbar*k/zs.bohrmag*np.sqrt(v0*v0 - 2*eta*atom.aslow*ze) - 2*np.pi*detu*1e6/zs.uprimehbar)
    x = zs.bideal(atom, ze, eta, v0, vf, detuning=detu)
    
    # # vx = (2*np.pi*zs.uprimehbar*bf/eta + detu*1e6*2*np.pi)/ k

    # y = 1/(-2*np.pi*detu*1e6 + k*v0 - zs.uprimehbar*x)

    filename = "trajectory03"
    fig = pl.figure(figsize=(fw, fh))

    # videal = (zs.uprimehbar*x+2*np.pi*detu*1e6)/k
    # pl.plot(ze, videal)

    tmax = 1.5/vf
    t = np.linspace(0, tmax, 1001)
    s = 2
    detuw = -2*np.pi*detu*1e6
    field = lambda z: lo.fieldcalc(z, setup)*I
    b = (atom, field, s, detuw)

    out = []
    # for v in np.linspace(v0, vf, 1):
    vlist = np.linspace(v0, vf, 15)
    vlist = np.linspace(v0, vf, 1)
    # vlist = np.linspace(1, 1, 1)
    # vlist = [v0/2]
    z0 = -0.3

    zv = np.linspace(z0, ze[-1], 1000)
    vfield = (zs.uprimehbar*field(zv)-detuw)/k
    pl.plot(zv, vfield, 'k-', linewidth=4, label='$\delta = 0$ velocity')


    for v in vlist:
        print "Velocity: %.1f" %(v)
        start = time.time()
        y0 = [-0.3, v]
        y = odeint(eqmotion, y0, t, args=(b,))
        # pl.plot(y[:,0], y[:, 1], label="v=%.1f" %(v))
        pl.plot(y[:,0], y[:, 1], linewidth=2)
        out += [y]
        print "Elapsed: %.1fs" %(time.time()-start) 

    pl.xlim([z0, ze[-1]])
    pl.ylim([0, v0+20])

    np.savez(filename, vlist, out)
    # pl.plot(ze, x)
    # pl.plot(ze, y)
    # pl.plot([0], [v0], 'x')
    pl.xlabel('z position (m)', fontsize=15)
    pl.ylabel('velocity (m/s)', fontsize=15);

    # pl.plot(ze, bf)
    # # pl.plot(ze[1:], np.diff(bf))
    pl.legend(loc="best")

    pl.title("Negative coil current scaling factor= %g, saturation parameter= %g" %(ratio, s), fontsize=15)
    pl.savefig("%s_%d.pdf" %(filename, series))
    pl.savefig("%s_%d.png" %(filename, series))
    pl.show()

