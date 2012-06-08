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

    filename = "trajectory10"
    # fig = pl.figure(figsize=(fw, fh))

    # videal = (zs.uprimehbar*x+2*np.pi*detu*1e6)/k
    # pl.plot(ze, videal)

    tmax = 1.5/vf
    t = np.linspace(0, tmax, 1001)
    s = 4
    detuw = -2*np.pi*detu*1e6
    field = lambda z: lo.fieldcalc(z, setup)*I
    field = lambda z: zs.bideal(atom, z, eta, v0, vf, detuning=detu)
    b = (atom, field, s, detuw)

    vfield = (zs.uprimehbar*field(ze)-detuw)/k
    # pl.plot(ze, vfield, 'k-', linewidth=3, label='0 detuning velocity')

    out = []
    # for v in np.linspace(v0, vf, 1):
    vlist = np.linspace(v0, vf, 15)
    # vlist = [v0/2]
    z0 = -0.3
    # for v in vlist:
    #     print "Velocity: %.1f" %(v)
    #     start = time.time()
    #     y0 = [-0.3, v]
    #     y = odeint(eqmotion, y0, t, args=(b,))
    #     # pl.plot(y[:,0], y[:, 1], label="v=%.1f" %(v))
    #     pl.figure(1, figsize=(fw, fh))
    #     pl.plot(y[:,0], y[:, 1])

    #     pl.figure(2, figsize=(fw, fh))
    #     pl.plot(y[:,0], -(zs.uprimehbar*field(y[:, 0])/k - detuw/k - y[:, 1]))

    #     out += [y]
    #     print "Elapsed: %.1fs" %(time.time()-start) 

    # pl.figure(2, figsize=(fw, fh))
    # pl.xlim([z0, ze[-1]])
    # # pl.ylim([0, v0+20])

    # pl.figure(1, figsize=(fw, fh))
    # pl.xlim([z0, ze[-1]])
    # pl.ylim([0, v0+20])

    # np.savez(filename, vlist, out)
    # # pl.plot(ze, x)
    # # pl.plot(ze, y)
    # # pl.plot([0], [v0], 'x')
    # pl.xlabel('z position (m)')
    # pl.ylabel('velocity (m/s)');

    # # pl.plot(ze, bf)
    # # # pl.plot(ze[1:], np.diff(bf))
    # pl.legend(loc="best")

    # pl.savefig("%s.pdf" %(filename))
    # pl.savefig("%s.png" %(filename))


    vl = np.linspace(v0*0.9, v0*1.1, 501)
    a0 = -zs.hbar*k*atom.G/(2*atom.m)
    eta = 0.5
    z = [0]
    dd = k*a0*s/(vl*(1+s*4/atom.G**2*(k*vl + detuw-zs.uprimehbar*field(z))**2))
    dd2 = dd - a0*eta/np.sqrt(v0**2 + 2*a0*eta*z) * k
    pl.plot(vl/v0, dd)
    pl.plot(vl/v0, dd2)
    pl.plot(vl/v0, vl*0)
    pl.show()


    # tmax = 1.5/vf
    # t = np.linspace(0, tmax, 51)
    # z0 = -0.3
    # s = 10
    # detuw = -2*np.pi*detu*1e6
    # field = lambda z: lo.fieldcalc(z, setup)*I
    # b = (atom, field, s, detuw)

    # TASKS = [([v, z0],
    #           t,
    #           b) for v in range(250, 40, -50)]


    # try:
    #     import multiprocessing as processing
    # except:
    #     import processing
    # NUMBER_OF_PROCESSES = processing.cpu_count()
    # pool = processing.Pool(processes=NUMBER_OF_PROCESSES)  

    # try:
    #     out = pool.map_async(fly, TASKS).get(999999)
    # except KeyboardInterrupt:
    #     pool.terminate()
    #     sys.exit(0)

    # print out


    # for i in range(len(layer)):
    #     if np.diff(segments[i]) != 0:
    #         print "%7.2f mm to %7.2f mm: %3d turns / %3d layers" %((looppos[segments[i][0]]-wire[0]/2)*1000, (looppos[segments[i][1]-1]+wire[0]/2)*1000, np.diff(segments[i]), layer[i])
    # sta = looppos[segments[1][0]]-wire[0]/2
    # stb = looppos[segments[-2][1]-1]+wire[0]/2  # The upper boundary is the start of the next layer, need -1
    # efflen = zs.slowerlength(atom.aslow, eta, v0, vf) * 1000
    # coillen = (stb - sta)*1000
    # print "Effective length: %.1f mm" %(efflen)
    # print "Slower coils' length: %.1f mm (%.1f mm extra)" %(coillen, (coillen-efflen))

    # layers = lo.getLayerNumber(setup)

    # fig = pl.figure(figsize=(fh, fw))
    # fig.text(0.5, 0.95, '%s: slower = %g m, Max B = %g G' %(filename, sl, max(bfield)*1e4 ), horizontalalignment='center')
    # pl.subplot(211)
    # pl.plot(looppos, layers[0], 'k.')
    # pl.xlim([looppos[0], looppos[-1]])
    # pl.ylim([0, max(layer)+1])
    # pl.xlabel('Position (m)')
    # pl.ylabel('Layer number')

    # pl.subplot(212)
    # pl.plot(ze, lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup), 'k-', linewidth=3)
    # pl.plot(ze, ze*0, 'k--', linewidth=1)
    # pl.xlim([looppos[0], looppos[-1]])
    # pl.xlabel('Position (m)')
    # pl.ylabel('Normalized B field')

    # pl.savefig("%s.pdf" %(filename))
    # pl.savefig("%s.png" %(filename))
    # pl.show()
