import numpy as np
import pylab as pl

import wires
import zeemanslower as zs
from layeroptimize import *

def dowire(wire, n1, n2, series=0, maxtry=100):
    R = 0.0383
    eta, v0, vf, detu = 0.7, 365, 20, 260
    # The field that we want to match
    nz = 61
    atom = zs.Rb85()
    sl = zs.slowerlength(atom.aslow, eta, v0, vf)
    z = np.append([-5.5*R], np.append(np.linspace(0, sl, nz-2), sl+5.5*R))
    bfield = zs.bideal(atom, z, eta, v0, vf, detu)

    # # Our coil
    nloops, loops, layer, csign, segments = createStruct((z[0],z[-1]), wire[0], (n1, n2))

    setup = {'looppos': loops,
             'layer': layer,
             'csign': csign,
             'segments': segments,
             'R': R,
             'd': wire[0],
             }

    newsetup = optimize(z, bfield, setup, maxtry=maxtry)
    # newsetup = setup
    # pl.plot(z, bfield, 'x')
    ze = np.linspace(z[0], z[-1], 201)
    pl.figure(figsize=(11.69, 8.27))
    pl.plot(ze, fieldcalc(ze, newsetup)/fieldcalc(0, newsetup), 'r-', linewidth=2, label='coil field')
    pl.plot(z, bfield/bfield[1], 'ko', markersize=5, label='target field')
    pl.title("%s, R: %g mm, v: %d-%d m/s, %d MHz, (%d, %d)" %(wire[2], R*1e3, v0, vf, detu, n1, n2))
    pl.xlabel('position (m)')
    pl.ylabel('normalized magnetic field')
    pl.legend(loc='best')
    # pl.plot(zz, normalize(fieldcalc(zz, newsetup), 10), 'r-')
    # print newsetup['segments']

    simulation = {'R': R,
                  'eta': eta,
                  'v0': v0,
                  'vf': vf,
                  'detu': detu,
                  'wire': wire,
                  'setup': newsetup
                  }

    savefile = "%d_%s" %(series, wire[2])
    np.savez('%s' %(savefile), simulation=simulation)
    pl.savefig('%s.png' %(savefile))
    # pl.show()
    
if __name__ == "__main__":
    maxtry = 20000
    for wire in wires.AWG[0::1]:
        for n in range(5, 11):
            starts = n
            print wire[2]
            dowire(wire, n, n, starts, maxtry)

    pl.show()
