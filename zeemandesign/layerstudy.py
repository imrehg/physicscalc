import numpy as np
import pylab as pl

import wires
import zeemanslower as zs
from layeroptimize import *

def dowire(inparam):
    params, wire, nlayer, simparam = inparam
    R = params['R']
    eta = params['eta']
    v0 = params['v0']
    vf = params['vf']
    sl = params['Ls']
    detu = params['detu']
    atom = params['atom']
    n1, n2 = nlayer, nlayer
    maxtry = simparam['maxtry']
    series = simparam['series']
    printprogress = simparam['printprogress']

    # The field that we want to match
    nz = 61
    z = np.append([-5.5*R], np.append(np.linspace(0, sl, nz-2), sl+5.5*R))
    bfield = zs.bidealLs(atom, z, eta, sl, vf, detu)

    # # Our coil
    nloops, loops, layer, csign, segments = createStruct((z[0],z[-1]), wire[0], (n1, n2))

    setup = {'looppos': loops,
             'layer': layer,
             'csign': csign,
             'segments': segments,
             'R': R,
             'd': wire[0],
             }
    print("Optimize: %s, %d layer" %(wire[2], nlayer))
    newsetup = optimize(z, bfield, setup, maxtry=maxtry, printprogress=printprogress)
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

    savefile = "%d_%s_%d" %(series, wire[2], nlayer)
    np.savez('%s' %(savefile), simulation=simulation)
    pl.savefig('%s.png' %(savefile))
    
if __name__ == "__main__":

    import multiprocessing as processing
    import itertools

    # # design parameters
    # atom = zs.Rb85()
    # R = 0.03 / 2 # larger diameter slower tube
    # eta = 0.5 # efficiency
    # Ls = 0.6 # set slower length
    # vf = 30
    # detu = 160
    # series = 3
    # ## Derived parameters
    # v0 = np.sqrt(2 * Ls * atom.aslow * eta + vf**2)
    # print("Max capture velocity: %g" %(v0))

    # design parameters
    atom = zs.K41()
    R = 2.75 * 2.54e-2 / 2 # 2.75" slower tube
    eta = 0.7 # efficiency
    Ls = 0.6 # set slower length
    vf = 30
    v0 = 476
    detu = 328
    series = 5
    ## Derived parameters
    Ls = zs.slowerlength(atom.aslow, eta, v0, vf)
    print("Slower length: %g" %(Ls))

    # Sim parameters
    maxtry = 30000

    ## Derived parameters
    v0 = np.sqrt(2 * Ls * atom.aslow * eta + vf**2)
    print("Max capture velocity: %g" %(v0))

    input_params = {'R': R,
                    'eta': eta,
                    'Ls': Ls,
                    'v0': v0,
                    'vf': vf,
                    'detu': detu,
                    'atom': atom
                    }
    simparam = {'maxtry': maxtry,
                'series': series,
                'printprogress': False,
                }

    wirelist = wires.AWG
    maxlayers = range(5, 11)
    TASKS = [(input_params, x[0], x[1], simparam) for x in itertools.product(wirelist, maxlayers)]

    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)

    try:
        out = pool.map_async(dowire, TASKS).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)

    pl.show()
