"""
Optimizing the Zeeman slower arrangement and save results
Set parameters inside the software for the slower
"""
import numpy as np
import pylab as pl
import sys

import wires
import zeemanslower as zs
from layeroptimize import *

def dowire(inparam):
    """ Run the wire simulation

    Input parameters:
    inparam = (params, wire, nlayer, simparam)

    * params is a dictionary with keys:
    'R': inner radius in m
    'eta': efficincy
    'v0': final velocity in m/s
    'vf': max capture velocity in m/s
    'sl': slower length in m
    'detu': red detuning in MHz (meaning e.g. -100MHz detuning fom resonance is detu=100)
    'atom': type of atom slowed, as in zeemanslower.py

    * wire: type of wire (just like in wires.py)

    * nlayer: max number of layers

    * simparam is a dictonary with keys:
    'maxtry': integer, maximum rounds of simulation, the total "time" of simulated annealing
    'series': integer, numerical counter for the simulation, used for the output name
    'printprogress': boolean, whether to show simulation rounds as they are done
    """
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
    ze = np.linspace(z[0], z[-1], 201)

    pl.figure(figsize=(11.69, 8.27))
    pl.plot(ze, fieldcalc(ze, newsetup)/fieldcalc(0, newsetup), 'r-', linewidth=2, label='coil field')
    pl.plot(z, bfield/bfield[1], 'ko', markersize=5, label='target field')
    pl.title("%s, R: %g mm, v: %d-%d m/s, %d MHz, (%d, %d)" %(wire[2], R*1e3, v0, vf, detu, n1, n2))
    pl.xlabel('position (m)')
    pl.ylabel('normalized magnetic field')
    pl.legend(loc='best')

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

    # design parameters
    atom = zs.Rb87()  # chosen atom
    R = 0.019 / 2  # inner radius
    eta = 0.5  # efficiency
    Ls = 0.59  # set slower length m
    vf = 30  # final velocity m/s
    detu = 175  # red detuning MHz
    series = 20  # name
    ## Derived parameters
    v0 = np.sqrt(2 * Ls * atom.aslow * eta + vf**2)
    print("Max capture velocity: %g" %(v0))

    # Sim parameters
    maxtry = 30000

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

    ## Previously used options
    # wirelist = wires.AWG
    # wirelist = wirelist[0:15]

    # Use enamel coated AWG12 wire
    wirelist = [(2.1e-3, 5.211e-3, "AWG12Coat")]

    maxlayers = range(5, 16)
    TASKS = [(input_params, x[0], x[1], simparam) for x in itertools.product(wirelist, maxlayers)]

    NUMBER_OF_PROCESSES = processing.cpu_count()
    pool = processing.Pool(processes=NUMBER_OF_PROCESSES)

    try:
        out = pool.map_async(dowire, TASKS).get(999999)
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(0)
