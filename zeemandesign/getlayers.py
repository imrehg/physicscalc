"""
Print and display the results of the design optimization of the Zeeman slower layers
Numerical values and plotting.
"""
import numpy as np
import pylab as pl

import layeroptimize as lo
import zeemanslower as zs
import sys

# Figure width and height
fw, fh = (11.69, 8.27)  # A4

if __name__ == "__main__":

    # Choose filename for input
    try:
        filename = sys.argv[1]
    except IndexError:
        import ourgui
        filename = ourgui.openFile('npz')

    ## Choose atom
    atomq = -1
    while atomq not in range(3):
        try:
            atomq = int(raw_input("Atom type: (0) Rb85, (1) Rb87, (2) K41 ? "))
        except ValueError:
            continue
    if atomq == 0:
        atom = zs.Rb85()
    elif atomq == 1:
        atom = zs.Rb87()
    elif atomq == 2:
        atom = zs.K41()
    else:
        atom = zs.Rb87()  # default
    ##

    sim = np.load(filename)['simulation'][()]
    wire = sim['wire']
    setup = sim['setup']
    eta = sim['eta']
    v0 = sim['v0']
    vf = sim['vf']

    looppos, segments, layer = setup['looppos'], setup['segments'], setup['layer']

    eta, v0, vf, detu, R = sim['eta'], sim['v0'], sim['vf'], sim['detu'], sim['R']
    sl = zs.slowerlength(atom.aslow, eta, v0, vf)
    z = np.array([0])
    bfield = zs.bideal(atom, z, eta, v0, vf, detu)
    z = np.append([-5.5*R], np.append(np.linspace(0, sl, 61), [sl+5.5*R]))
    ze = np.linspace(z[0], z[-1], 201)

    ## Print winding data
    for i in range(len(layer)):
        if np.diff(segments[i]) != 0:
            print "%7.2f mm to %7.2f mm: %3d turns / %3d layers" %((looppos[segments[i][0]]-wire[0]/2)*1000, (looppos[segments[i][1]-1]+wire[0]/2)*1000, np.diff(segments[i]), layer[i])
    sta = looppos[segments[1][0]]-wire[0]/2
    stb = looppos[segments[-2][1]-1]+wire[0]/2  # The upper boundary is the start of the next layer, need -1
    efflen = zs.slowerlength(atom.aslow, eta, v0, vf) * 1000
    coillen = (stb - sta)*1000
    print "Effective length: %.1f mm" %(efflen)
    print "Slower coils' length: %.1f mm (%.1f mm extra)" %(coillen, (coillen-efflen))

    layers = lo.getLayerNumber(setup)

    fig = pl.figure(figsize=(fh, fw))
    fig.text(0.5, 0.95, '%s: slower = %g m, Max B = %g G' %(filename, sl, max(bfield)*1e4 ), horizontalalignment='center')
    ## Layers plot
    pl.subplot(211)
    pl.plot(looppos, layers[0], 'k.')
    pl.xlim([looppos[0], looppos[-1]])
    pl.ylim([0, max(layer)+1])
    pl.xlabel('Position (m)')
    pl.ylabel('Layer number')

    ## Field plot
    pl.subplot(212)
    pl.plot(ze, lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup), 'k-', linewidth=3)
    pl.plot(ze, ze*0, 'k--', linewidth=1)
    pl.xlim([looppos[0], looppos[-1]])
    pl.xlabel('Position (m)')
    pl.ylabel('Normalized B field')

    pl.savefig("%s.pdf" %(filename))
    pl.savefig("%s.png" %(filename))
    pl.show()
