""" Get the (ideal) field data from the simulation file """
import numpy as np
import pylab as pl

import layeroptimize as lo
import zeemanslower as zs

fw, fh = (11.69, 8.27)

if __name__ == "__main__":
    filename = "9_AWG12_10.npz"
    # filename = "9_AWG18_10.npz"
    # filename = "9_AWG18_15.npz"
    atom = zs.K41()

    # filename = "3_AWG10_8.npz"
    # atom = zs.Rb85()

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

    # z = np.append([-5.5*R], np.append(np.linspace(0, sl, 61), [sl+5.5*R]))
    z = np.append([-15.5*R], np.append(np.linspace(0, sl, 61), [sl+17.5*R]))
    ze = np.linspace(z[0], z[-1], 201)

    for i in range(len(layer)):
        if np.diff(segments[i]) != 0:
            print "%7.2f mm to %7.2f mm: %3d turns / %3d layers" %((looppos[segments[i][0]]-wire[0]/2)*1000, (looppos[segments[i][1]-1]+wire[0]/2)*1000, np.diff(segments[i]), layer[i])
    sta = looppos[segments[1][0]]-wire[0]/2
    stb = looppos[segments[-2][1]]+wire[0]/2
    efflen = zs.slowerlength(atom.aslow, eta, v0, vf) * 1000
    coillen = (stb - sta)*1000
    print "Effective length: %.1f mm" %(efflen)
    print "Slower coils' length: %.1f mm (%.1f mm extra)" %(coillen, (coillen-efflen))

    layers = lo.getLayerNumber(setup)

    fig = pl.figure(figsize=(fh, fw))
    fig.text(0.5, 0.95, '%s: slower = %g m, Max B = %g G' %(filename, sl, max(bfield)*1e4 ), horizontalalignment='center')
    pl.subplot(211)
    pl.plot(looppos, layers[0], 'k.')
    pl.xlim([looppos[0], looppos[-1]])
    pl.ylim([0, max(layer)+1])
    pl.xlabel('Position (m)')
    pl.ylabel('Layer number')

    pl.subplot(212)
    pl.plot(ze, lo.fieldcalc(ze, setup)/lo.fieldcalc(0, setup), 'k-', linewidth=3)
    pl.plot(ze, ze*0, 'k--', linewidth=1)
    pl.xlim([looppos[0], looppos[-1]])
    pl.xlabel('Position (m)')
    pl.ylabel('Normalized B field')

    outdata = zip(ze, lo.fieldcalc(ze, setup))
    outfile = "".join([filename, ".txt"])
    np.savetxt(outfile, outdata)

    # pl.savefig("%s.pdf" %(filename))
    # pl.savefig("%s.png" %(filename))
    pl.show()
