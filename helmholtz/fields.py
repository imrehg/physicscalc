#!/usr/bin/python2
"""
Simulation of Helmholtz coils to find optimal geometry for
vapour cell coil design.
"""
# import multiprocessing as mp
import numpy as np
import pylab as pl
import scipy.integrate as spinteg
from time import strftime

def veclen(vec):
    return np.sqrt(np.dot(vec,vec))

def biotsavart(phi, r, coil, dim):
    """ Implement the integrand of the Biot Savart law"""
    dr = r - coil.getrim(phi)
    dl = coil.gettangent(phi)
    bs = np.cross(dl, dr) / veclen(dr)**3
    return bs[dim]

class Coil(object):
    """ Coil with normal vector along z """

    def __init__(self, centrez):
        self.centrez = centrez

    def getrim(self, phi):
        return np.array([np.cos(phi),
                         np.sin(phi),
                         self.centrez]
                        )

    def gettangent(self, phi):
        """ Normalized tangent as a function of turn """
        vec = np.array([-np.sin(phi),
                        np.cos(phi),
                        0])
        return vec

    def getfield(self, points, dims=[1,2]):
        """ points are [(x, y, z)] format """
        field = []
        for r in points:
            field += [np.array([spinteg.quad(biotsavart, 0, 2*np.pi, args=(r, self, i))[0] for i in dims])]
        return np.array(field)

def getval(args):
    return args[0].getfield(args[1])

class HelmholtzCoils(object):
    """ A Helmholtz coil pair """

    def __init__(self, centrez, width, turns):
        coilspaces = np.linspace(-width, width, turns) + centrez
        self.coils = []
        for pos in coilspaces:
            self.coils += [Coil(pos)]

    def getfield(self, yparam, zparam):
        maxy, numy = yparam
        ypoints = np.linspace(0, maxy, numy)
        maxz, numz = zparam
        zpoints = np.linspace(0, maxz, numz)
        Y, Z = np.meshgrid(ypoints, zpoints)
        BY = np.zeros((numz, numy))
        BZ = np.zeros((numz, numy))
        # Careful with the indeces!
        for i in xrange(numy):
            for j in xrange(numz):
                for coil in self.coils:
                    points = [(0, Y[j, i], Z[j, i]),
                              (0, Y[j, i], -Z[j, i])]
                    fields = coil.getfield(points, dims=[1,2])
                    # BY[i, j] += fields[0][0] - fields[1][0]
                    # BZ[i, j] += fields[0][1] + fields[1][1]
                    BY[j, i] += fields[0][0] - fields[1][0]
                    BZ[j, i] += fields[0][1] + fields[1][1]
        return Y, Z, BY, BZ

if __name__ == "__main__":
    # """ Command line interface """

    print "Coil parameters"
    print " => Distances are scaled such that coil radius R=1"
    coilpos = float(raw_input("Coil position: "))
    coilwidth = float(raw_input("Coil width: "))
    nturns = float(raw_input("Number of turns "))
    print "Analysis parameters"
    maxy = float(raw_input("Maximum Y: "))
    numy = int((raw_input("Number of Y points: ")))
    maxz = float(raw_input("Maximum Z: "))
    numz = int((raw_input("Number of Z points: ")))
    
    

    print
    print "Starting..."
    hhcoil = HelmholtzCoils(coilpos, coilwidth/2, nturns)
    Y, Z, BY, BZ =  hhcoil.getfield((maxy, numy), (maxz, numz))

    filename = "helmholtz_%s" %(strftime("%y%m%d_%H%M%S"))

    np.savez(filename, Y=Y, Z=Z, BY=BY, BZ=BZ, coilpos=coilpos, coilwidth=coilwidth, nturns=nturns)

    print "Done / Saved"
    import analyze
    analyze.runanalysis()
