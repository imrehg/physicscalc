#!/usr/bin/python2
"""
Simulation of Helmholtz coils to find optimal geometry for
vapour cell coil design.
"""
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

class Solenoid(object):
    """ Solenoid with fixed pitch """

    def __init__(self, pitch, windnum):
        self.pitch = pitch
        self.windnum = windnum
        self.philim = np.array([-self.windnum/2.0, self.windnum/2.0])*2*np.pi

    def showsolenoid(self):
        p = self.windnum * 200
        phi = np.linspace(-self.windnum/2.0, self.windnum/2.0, p)*2*np.pi
        x = np.cos(phi)
        y = np.sin(phi)
        z = phi / (2*np.pi) * self.pitch
        pl.plot(z, y)
        pl.show()

    def getrim(self, phi):
        return np.array([np.cos(phi),
                         np.sin(phi),
                         phi / (2*np.pi) * self.pitch,
                         ])

    def gettangent(self, phi):
        """ Normalized tangent as a function of turn """
        vec = np.array([-np.sin(phi),
                         np.cos(phi),
                         1.0 / (2*np.pi) * self.pitch,
                         ])
        return vec / veclen(vec)

    def getfield(self, points, dims=[0,1,2]):
        return [spinteg.quad(biotsavart, self.philim[0], self.philim[1], args=(points, self, i))[0] for i in dims]



if __name__ == "__main__":
    # """ Command line interface """

    print "Solenoid  parameters"
    pitch = float(raw_input("Pitch (dz/turn): "))
    windnum = int(raw_input("Wind number: "))
    dz = float(raw_input("Z displacement: "))
    nump = int(raw_input("Number of points / axis: "))
    print
    print "Starting..."

    sole = Solenoid(pitch, windnum)

    ypoints = np.linspace(-1, 1, nump)
    zpoints = np.linspace(-1, 1, nump)
    phi = np.linspace(0, np.pi*2, 100)

    X, Y = np.meshgrid(ypoints, zpoints)
    BX = np.zeros((nump, nump))
    BY = np.zeros((nump, nump))
    BZ = np.zeros((nump, nump))
    for i in xrange(nump):
        for j in xrange(nump):
            if (X[i, j]**2 + Y[i, j]**2) < 1:
                bx, by, bz = sole.getfield(np.array([X[i, j], Y[i, j], dz]))
                BX[i, j] = bx
                BY[i, j] = by
                BZ[i, j] = bz

    filename = "solenoid_%s" %(strftime("%y%m%d_%H%M%S"))
    np.savez(filename, X=X, Y=Y, BX=BX, BY=BY, BZ=BZ, pitch=pitch, windnum=windnum, dz=dz)

    print "Done / Saved"
    import soleanalyze
    soleanalyze.runanalysis()
