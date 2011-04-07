#!/usr/bin/python2
"""
Simulation of Helmholtz coils to find optimal geometry for
vapour cell coil design.
"""
import multiprocessing as mp
import numpy as np
import pylab as pl
import scipy as sp
import scipy.integrate as spinteg

class Point3D(object):
    """ A positition in 3D space """

    def __init__(self, pos):
        self.pos = np.array(pos)

    def get(self):
        return self.pos

    def x(self):
        return self.pos[0]

    def y(self):
        return self.pos[1]

    def z(self):
        return self.pos[2]

    def __str__(self):
        return "Point: %e,%e,%e" %(self.pos[0], self.pos[1], self.pos[2])

def veclen(vec):
    return np.sqrt(np.dot(vec,vec))

def biotsavart(phi, r, coil, dim):
    """ Implement the integrand of the Biot Savart law"""
    dr = r.get() - coil.getrim(phi).get()
    dl = coil.gettangent(phi)
    bs = np.cross(dl, dr) / veclen(dr)**3
    return bs[dim]

    # R = coil.radius
    # dx = r.x() - (coil.centre.x() + R * np.cos(phi))
    # dy = r.y() - (coil.centre.y() + R * np.sin(phi))
    # dz = r.z() - coil.centre.z()
    # dr = (dx, dy, dz)
    # dl = (-R * np.sin(phi), R * np.cos(phi))
    # return dr



class Coil(object):
    """ Coil with normal vector along z """

    def __init__(self, centre, radius=1):
        self.centre = centre
        self.radius = radius

    def __str__(self):
        return "Coil centered at (%g, %g, %g) with radius %g" %(self.centre.x(), self.centre.y(), self.centre.z(), self.radius)

    def getrim(self, phi):
        return Point3D(self.centre.get() +
                       np.array([self.radius * np.cos(phi),
                                 self.radius * np.sin(phi),
                                 0]
                                )
                       )

    def gettangent(self, phi):
        """ Normalized tangent as a function of turn """
        vec = np.array([self.radius * -np.sin(phi),
                        self.radius * np.cos(phi),
                        0])
        vec = vec / veclen(vec)
        return vec

    def getfield(self, points):
        field = []
        for r in points:
            field += [np.array([spinteg.quad(biotsavart, 0, 2*np.pi, args=(r, self, i))[0] for i in range(3)])]
        return np.array(field)

def getval(args):
    return args[0].getfield(args[1])

class HelmholtzCoils(object):
    """ A Helmholtz coil pair """

    def __init__(self, centre1, centre2, width, turns):
        coilspaces = np.linspace(0, width, turns)
        self.coils = []
        for pos in coilspaces:
            coilc = Point3D(centre1.get() + [0, 0, pos])
            self.coils += [Coil(coilc)]

            coilc = Point3D(centre2.get() + [0, 0, -pos])
            self.coils += [Coil(coilc)]

    def getfield(self, points):
        pool = mp.Pool(4)
        result = pool.map(getval, [(coil, points) for coil in self.coils])
        return sum(result)

    def getfieldgrid(self, grid):
        # x, y, z = grid
        # test = []
        # for xi in x:
        #     for yi in y:
        #         for zi in z:
        #             test += [Point3D((xi, yi, zi))]

        result = []
        for coil in self.coils:
            result += [coil.getfield(grid)]
        print sum(result)
        # pool = mp.Pool(4)
        # result = pool.map(getval, [(coil, grid) for coil in self.coils])
        # return sum(result)
        

def zcircle(phi, R, z):
    """ A circle for integration """
    return Point3D((R*np.cos(phi), R*np.sin(phi), z))


if __name__ == "__main__":
    """ Command line interface """
    # centre = Point3D((0, 0, 0))
    # # coil = Coil(centre, 1)
    # # # print coil.getrim(np.pi/2)
    # # # print coil.gettangent(np.pi/2)
    # # # print biotsavart(0, centre, coil)
    # # # print coil.gettangent(np.pi/4)
    # zline = np.linspace(-3, 3, 31)
    # points = [Point3D((0, 0, z)) for z in zline]
    # # field = coil.getfield(points)
    # # pl.plot(zline, field[:,2])


    x = 3
    centre1 = Point3D((0, 0, x))
    centre2 = Point3D((0, 0, -x))
    hhcoil = HelmholtzCoils(centre1, centre2, 0.2, 2)
    # # hhcoil = HelmholtzCoils(coil1, coil2)
    # hfield = hhcoil.getfield(points)
    # pl.plot(zline, hfield[:,2])
    # pl.show()


    n = 10
    mz = 0.7
    my = 1.2
    y, z = np.linspace(-my, my, n), np.linspace(-mz, mz, n)
    # grid = []
    # for yi in y:
    #     for zi in z:
    #         grid += [Point3D((0, yi, zi))]
    # res = hfield = hhcoil.getfieldgrid(grid)
    # print res
    Y,Z = np.meshgrid( y, z)
    BY = np.zeros((len(y), len(z)))
    BZ = np.zeros((len(y), len(z)))

    points = []
    for i in xrange(len(y)):
        for j in xrange(len(z)):
            field = hhcoil.getfield([Point3D((0, Y[i,j], Z[i,j]))])
            BY[i, j] = field[0][1]
            BZ[i, j] = field[0][2]

    data = []
    for i in xrange(len(y)):
        for j in xrange(len(z)):
            data += [[Y[i, j], Z[i, j], BY[i, j], BZ[i, j]]]
    data = np.array(data)
    np.savetxt('test.csv', data)

    pl.quiver(Z, Y, BZ, BY)
    pl.show()

    # TODO:
    # Easy setting geometry (from command line?)
    # Save to text file
    # Create analytics separately
