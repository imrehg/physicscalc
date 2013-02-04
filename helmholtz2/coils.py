#!/usr/bin/env python2
"""
Calculating the magnetic field by the Anti-Helmholtz coils
"""
import numpy
import pylab as pl

u0 = 4e-7 * numpy.pi * 1e4  # Vacuum permeability (G*m/A)

def bz(R, z):
    """ Magnetic field on axis of a current loop, for unit current
    R: loop radius (m)
    z: distance from the centre along the axis (m)
    """
    return u0 * R**2 / (2 * (z**2 + R**2)**1.5)

def loop(params, pos, simtitle=""):
    """ Loop calculation for wires:
    params:
    R0 : inner radius
    L0 : inner distance to centre
    d : wire diameter
    n : radial wind number
    m : axial wind number
   
    pos: positions to simulate (should be around the centre, that is L0+-delta)
    simtitle: Simulation title for plot
    """
    R0, L0, d, n, m = params

    # Calculate the field created by all wires
    field = numpy.zeros(len(pos))
    for i in range(n):
        for j in range(m):
            rij = R0 + d*(i + 0.5)  # positive offset going out
            lij = pos - d*(j + 0.5)  # negative offset going away from the centre
            field += bz(rij, lij)

    # Anti-Helmholtz configuration
    field2 = field[::-1]
    totalfield = field-field2

    # Gradient fit
    p = numpy.polyfit(pos-L0, totalfield, 1)
    gradcm = p[0] / 100  # Gradient G/cm

    # plot
    pl.plot((pos-L0)*100, totalfield, 'ko', label="Calculated")
    pl.plot((pos-L0)*100, numpy.polyval(p, pos-L0), 'r-', label="Linear fit", linewidth=2)
    pl.xlabel("Distance from centre (cm)")
    pl.ylabel("Magnetic field (G)")
    pl.legend(loc="best")
    pl.title(simtitle+" gradient: %.2f (G/cm)/A" %(gradcm))
    pl.grid()

    return gradcm

if __name__ == "__main__":

    # General parameters
    R0 = 55e-3  # inner radius
    L0 = 29.26e-3  # centre offset

    # Build settings
    d = 4e-3
    n, m = 3, 6
    simtitle = "Hollow core"

    # Calculation
    params = (R0, L0, d, n, m)
    pos = numpy.linspace(L0-0.01, L0+0.01, 41)
    gradcm = loop(params, pos, simtitle)
    print "Gradient: %.2f (G/cm)/A" %(gradcm)
    pl.show()
