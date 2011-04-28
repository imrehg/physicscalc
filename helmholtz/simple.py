#!/usr/bin/env python
"""
Back to basics Helmholtz coil calculation to verify some
of the results from the more complex calculations
"""
import numpy as np
import pylab as pl

def axisfield(rad, xpos):
    """ Simple field calculation 

    rad: coil radius
    xpos: position along the axis (coil centered on 0)
    """
    return rad**2/(2*(rad**2 + xpos**2)**(1.5))

def helmholtz(pos, width, radius, coils, xpos):
    """ Helmholtz coil field calculation on the axis """

    # Position of single coils, both sides
    coilposright = np.linspace(pos-width/2.0, pos+width/2.0, coils)
    coilpos = np.append(coilposright, -coilposright)

    test = lambda a: sum(axisfield(radius,
                                   np.array([a-zpos for zpos in coilpos])))
    return [test(x) for x in xpos]

def scanradius(pos, width, radius, coils, xpos):
    """ Doing simulation for a range of different radii coils """
    scanfunc = lambda Rthis: helmholtz(pos, width, Rthis, coils, xpos)
    return [scanfunc(r) for r in radius]

def singlecoil(radius, offset, xpos):
    """ Single coil pair calculation
    
    radius : coil pair radius (can be array)
    offset : distance of coils from centre
    xpos : position to check Bz(x)/Bz(0)
    """
    return (axisfield(radius, (xpos - offset)) + \
                axisfield(radius, (xpos + offset)))/ \
        (axisfield(radius, -offset) + axisfield(radius, offset))


if __name__ == "__main__":
    pl.figure(figsize=(8.27, 11.69), dpi=100)

    # Simple setting
    RADII = np.logspace(-2, 1, 201)
    CPOS = 1
    MAXZ = 0.25
    RES = singlecoil(RADII, CPOS, MAXZ)
    pl.subplot(2, 1, 1)
    pl.semilogx(RADII, RES, linewidth=2)
    pl.xlabel('Coil-pair radius')
    pl.ylabel("Ratio Bz(%g)/Bz(0)" %(MAXZ))
    pl.title("Single coil pair at z=+-%g, field checked at z=%g" %(CPOS, MAXZ))
    pl.plot([RADII[0], RADII[-1]], [1, 1])

    # Detailed calculation
    pl.subplot(2, 1, 2)
    POS = np.linspace(0, 0.5, 31)
    COILPOS = 1
    WIDTH = 1
    NCOIL = 40
    RADII = np.linspace(1, 2, 11)
    SCANRAD = scanradius(COILPOS, WIDTH, RADII, NCOIL, POS)
    for i, scan in enumerate(SCANRAD):
        pl.plot(POS, scan/scan[0], label="R=%g" %(RADII[i]))
    pl.xlabel('z')
    pl.ylabel('Bz(z) / Bz(0)')
    pl.title('Helmholtz coil: pos=+-%g, width=%g, turns=%d'
             %(COILPOS, WIDTH, NCOIL))
    pl.legend(loc='best')
    pl.show()
