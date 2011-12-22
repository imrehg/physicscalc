"""
Library to do Zeeman slower specific calculations
"""
from __future__ import division
import numpy as np

#### Dimensioned
mU = 1.667e-27 # Mass: Atomic Mass Unit
h = 6.626e-34
hbar = h / 2 / np.pi
# uprimehbar = 1.399e10 * 2 * np.pi
bohrmag = 9.27400915e-24
uprimehbar = bohrmag/hbar

#### Atom classes

class Rb85:
    """ Rubidium-85 atom parameters """
    def __init__(self):
        self.G = 38.117e6 # 2pi x 6.07MHz
        self.k = 1281654.9389 * 2 * np.pi # 1/m
        self.m = 85 * mU  # Rb85
        self.aslow = hbar*self.k*self.G/(2*self.m)

class Rb87:
    """ Rubidium-87 atom parameters """
    def __init__(self):
        self.G = 38.117e6 # 2pi x 6.07MHz
        self.k = 1281654.9389 * 2 * np.pi # 1/m
        self.m = 86.909180 * mU  # Rb85
        self.aslow = hbar*self.k*self.G/(2*self.m)

class Cs133:
    """ Cesium-133 atom parameters """
    def __init__(self):
        self.G = 32.889e6 # 2pi x 5.234MHz
        self.k = 1173230.7104 * 2 * np.pi # 1/m
        self.m = 132.905451 * mU  # Cs133
        self.aslow = hbar*self.k*self.G/(2*self.m)

class K41:
    """ Potassium-41 atom parameters """
    def __init__(self):
        self.G = 37.919e6 # 2pi x 6.035MHz
        self.k = 13042.903375 * 2 * np.pi # 1/m (D2 line)
        self.m = 40.96182576 * mU  # K41
        self.aslow = hbar*self.k*self.G/(2*self.m)

#### Helper functions

def slowerlength(aslow, eta, v0, vf = 0):
    """
    Total length of Zeeman slower to slow atoms from a given initial
    velocity to a final one.

    aslow: constans decelartion
    eta: goodness factor (within [0, 1])
    v0: maximum capture velcoity
    vf: final velocity (optional, defaults to 0)
    """
    return (v0*v0 - vf*vf) / (2 * np.abs(aslow) * eta)


def bideal(atom, z, eta, v0, vf=0, detuning=0):
    """
    Calculate the ideal field of the zeeman slower

    atom: slowed atom, use the classes defined here, or something similar
    z: positions of query [m]
    v0: max captured velocity [m/s]
    vf: final velocity (defaults to 0) [m/s]
    detuning: laser detuning used [MHz], understood as red detuning

    Output:
    z < 0 : 0
    z > slowerlength : 0
    0 <= z <= slowerlength: calculated ideal slowing magnetic field 
    """
    # Can accept Numpy array, python list and single number, convert everything to numpy
    if type(z) != type(np.array([1])):
        if type(z) == type(1):
            z = np.array([z])
        elif type(z) == type([1]):
            z = np.array(z)
        else:
            raise TypeError('Input positions should be number, Python list or Numpy array')
    if not (0 < eta <= 1):
        raise ValueError('Efficiency has to be between 0 and 1')


    sl = slowerlength(atom.aslow, eta, v0, vf)
    detu = 2*np.pi*detuning*1e6 / uprimehbar

    # Filter z locations into three relevant regions
    z_pre = z<0
    z_post = z>sl
    z_mid = ~z_pre & ~z_post


    # Calculate field for the three regions
    b_pre = z[z_pre]*0
    b_post = z[z_post]*0
    b_mid = hbar*atom.k/bohrmag*np.sqrt(v0*v0 - 2 * eta * atom.aslow * z[z_mid]) - detu

    # Combine output
    res = np.append(b_pre, np.append(b_mid, b_post))
    return res

def bidealLs(atom, z, eta, sl, vf=0, detuning=0):
    """
    Calculate the ideal field of the zeeman slower

    atom: slowed atom, use the classes defined here, or something similar
    z: positions of query [m]
    v0: max captured velocity [m/s]
    vf: final velocity (defaults to 0) [m/s]
    detuning: laser detuning used [MHz], understood as red detuning

    Output:
    z < 0 : 0
    z > slowerlength : 0
    0 <= z <= slowerlength: calculated ideal slowing magnetic field
    """
    # Can accept Numpy array, python list and single number, convert everything to numpy
    if type(z) != type(np.array([1])):
        if type(z) == type(1):
            z = np.array([z])
        elif type(z) == type([1]):
            z = np.array(z)
        else:
            raise TypeError('Input positions should be number, Python list or Numpy array')
    if not (0 < eta <= 1):
        raise ValueError('Efficiency has to be between 0 and 1')

    detu = 2*np.pi*detuning*1e6 / uprimehbar

    v0 = np.sqrt(2 * sl * atom.aslow * eta + vf**2)

    # Filter z locations into three relevant regions
    z_pre = z<0
    z_post = z>sl
    z_mid = ~z_pre & ~z_post


    # Calculate field for the three regions
    b_pre = z[z_pre]*0
    b_post = z[z_post]*0
    b_mid = hbar*atom.k/bohrmag*np.sqrt(v0*v0 - 2 * eta * atom.aslow * z[z_mid]) - detu

    # Combine output
    res = np.append(b_pre, np.append(b_mid, b_post))
    return res


if __name__ == "__main__":
    import pylab as pl
    atom = Rb85()
    print slowerlength(atom.aslow, 0.7, 365, 20)
    zl = np.linspace(-0.2, 1, 40001)
    bz = bideal(atom, zl, 0.7, 365, 20, 252)
    pl.plot(zl, bz*1e4)
    pl.title('Ideal Zeeman-slower field strength (Bell2010)')
    pl.xlabel('position (m)')
    pl.ylabel('Magnetic field (G)')
    pl.show()
