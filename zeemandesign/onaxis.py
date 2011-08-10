"""
Zeeman slower design

Basic on-axis calculation.
"""
from __future__ import division
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import scipy.integrate as integ
from scipy.interpolate import interp1d
import pylab as pl

##### Units
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
tU = 1/GRb # Time: 1/linewidth
xU = 1/kRb # Distance: 1/wavenumber
# Planck constant rescaled
hU = 6.626e-34 / mU / xU / xU * tU
hbar = hU / 2 / np.pi
# Bohr-magneton rescaled
uB = hU * 1.399e6 * tU
BU = 1/uB # Magnetic field unit

# Dimensionless numbers
mRb = 86.909


# Some variables
k = 1.0
G = 1.0
s0 = 10
det0 = -2
Bb = 1 / BU
Bt = 3 / BU
gP3 = 1.3362
gS1 = 2.002
me = -4
mg = -3
uprime = uB * (gP3*me - gS1*mg)

def bfield(z, z0, Bb, Bt):
    """ Currently non-functional """
    return z * 0
    # return Bb - Bt * np.sqrt(1 - z / z0)

def decelerate(y, t):
    """ Decelartion function """
    x, v = y[0], y[1]
    return [v, -hbar / (2*mRb) * s0 / (1 + s0 + 4*(det0 + v + bfield(x, Ls, Bb, Bt)*uprime/hbar)**2)]

# #### From another source: ##########
# ## Constants:
# mu = 1.26e-6*1.5e-3
# def bzfield(p, params, z):
#     cord = params['cord']
#     cp = params['cp']
#     cn = params['cn']
#     cpnew = [cp[i]*(cord-i) for i in range(0,cord)]
#     cnnew = [cn[i]*(cord-i) for i in range(0,cord)]
#     R = params['R']
#     I = params['I']
#     upperp = np.polyval(cpnew, p)
#     lowerp = (R**2 + (np.polyval(cp[cord+1:], p)-z)**2)**(3/2)
#     uppern = np.polyval(cnnew, p)
#     lowern = (R**2 + (np.polyval(cn[cord+1:], p)-z)**2)**(3/2)
#     return mu*I/(4*np.pi)*(upperp/lowerp + uppern/lowern)

# ## Coil parameters
# cp = [0, 1.58e-3, 0, -0.302, 8.698, 1.999, 0, -0.106, 0.581]
# cn = [-0.0012, 0.0035, 0.02, -0.23, -1.257, -2, 3.14, 0.043, 0.62]
# cord = 6
# ## Input parameters
# R = 0.0383
# I = 110
# params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
# nump = 51
# z = np.linspace(-0.3, 1.2, nump)
# res = np.array([integ.quad(bzfield, 0, 2*np.pi, args=(params, zz))[0] for zz in z])*1e4 / BU  # magnetic field, T->G, then dimensionless 

# ## End result: interpolated field for the slower
# bfinter = interp1d(z, res)
# ############### 

# detu = -260e6 * 2*np.pi * tU
# v = 365 / xU * tU
# print detu, v

# # print bfinter(0)*uprime/hbar
# z = np.linspace(-0.3, 1.2, 401)
# field = (bfinter(z) * uprime / hbar) / tU
# pl.plot(z, field)
# pl.plot(0, (detu+k*v)/tU, 'x')
# pl.ylabel('detuning [1/s]')
# pl.show()

b = 0.015 * 1e4 * uprime / hbar
print b



# Ls = 1 / xU
# v0 = 1 / xU * tU

# t = sp.linspace(0, 2e6, 100000)

# # # Btlist = xrange(0, 100, 1000)
# # # Btlist = [(0, 0)]
# # for B in Btlist:
# #     Bb, Bt = B
# #     y0 = [0, 100*v0]
# #     print y0
# #     y = odeint(decelerate, y0, t)

# #     # pl.subplot(211)
# #     # pl.plot(t, y[:, 0]/Ls)
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Position')
# #     # pl.ylim([0, 1])

# #     # pl.subplot(212)
# #     # pl.plot(t, y[:, 1])
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Velocity')

# #     # pl.plot(t, det0 + y[:, 1])
# #     pl.plot(t, det0 + y[:, 1] + bfield(y[:, 0], Ls, Bb, Bt)*uprime/hbar)
# #     # pl.plot(y[:, 0], y[:, 0]/Ls)
# #     # pl.xlabel('Time')
# #     # pl.ylabel('Velocity')

# vlist = range(10, 200, 10)
# vfinal = []
# for v in vlist:
#     y0 = [0, v*v0]
#     y = odeint(decelerate, y0, t)

#     xi = y[:, 0] / Ls
#     vi = y[:, 1]

#     if xi[-1] < 1:
#         vfinal += [0]
#     else:
#         vfinal += [vi[xi <= 1][-1]]

#     # # filter negative direction
#     # xi = xi[vi > 0]
#     # t = t[vi > 0]
#     # vi = vi[vi > 0]

#     # # filter slower length
#     # t = t[xi < 1]
#     # vi = vi[xi < 1]
#     # xi = xi[xi < 1]

#     # filter negative direction
#     # xi = xi[vi > 0]
#     # t = t[vi > 0]
#     # vi = vi[vi > 0]
#     # pl.plot(t[vi>0], xi[vi>0])

#     # pl.subplot(211)
#     # # pl.plot(t*tU, xi)
#     # # pl.plot(xi)
#     # pl.xlabel('Time (s)')
#     # pl.ylabel('Position (1 = slower length)')
#     # # pl.legend(loc='best')
#     # # pl.ylim([0, 1])

#     # pl.subplot(212)
#     # pl.plot(t*tU, vi)
#     # pl.xlabel('Time (s)')
#     # pl.ylabel('Velocity')
#     pl.plot(xi[vi>0], det0 + vi[vi>0])
#     pl.xlim([0, 1])

# # pl.plot(vlist, vfinal)
# pl.show()


