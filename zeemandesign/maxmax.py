import numpy as np
import flux
import scipy.integrate as integ
import pylab as pl

#### Dimensioned
GRb = 38.117e6 # 2pi x 6.07MHz
kRb = 1281654.9389 * 2 * np.pi # 1/m
mU = 1.667e-27 # Mass: Atomic Mass Unit
# mRb = 86.909 * mU
mRb = 85 * mU  # Rb85
h = 6.626e-34
hbar = h / 2 / np.pi
uprimehbar = 1.399e10 * 2 * np.pi
uBr = 9.27e-24
kB = 1.3806488e-23


T = 353

u = np.sqrt(2 * kB * T / mRb)
vdist3 = lambda v: v**3 * np.exp(-v*v/u/u)
vdist2 = lambda v: v**2 * np.exp(-v*v/u/u)
# total = integ.quad(vdist, 0, maxv)[0]

# print total
v  = np.linspace(0, 800, 100)

v3 = vdist3(v)
s3 = max(v3)
v3 /= s3
v2 = vdist2(v)
v2 /= max(v2)

vmax = 253
vlist = np.linspace(0, vmax, 50)

pl.plot(v, v2, 'r--')
pl.plot(v, v3, 'k-')
pl.fill_between(vlist, vdist3(vlist)/s3)
pl.show()
