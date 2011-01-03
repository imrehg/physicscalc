#!/usr/bin/env python
"""
Effect of atomic dispersion on the observed line-shape in a cavity.

Using some values from the experiment and trying to understand why would
or wouldn't one see asymmetric or otherwise changed lineshapes.
Ultimate aim is to fit a fluorescence data and match it to this theory.

Greg, Mon Jan  3 17:24:28 CST 2011
"""
import numpy as np
import pylab as pl

def lorentz(p, x):
    return p[2]/((x - p[0])**2 + p[1]**2)

def dispersion(p, x):
    return -2*p[2]*(x-p[0])/((x - p[0])**2 + p[1]**2)**2 + 1

def cavityfreq(p, natom, Lair):
    """
    Parameters:
    constans params:
    p = [Mode_number, n_air, cell_length]
    natom = atomic index of refraction
    Lair = beam path in air
    """
    N = p[0]
    nair = p[1]
    Lcell = p[2]
    return N/(2*(nair * Lair  + natom * Lcell))

def findf(param, L0, dL):
    n = 1
    cavityparam = param[3:6]
    f0 = cavityfreq(cavityparam, 1, L0)
    for i in xrange(20):
        f = cavityfreq(cavityparam, n, L0+dL) - f0
        n = dispersion(param, f)
    return (f, n)

### Dimensionless variables:
## Scaling: Gamma = 1, c = 1 
v0 = 299792458
t0 = 1/(2*np.pi*2e6)
d0 = v0*t0
f0 = 1/(t0*2*np.pi)

## From the experiment
Lair = 0.03 / d0
Lcell = 0.04 / d0
nair = 1
N = 170000

# p = [N, nair, Lcell]
# fc = cavityfreq(p, 1, Lair)
# dl = 2e-9 /d0
# df = fc-cavityfreq(p, 1, Lair+dl)
# print df*f0/1e6

# p = [0, 1, 1]
# xlist = np.linspace(-4, 4, 101)
# y = lorentz(p, xlist)
# pl.plot(xlist, y)


nodisp = [0, 1, 0] + [N, nair, Lcell]
disp = [0, 1, 1e-8] + [N, nair, Lcell]
fcent, ncent = findf(nodisp, Lair, 0)

dllist = np.linspace(-4, 4, 101)*1e-9 / d0
fdisp = np.zeros(len(dllist))
fnodisp = np.zeros(len(dllist))
for i, dl in enumerate(dllist):
    f1, n1 = findf(nodisp, Lair, dl)
    fnodisp[i] = f1
    f2, n2 = findf(disp, Lair, dl)
    fdisp[i] = f2

pl.figure()
pl.plot(fnodisp, (dispersion(disp, fnodisp)-1)*1e8)
pl.xlabel('Laser frequency [$\Gamma$]')
pl.ylabel(r'Refractive index $(n_{atom} - 1) \times 10^8$')

pl.figure()
pl.plot(fnodisp, fnodisp, label='No dispersion')
pl.plot(fnodisp, fdisp, label='Dispersion')
pl.legend(loc="best")
pl.xlabel('Not dispersed resonance frequency [$\Gamma$]')
pl.ylabel('Resonance frequency [$\Gamma$]')
pl.title('Resonance frequency when tuning cavity length')

# pl.plot(fnodisp, fdisp-fnodisp)
plor = [0, 1, 1]
pl.figure()
pl.plot(fnodisp, lorentz(plor, fnodisp), label='No dispersion')
pl.plot(fnodisp, lorentz(plor, fdisp), label='Dispersion')
pl.xlabel(r'Not dispersed frequency [$\Gamma$]')
pl.ylabel(r'Fluorescence [a.u.]')
pl.legend(loc="best")
pl.show()
