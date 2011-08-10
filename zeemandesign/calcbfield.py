from __future__ import division
import numpy as np
import scipy.integrate as integ



#### From another source: ##########
## Constants:
mu = 1.26e-6*1.5e-3
def bzfield(p, params, z):
    cord = params['cord']
    cp = params['cp']
    cn = params['cn']
    cpnew = [cp[i]*(cord-i) for i in range(0,cord)]
    cnnew = [cn[i]*(cord-i) for i in range(0,cord)]
    R = params['R']
    I = params['I']
    upperp = np.polyval(cpnew, p)
    lowerp = (R**2 + (np.polyval(cp[cord+1:], p)-z)**2)**(3/2)
    uppern = np.polyval(cnnew, p)
    lowern = (R**2 + (np.polyval(cn[cord+1:], p)-z)**2)**(3/2)
    return mu*I/(4*np.pi)*(upperp/lowerp + uppern/lowern)

## Coil parameters
cp = [0, 1.58e-3, 0, -0.302, 8.698, 1.999, 0, -0.106, 0.581]
cn = [-0.0012, 0.0035, 0.02, -0.23, -1.257, -2, 3.14, 0.043, 0.62]
cord = 6
## Input parameters
R = 0.0383
I = 1
params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
nump = 501
z = np.linspace(-0.5, 1.5, nump)
res = np.array([integ.quad(bzfield, 0, 2*np.pi, args=(params, zz))[0] for zz in z])

# ## End result: interpolated field for the slower
# bfinter = interp1d(z, res)
# ############### 

np.savetxt('field.csv', zip(z, res))
import pylab as pl
pl.plot(z, res)
pl.show()
