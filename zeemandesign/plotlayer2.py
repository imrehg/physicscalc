import numpy as np
import numpy.random as random
import pylab as pl
import scipy.odr as odr
from sys import exit
import layeroptimize as layer

# import ourgui

# filename = ourgui.openFile('npz')
# if not filename:
#     exit(0)
# basefile = filename.split('/')[-1]
# # filename = "layers_110830_153517.npz"

filename = "15_AWG12Coat_9.npz"

data = np.load(filename)
print data
newsetup = data['setup'][()]
# zl, ideal = zip(*data['ideal'])
# print zl, ideal
zl = data['ideal'][0]
ideal = data['ideal'][1]
d0 = data['d0']

# print data['setup']
# print newsetup['layer']
# print data.files

zeropos = 31
zz = np.append(np.linspace(-0.3, 0.01, zeropos), np.linspace(0, 1.1, 201))
curr = layer.bideal(0)/layer.fieldcalc(zz, newsetup, d=d0)[zeropos]
print 'Current = %.2f A' %(curr)

wirelen = layer.totalcoillength(newsetup, d=d0)
print 'Wire length = %.3f m ' %(wirelen)

resist = 0.334 * wirelen / 100.0
power = resist * curr*curr
print 'Minimum power = %.2f W' %(power)


pl.figure(figsize=(11.69, 8.27), dpi=100)
pl.plot(zl, layer.normalize(ideal, 1), 'ko', markersize=6, label='fitted field')
pl.plot(zz, layer.normalize(layer.fieldcalc(zz, newsetup, d=d0), zeropos), 'r-', linewidth=2.4, label='coil field')

x = layer.getLayerNumber(newsetup)
pl.plot(newsetup['looppos'], np.array(x[0])/10.0, 'g--', label='coil wind number / 10')
pl.xlabel('Position')
pl.ylabel('Field / Coil wind number')
pl.ylim([-1.1, 1.1])
pl.title(r'pitch=%.1fmm, %.1fA, >%.1f$\Omega$, >%.1fW, %s' %(data['d0']*1000, curr, resist, power, basefile))
pl.legend(loc='best')


pl.savefig('%s.png' %basefile)
pl.savefig('%s.pdf' %basefile)

pl.show()
