import numpy as np
import pylab as pl
import scipy as sp
from scipy.special import jv

maxorder = 2
moddepth = np.linspace(0, 10, 100001)

markers = ['-',':','--']
name = ['0th', '1st', '2nd']

olist = []
for order in xrange(0, maxorder+1):
    thisorder = jv(order, moddepth)**2
    olist += [thisorder]
    pl.plot(moddepth, thisorder, markers[order], label=name[order])

text = ''
for o1 in xrange(0, maxorder):
    for o2 in xrange(o1+1, maxorder+1):
        bdiff = olist[o1] - olist[o2]
        zval = sp.where(bdiff[:-1] * bdiff[1:] < 0)
        crossv = moddepth[zval]
        print "=== Equal amplitude for orders %d - %d ===" %(o1, o2)
        print crossv
        text += "Orders: %d -> %d\n%s\n\n" %(o1, o2, crossv)
        pl.plot(crossv, olist[o1][zval], '.', markersize=10)

pl.text(3, 0.3, text, size=12)
pl.xlabel('Modulation depth')
pl.ylabel('Amplitude (normalized)')
pl.legend(loc="upper right")
pl.show()
