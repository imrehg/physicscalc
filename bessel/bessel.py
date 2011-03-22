import numpy as np
import pylab as pl
import scipy as sp
from scipy.special import jv

maxorder = 3
moddepth = np.linspace(0, 10, 100001)

pl.figure(num=None, figsize=(8.27, 11.69), dpi=120, facecolor='w', edgecolor='k')

### Comparison of orders
pl.subplot(2,1,1)

markers = ['-',':','--','-.']
name = ['0th', '1st', '2nd', '3rd', '4th']

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
        text += "Orders: %d -> %d\n%s\n" %(o1, o2, crossv)
        pl.plot(crossv, olist[o1][zval], '.', markersize=10)

pl.text(1.5, 0.35, text, size=12)
pl.xlabel('Modulation depth')
pl.ylabel('Amplitude (normalized)')
pl.legend(loc="upper right")

### Total
pl.subplot(2,1,2)
for i in xrange(len(olist)):
    olist[i] *= (-1)**i
total = sum(olist)**2
tdiff = np.diff(total)
mins = sp.where(tdiff[:-1] * tdiff[1:] < 0)
mstep = moddepth[1] - moddepth[0]
pl.plot(moddepth, total)
pl.plot(moddepth[mins][0::2]+mstep/2, np.zeros(len(mins[0][0::2])), '.', markersize=10)
text = "Min: %s\nMax: %s\nMaxval: %s" %(moddepth[mins][0::2], moddepth[mins][1::2], total[mins][1::2])
pl.text(1.2, 0.3, text, size=12)
pl.xlabel('Modulation depth')
pl.ylabel('Predicted interraction')
pl.title('Up to order %d' %(maxorder))
print "== Total effect: extremas =="
print moddepth[mins]+mstep/2

pl.show()
