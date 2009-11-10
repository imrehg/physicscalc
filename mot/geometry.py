import scipy as sp
import pylab as pl

a = 30
th = sp.linspace(20,45,100) / 180 * sp.pi
rmax = a / 2 * sp.sin(th)
thdeg = th / sp.pi * 180
grad = (rmax[-1] - rmax[0]) / (thdeg[-1] - thdeg[0])
print "Max visible laser beam radius: %f mm + %f mm / deg" %(rmax[0],grad)
pl.plot(thdeg, rmax,label='exact')
pl.plot((thdeg[0],thdeg[-1]),(rmax[0],rmax[0] + grad*(thdeg[-1]-thdeg[0])),'--',label='linear approximation')
pl.legend(loc="lower right")
pl.xlabel('view angle (deg)')
pl.ylabel('max visible laser beam radius')
pl.show()

