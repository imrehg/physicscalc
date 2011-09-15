import numpy as np
import pylab as pl

import pinholes


# pl.figure(1, figsize=(8.27/1.5, 11.69/1.5))
fy = 5
fx = fy * 1.618
savefig = True

pl.figure(1, figsize=(fx, fy))
styles = ['k-', 'r--', 'b:']
for i, r1 in enumerate([0.5, 1, 2]):
    r2 = 1
    l = 10
    maxth = 0.34
    th = np.linspace(0, maxth, 1001)
    area = []
    for thval in th:
        area += [ pinholes.angleArea(thval, r1, r2, l) ]
    area = np.array(area) / (r1*r1*np.pi) * np.cos(th)
    pl.plot(th, area, styles[i], linewidth=3, label=r"$r_1$ =  %g" %(r1))
pl.xlim(0, maxth)
pl.ylim(0, 1.05)
pl.legend(loc='best')
pl.title(r"Effect of entry pinhole size, $r_2$ = %g, $L$ = %g" %(r2, l))
pl.xlabel(r'Escape angle $\theta$', fontsize=14)
pl.ylabel('Proportion of atoms escaping', fontsize=14)
if savefig:
    pl.savefig('pinhole_entrysize.pdf')

# pl.figure(2, figsize=(11.69/1.5, 8.27/1.5))
pl.figure(2, figsize=(fx, fy))
styles = ['k-', 'r--', 'b:']
for i, l in enumerate([5, 10, 20]):
    r1, r2 = 1, 1
    maxth = 0.34
    th = np.linspace(0, maxth, 1001)
    area = []
    for thval in th:
        area += [ pinholes.angleArea(thval, r1, r2, l) ]
    area = np.array(area) / (r1*r1*np.pi) * np.cos(th)
    pl.plot(th, area, styles[i], linewidth=3, label=r"$L$ =  %g" %(l))
pl.xlim(0, maxth)
pl.ylim(0, 1.05)
pl.legend(loc='best')
pl.title(r"Effect of aspect ratio, $r_1 = r_2 =$ %g" %(r1))
pl.xlabel(r'Escape angle $\theta$', fontsize=14)
pl.ylabel('Proportion of atoms escaping', fontsize=14)
if savefig:
    pl.savefig('pinhole_aspect.pdf')


pl.show()
