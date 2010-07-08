from __future__ import division
from numpy import *
from pylab import *

gn = 1800
d = 1e-3/gn
thi = linspace(0, pi/2, 200)
m = 1
l = 840e-9
dosave = False

def grate(thi, d, m, l):
    stho = m*l/d + sin(thi)
    thi2 = thi[abs(stho) < 1]
    tho2 = arcsin(stho[abs(stho) < 1])
    return (thi2, tho2)

thi0, tho0 = grate(thi, d, m, l)
thi1, tho1 = grate(thi, d, m, l-15e-9)
thi2, tho2 = grate(thi, d, m, l+15e-9)

figure(1)
plot(thi1/pi*180, tho1/pi*180, ':', label="short")
plot(thi0/pi*180, tho0/pi*180, '-', label="centre")
plot(thi2/pi*180, tho2/pi*180, '--', label="long")
xlabel("Incidence angle (deg)")
ylabel("Diffraction angle (deg)")
title('Groove number: %d' %gn)
legend(loc="best")
if dosave:
    savefig("diff1_%d.eps" %(gn))

figure(2)
plot(thi1/pi*180, (thi1+tho1)/pi*180, ':', label="short")
plot(thi0/pi*180, (thi0+tho0)/pi*180, '-', label="centre")
plot(thi2/pi*180, (thi2+tho2)/pi*180, '--', label="long")
xlabel("Incidence angle (deg)")
ylabel("Angle between incidence and diffraction(deg)")
title('Groove number: %d' %gn)
legend(loc="best")
if dosave:
    savefig("diff2_%d.eps" %(gn))


# figure(2)
# plot(thi0/pi*180, (tho0 - tho2)/pi*180)
# plot(thi0/pi*180, -1*(tho0 - tho1)/pi*180)

figure(3)
ll = linspace(l-15e-9, l+15e-9)
thi = 45.0/180*pi
f = 0.1
x = tan(arcsin(m*ll/d - sin(thi))-arcsin(m*l/d - sin(thi)))*f
plot(ll/1e-9, x/1e-3)
xlabel('Wavelength component (nm)')
ylabel('Distance from centre frequency position (mm)')
title('Dispersion at 10cm from grating, groove number: %d' %gn)
if dosave:
    savefig("diff3_%d.eps" %(gn))

show()
