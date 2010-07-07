from __future__ import division
from numpy import *

thi = 45/180*pi
thd = 45/180*pi
l = 840e-9
f = 0.10
wi = 0.5e-3
c = 3e8
d = 1e-3/1200
w = 2*pi*c / l

w0 = cos(thi)/cos(thd) * (f * l / (pi*wi))
print "Focused beam size  ~ %.0f um " %(w0*1e6)

alpha = l*l*f/(2*pi*c*d*cos(thd))
print "Spatial dispersion ~ %e mm (rad / s)^-1" %(alpha)

# x = (w0/(sqrt(8)*alpha))
# t = linspace(-50,50,101)*1e-15
# g = exp(-(x*t)**2)
# plot(t, g-1)
# show()

T = 4*alpha*sqrt(log(2))/w0
print "Temporal window ~ %e s" %T

# df = 0.44/T
# print "Smalles spectral feature ~ %e" %(df)

w2 = 2*pi*c / (l+15e-9)
w1 = 2*pi*c / (l-15e-9)

df = (w1 - w2)/128
print "Frequency change per pixel: %e" %(df/1e12)
