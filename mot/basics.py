import scipy as sp


hbar = 1.05457148e-34
kB = 1.3806503e-23
amu = 1.66053886e-27
M = 133 * amu
c = 299792458
l0 = 852e-9
G = 2 * sp.pi * 5.22e6
kappa0 = 3.3e-19

#def wave(T,M):
#    ''' Thermal de Broglie wavelength '''
#    return hbar * sp.sqrt(2 * sp.pi) / sp.sqrt(M * kB * T)
#
#T = 50e-6
#L = wave(50e-6, M)
#n = 2.612 / L**3 / 100**3
#print "Bose-Einstein condensation density for T = %.1f uK : %.2e atoms/cm^3" %(T*1e6,n)

###
#T = 100e-6
#v = sp.sqrt(3 * kB * T / M)
#
#print "Velocity at %f K : %f m/s" %(T,v)
#df = v / l0 / 1e3
#print "Doppler shift equivalent to velocity: %.2f kHz" %(df)

T = 10e-6
d = G * 3
b = 5
kappa = kappa0 * G / d * b
print kappa
r = sp.sqrt( kB * T / kappa)
print r/1e-6