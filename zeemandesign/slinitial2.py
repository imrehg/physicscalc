from __future__ import division
import numpy as np
import pylab as pl
import zeemanslower as zs
from scipy.integrate import trapz
import vpitchcoil as vp
from scipy.interpolate import interp1d

# mu = 1.26e-6*1.5e-3
mu = 1.2566370614e-6

# Generate only positive or only negative part of the field

# Calculate z(p)
# Calculate B(z)
# Calculate theta(p) ~ dz/dp*B(z(p)-R)
# Compare with real parameters
# Use it to make initial parameters for dimensions less than 6 (e.g. c0-c4, c0-c3)

def copyprint(thisarray):
    print "[%s]" %(", ".join(["%g"%val for val in thisarray]))


def getfield(prange, zpr, c6, cz, R, I, mu):
    cord = len(c6)-1
    cnew = [c6[i]*(cord-i) for i in range(0,cord)]
    testfield = np.zeros(len(zpr))
    for i, z in enumerate(zpr):
        upperp = np.polyval(cnew, prange)
        lowerp = (R**2 + (np.polyval(cz, prange)-z)**2)**(3/2)
        fields = mu*I/(4*np.pi)*upperp/lowerp
        testfield[i] = trapz(fields, prange)
    return testfield


def getcoeff(z, b, mu=1.2566370614e-6, I=1, R=1, shifted=0):

    pzp = np.polyfit(z, b/mu/I, 10)
    p8 = z[0]
    # ztop = z[-1]+np.sign(z[-1]-z[-2])*2*R  # Make things a little bit longer
    ztop = z[-1]
    p7 = (ztop-p8)/2.0/np.pi
    prange = np.linspace(0, 2*np.pi, 101)
    zpr = np.polyval([p7, p8], prange)
    dtheta = np.sign(p7)*np.diff(zpr)*np.polyval(pzp, zpr)[1:]*2*np.pi
    cumtheta = [sum(dtheta[:(i+1)]) for i in range(len(dtheta))]
    cumtheta = np.append(np.array([0]), cumtheta)

    c6 = np.polyfit(prange, cumtheta, 6)
    cz = [p7, p8]

    extrarange = 1.25
    extrarange = 1
    ztop2 = z[-1]+np.sign(z[-1]-z[-2])*extrarange*R  # Make things a little bit longer
    p72 = (ztop2-p8)/2.0/np.pi
    cz = [p72, p8]
    pmax = (ztop2 - p8) / p7
    print pmax/2/np.pi
    newcumtheta = np.polyval(c6, np.linspace(0, pmax, len(prange)))
    c6 = np.polyfit(prange, newcumtheta, 6)

    zpr = np.polyval(cz, prange)


    cp = np.append(c6, cz)

    mu2 = 1.26e-6*1.5e-3
    # # print 'CP ', cp

    # # print cz
    # cp = [-0.000132412, 0.00181674, -0.0100659, 0.0179691, -0.220672, -0.472009, -0.000323899, 0.0113961, 0.778152]
    # c6 = cp[0:7]
    # cz = cp[7:]
    # thtotal = np.polyval(c6, 2*np.pi) / 2/np.pi
    # print "Total turns:", thtotal

    # testfield = getfield(prange, zpr, c6, cz, R, I, mu)

    # cord = len(c6)-1
    # # params = {'R': R, 'cp': cp, 'cn': [0 for val in cp], 'cord': cord, 'I': I}  # only positive coil
    # params = {'R': R, 'cp': [0 for val in cp], 'cn': cp, 'cord': cord, 'I': I}  # only negative coil
    # otherfield = [trapz(vp.bzfield(prange, params, zv), prange) for zv in zpr]
    # # s1 = sum(vp.bzfield(prange, params, 0))
    # # s2 = sum(getfield(prange, zpr, c6, cz, R, I, mu))

    # for i in xrange(1):
    #     c6 = np.polyfit(prange, cumtheta, 6)

    #     cz = [p7, p8]
    #     cp = np.append(c6, cz)

    #     mu2 = 1.26e-6*1.5e-3

    #     testfield = getfield(prange, zpr, c6, cz, R, I, mu2)
    #     target = np.polyval(pzp, zpr)*mu*I
    #     cumtheta = cumtheta * target/testfield 

    # pl.figure()

    # pl.plot(zpr, cumtheta)

    # pl.figure()

    # c6 = np.polyfit(prange, cumtheta, 6)
    testfield = getfield(prange, zpr, c6, cz, R, I, mu2)
    pl.plot(zpr, testfield, 'r-')
    pl.plot(z, b, '--')
    # pl.plot(zpr, np.polyval(pzp, zpr)*mu*I, 'kx-')

    # pl.plot(prange, cumtheta)
    # cord = len(c6)-1
    # # params = {'R': R, 'cp': cp, 'cn': [0 for val in cp], 'cord': cord, 'I': I}  # only positive coil
    # params = {'R': R, 'cp': [0 for val in cp], 'cn': cp, 'cord': cord, 'I': I}  # only negative coil
    # otherfield = [trapz(vp.bzfield(prange, params, zv), prange) for zv in zpr]
    # s1 = sum(vp.bzfield(prange, params, 0))
    # s2 = sum(getfield(prange, zpr, c6, cz, R, I, mu))

    # print s1
    # print s2
    # print testfield
    # print testfield

    # pl.figure()
    # # pl.plot(prange, cumtheta)
    # # pl.plot(prange, testfield)
    # pl.plot(zpr, testfield, 'o')
    # pl.plot(z, b, '--')
    # pl.plot(zpr, otherfield, 'rx-')
    # pl.plot(zpr, cumtheta/3000, 'g:')

    # pl.figure()
    # # pl.plot(z, b)
    # # pl.plot(z, b/mu/I)
    # pl.plot(zpr[1:], np.diff(zpr))
    # pl.plot(zpr)
    # pl.title('dtheta')
    # # pl.plot(z, b/mu/I)
    # # pl.plot(zpr, np.polyval(pzp, zpr))
    # pl.plot(prange, cumtheta)    
    # pl.plot
    # pl.title('theta vs p')

    # print np.polyval(c6, 2*np.pi)
    cz[-1] += shifted
    cout = np.append(c6, cz)
    return cout




#### Fitting the positive part of the field

R = 0.0383
I, eta, v0, vf, detu = 110, 0.7, 365, 20, 252

I, eta, v0, vf, detu = 110, 0.7, 230, 20, 150

I, eta, v0, vf, detu = 110, 0.5, 257, 30, 160

params = {'R' : R}
nz = 61
atom = zs.Rb85()
sl = zs.slowerlength(atom.aslow, eta, v0, vf)
z = np.append([-5.5*params['R']], np.append(np.linspace(0, sl, nz-2), sl+5.5*params['R']))
bfield = zs.bideal(atom, z, eta, v0, vf, detu)

zz = z[bfield > 0]
p8 = zz[-1]
# zz = zz[-1]-zz
bz = bfield[bfield > 0]

# ############ Testing
# cp = [0, 1.58e-3, 0, -0.30231, 8.69812, 2.99942, 0, -0.10626, 0.58149]
# c6 = cp[0:7]
# cz = cp[7:]
# prange = np.linspace(0, 2*np.pi, 101)
# zpr = np.polyval(cz, prange)
# theta = np.polyval(c6, prange)
# # cumtheta = [sum(dtheta[:(i+1)]) for i in range(len(dtheta))]
# # cumtheta = np.append(np.array([0]), cumtheta)
# dtheta = np.diff(theta)[::-1]
# zdth = zpr[1:][::-1]
# dth = interp1d(zdth, dtheta)
# pl.plot(zdth, dtheta / dth(0) )
# pl.plot(zz, bz/bz[0])
# pl.plot([0, 0], [0, 1], 'k-')
# pl.plot([-R, -R], [0, 1], 'r-')
# pl.plot([-2*R, -2*R], [0, 1], 'r-')

# # cord = len(c6)-1
# # params = {'R': R, 'cp': cp, 'cn': [0 for val in cp], 'cord': cord, 'I': I}  # only positive coil
# # otherfield = [trapz(vp.bzfield(prange, params, zv), prange) for zv in zpr]
# ############


pzp = np.polyfit(zz, bz/mu/I, 5)

# pl.plot(z, bfield, 'o')

p7 = -p8/2.0/np.pi
prange = np.linspace(0, 2*np.pi, 1001)
zpr = np.polyval([p7, p8], prange)
# pl.plot(prange, np.polyval(pzp, zpr))
theta = trapz(np.polyval(pzp, zpr), zpr)
print theta
dtheta = -np.diff(zpr)*np.polyval(pzp, zpr)[1:]
# print dtheta
print sum(dtheta)
# dt2 = np.append(np.array([0]), dtheta)
# print dt2
cumtheta = [sum(dtheta[:(i+1)]) for i in range(len(dtheta))]
cumtheta = np.append(np.array([0]), cumtheta)

# # pl.figure()
# # pl.plot(prange, cumtheta)
# # pl.title('theta - p')
c6 = np.polyfit(cumtheta, prange, 6)
# print c6
cnew = getcoeff(zz[::-1], bz[::-1], I=I, R=R, shifted=-0.5*R)
print "Func:", cnew

######## Turn this on when done
zn = z[bfield < 0]
bn = bfield[bfield < 0]
cnew2 = getcoeff(zn, bn, I=I, R=R, shifted=0.5*R)

print "Totalcoil"
print "[%s]" %(", ".join(["%g"%val for val in np.append(cnew, cnew2)]))
########
####################
 

# pl.plot(z, bfield)


# p5 = np.polyfit(zz, bz, 5)
# zd = np.linspace(zz[0], zz[-1], 301)
# print p5

# pl.plot(zz, bz, 'o')
# # pl.plot(zd, np.polyval(p5, zd))

# p7 = -p8/2.0/np.pi
# print p7, p8, p7*2*np.pi+p8-R
# prange = np.linspace(0, 2*np.pi)
# zpr = np.polyval([p7, p8], prange)
# pl.plot(zpr, np.polyval(p5, zpr))
# cnew = np.polyfit(prange, 1/np.polyval(p5, zpr), 6)
# cnew = np.append(cnew, [p7, p8-R])
# print "Positive coil"
# print cnew


# pl.figure()
# zn = z[bfield < 0]
# p8 = zn[0]
# # zn = zn-zn[0]
# bn = bfield[bfield < 0]
# p5n = np.polyfit(zn, bn, 6)
# zdn = np.linspace(zn[0], zn[-1], 301)
# p7 = (zn[-1]-p8)/2.0/np.pi
# # pl.plot(zn, bn, 'o')
# # pl.plot(zdn, np.polyval(p5n, zdn))

# prange = np.linspace(0, 2*np.pi)
# zpr = np.polyval([p7, p8], prange)

# pl.plot(zpr, np.polyval(p5n, zpr))
# pl.plot(zn, bn, '.')
# cnew2 = np.polyfit(prange, 1/np.polyval(p5, zpr), 6)
# print "Negative coil"
# print cnew2
# cnew2 =  np.append(cnew2, [p7, p8+R])
# print cnew2

# print "Totalcoil"
# print "[%s]" %(", ".join(["%g"%val for val in np.append(cnew, cnew2)]))





pl.show()
