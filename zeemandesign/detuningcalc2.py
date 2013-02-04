"""
Calculations about atom-photon interaction during Zeeman slowing,
focusing on the behaviour of the detuning as a function of position
"""
import numpy as np
import pylab as pl

# Own modules
import layeroptimize as lo
import zeemanslower as zs

# Solve cubic equation
## http://adorio-research.org/wordpress/?p=4493
## Check whether it is actually working well
from math import *
def polyCubicRoots(a,b, c):
    # print "input=", a,b,c
    aby3 = a / 3.0
    p = b - a*aby3
    q = (2*aby3**2- b)*(aby3) + c
    X = (p/3.0)**3
    Y = (q/2.0)**2
    Q = X + Y
    D = 18*a*b*c-4*a**3*c+a**2*b**2-4*b**3-27*c**2
    # print "Q=", Q, D
    if Q >= 0:
       sqQ = sqrt(Q)
       A0 = -q/2.0 + sqQ
       A = A0**(1/3.0) if A0 > 0 else -(-A0)**(1/3.0)
       B0 = -q/2.0 - sqQ
       B = B0**(1/3.0) if B0 > 0 else -(-B0)**(1/3.0)
       r1 = A + B - aby3
       re = -(A+B)/2.0-aby3
       im = sqrt(3.0)/2.0*(A-B)
       r2 = re+im*1j
       r3 = re-im*1j
    else:
       # This part has been tested.
       p3by27= sqrt(-p**3/27.0)
       costheta = -q/2.0/ p3by27
       # print "@@@ costheta=", costheta
       alpha = acos(costheta)
       mag = 2 * sqrt(-p/3.0)
       alphaby3 = alpha/3.0
       r1 = mag  * cos(alphaby3) - aby3
       r2 = -mag * cos(alphaby3+ pi/3)-aby3
       r3 = -mag * cos(alphaby3- pi/3) -aby3
    return r1, r2, r3


def getzeros(params):
    atom, detu, s, B, dBdz = params
    A = zs.uprimehbar*B
    X = A - detu
    Y = 1 + s
    Z = 4 / atom.G**2
    W = -zs.hbar*atom.k**3*atom.G*s/(2*atom.m)
    k = atom.k

    a = X
    b = Y/Z
    c = X*Y/Z+W/(A*Z)
    delta = np.linspace(-5, 3.5, 100)*atom.G

    def dddz(d):
        # d *= atom.G
        # v = (d - detu + zs.uprimehbar * B)/k
        # a = -zs.hbar*k*atom.G/(2*atom.m)*(s/(1+s+4*d**2/atom.G**2))
        # out = k*a/v-zs.uprimehbar*dBdz
        # out = -zs.hbar*k**3*atom.G*s/(2*atom.m)/((1+s+4*d**2/atom.G**2)*(d - detu + zs.uprimehbar*B)) - zs.uprimehbar*dBdz

        W = -zs.hbar*k**3*atom.G*s/(2*atom.m)
        A = zs.uprimehbar*dBdz
        X = -detu + zs.uprimehbar*B
        Y = 1 + s
        Z = 4/atom.G**2
        # out = W/((d + X)*(Y+Z*d**2)) - A
        a, b, c = X, Y/Z, X*Y/Z-W/A/Z
        p = [1, X, Y/Z, X*Y/Z-W/A/Z]
        # out = np.polyval(p, d)
        res = polyCubicRoots(a,b,c)
        out = np.array([x/atom.G for x in res])
        return out

    roots = dddz(delta)
    return roots

if __name__ == "__main__":
    atom = zs.Rb87()
    
    v0 = 254
    vf = 30
    eta = 0.5
    # z0 = 0
    # vmin, vmax = 230, 270

    z0 = zs.slowerlength(atom.aslow, eta, v0, vf)
    z0 = 0.5683
    vmin, vmax = 3, 130


    detu = -175*1e6*2*np.pi
    B0 = zs.bideal(atom, [z0], eta, v0, vf, detu/(-2e6*np.pi))
    vrange = np.linspace(vmin, vmax, 50000)
    delta = atom.k*vrange + detu - zs.uprimehbar * B0
    s = eta/(1-eta)
    s = 2
    print s
    a = -zs.hbar*atom.k*atom.G/(2*atom.m) * s/(1 + s + 4*delta**2 / atom.G**2)
    dz = 0.001
    B = zs.bideal(atom, [z0], eta, v0, vf, detu/(-2e6*np.pi))
    dBdz = np.diff(zs.bideal(atom, [z0, z0+dz], eta, v0, vf, detu/(-2e6*np.pi)))/dz
    # dBdz = 0
    # dBdz = zs.hbar*atom.k*atom.G/(2*atom.m) * s/(1 + s)
    print B, dBdz
    ddeltaz = atom.k*a/vrange - zs.uprimehbar*dBdz
    # pl.plot(vrange, delta)
    i = np.argmin(ddeltaz)
    xx = delta[i]/atom.G
    xxval = ddeltaz[i]/atom.G
    print("Min value at %g" %(xx))
    fig = pl.figure(figsize=(11.69, 8.27))
    ax1 = fig.add_subplot(111)
    ax1.plot(delta/atom.G, ddeltaz/atom.G, 'k-', linewidth=5)
    ax1.plot([delta[0]/atom.G, delta[-1]/atom.G], [0, 0], 'k--', linewidth=2)
    # ax1.set_xlim([delta[0]/atom.G, delta[-1]/atom.G])
    ax1.set_xlabel("Detuning $\delta$ (units of $\Gamma$)", fontsize=16)
    ax1.set_ylabel("Gradient $\partial \delta / \partial z$ (units of $\Gamma$/m)", fontsize=16)
    ax1.tick_params(labelsize='large')
    # ax1.set_ylim([-50, 50])

    params = atom, detu, s, B, dBdz
    roots = getzeros(params)
    # for x in roots:
    #     if x > 0:
    #         dx = 50
    #     else:
    #         dx = -50
    #     pl.annotate("%.3f" %(x), (x, -3),  xycoords='data', xytext=(dx, -20), textcoords='offset points', arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=90,rad=10"), ha='center', fontsize=14)
    #     # pl.plot(x, 0, 'bs', markersize=10)


    # pl.arrow(-1.4, 45, 0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # pl.arrow(0.4, 45, -0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # pl.arrow(0.6, 45, 0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # x1 = roots[0]
    # x2 = roots[1]
    # pl.plot([x1, x1], [-50, 50], 'k:')
    # pl.plot([x2, x2], [-50, 50], 'k:')



    ax2 = ax1.twiny()
    ax2.set_xlim([vrange[0], vrange[-1]])
    ax2.set_xlabel("velocity ($m/s$)", fontsize=16)
    ax2.tick_params(labelsize='large')

    # pl.savefig("detu3.pdf")

    pl.show()
    
