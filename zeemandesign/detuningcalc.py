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

# Predict dDelta/dZ = 0 points



if __name__ == "__main__":
    atom = zs.Rb87()
    
    v0 = 254
    vf = 30
    eta = 0.5
    z0 = 0
    detu = -175*1e6*2*np.pi
    B0 = zs.bideal(atom, z0, eta, v0, vf, detu/(-2e6*np.pi))
    vmin, vmax = 230, 270
    vrange = np.linspace(vmin, vmax, 50000)
    delta = atom.k*vrange + detu - zs.uprimehbar * B0
    s = eta/(1-eta)
    # s = 2
    print s
    a = -zs.hbar*atom.k*atom.G/(2*atom.m) * s/(1 + s + 4*delta**2 / atom.G**2)
    dz = 0.001
    dBdz = np.diff(zs.bideal(atom, [z0, z0+dz], eta, v0, vf, detu/(-2e6*np.pi)))/dz
    # dBdz = 0
    # dBdz = zs.hbar*atom.k*atom.G/(2*atom.m) * s/(1 + s)
    print dBdz
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
    ax1.set_xlim([delta[0]/atom.G, delta[-1]/atom.G])
    ax1.set_xlabel("Detuning $\delta$ (units of $\Gamma$)", fontsize=16)
    ax1.set_ylabel("Gradient $\partial \delta / \partial z$ (units of $\Gamma$/m)", fontsize=16)
    ax1.tick_params(labelsize='large')
    ax1.set_ylim([-50, 50])


    pl.arrow(-1.4, 45, 0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # pl.arrow(0.4, 45, -0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # pl.arrow(0.6, 45, 0.7, 0, color='k', width=1.0, head_width=3.0, head_length=0.1)
    # x1 = -0.5
    # x2 = 0.5
    # pl.plot([x1, x1], [-50, 50], 'k:')
    # pl.plot([x2, x2], [-50, 50], 'k:')



    ax2 = ax1.twiny()
    ax2.set_xlim([vrange[0], vrange[-1]])
    ax2.set_xlabel("velocity ($m/s$)", fontsize=16)
    ax2.tick_params(labelsize='large')

    # pl.savefig("detu2.pdf")

    pl.show()
    
