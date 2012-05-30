#!/usr/bin/env python
"""
Simulation and data fitting to design a Zeeman slower.
Follows in the footsteps of Bell et.al., Rev. Sci. Inst 81 (2010) 013105
"""
from __future__ import division
import numpy as np
import pylab as pl
import scipy.integrate as integ
from mpl_toolkits.mplot3d import Axes3D
import scipy.odr as odr
from time import time

import zeemanslower as zs

## Constants:
mu = 1.26e-6*1.5e-3
hbar = 1.05457148e-34

fw, fh = 11.69, 8.27
savefigs = True
figname = "vpitch01"

def copyprint(thisarray):
    print "[%s]" %(", ".join(["%g"%val for val in thisarray]))

def pos(p, params):
    """ parameteric coil """
    cord = params['cord']
    cp = params['cp']
    cn = params['cn']
    R = params['R']
    thetap = np.polyval(cp[0:(cord+1)], p)
    xp = R * np.cos(thetap)
    yp = R * np.sin(thetap)
    zp = np.polyval(cp[(cord+1):(cord+3)], p)
    thetan = np.polyval(cn[0:(cord+1)], p)
    xn = R * np.cos(thetan)
    yn = R * np.sin(thetan)
    zn = np.polyval(cn[(cord+1):(cord+3)], p)
    return (xp, yp, zp, xn, yn, zn)

def deriv(p, params):
    """ dr/dp """
    c = params['c']
    # The coefficients of the polynomial after derivation
    cnew = [c[i]*(6-i) for i in range(0,6)]
    R = params['R']
    theta = np.polyval(c[0:7], p)
    cpderiv = np.polyval(cnew, p)
    dx = -R * np.sin(theta) * cpderiv
    dy = R * np.cos(theta) * cpderiv
    dz = c[8] + p*0
    return zip(dx, dy, dz)

def bzfield(p, params, z):
    cord = params['cord']
    cp = params['cp']
    cn = params['cn']
    cpnew = [cp[i]*(cord-i) for i in range(0,cord)]
    cnnew = [cn[i]*(cord-i) for i in range(0,cord)]
    R = params['R']
    I = params['I']
    upperp = np.polyval(cpnew, p)
    lowerp = (R**2 + (np.polyval(cp[cord+1:], p)-z)**2)**(3/2)
    uppern = np.polyval(cnnew, p)
    lowern = (R**2 + (np.polyval(cn[cord+1:], p)-z)**2)**(3/2)
    return mu*I/(4*np.pi)*(upperp/lowerp + uppern/lowern)

def plotcoils(params):
    prange = np.linspace(0, 2*np.pi, 40001)
    xp, yp, zp, xn, yn, zn = pos(prange, params)
    pl.figure(figsize=(fw, fh))
    pl.plot(zp, yp, label="Positive coil")
    pl.plot(zn, yn, label="Negative coil")
    pl.xlabel('z position (m)', fontsize=16)
    pl.ylabel('y position (m)', fontsize=16)
    pl.text(0, 0.03, 
            "Positive coil", 
            fontsize=17, 
            horizontalalignment='center',
            verticalalignment='center',
            bbox=dict(facecolor='white', alpha=0.8)
            )
    pl.text(0.55, 0.03, 
            "Negative coil", 
            fontsize=17,
            horizontalalignment='center',
            verticalalignment='center',
            bbox=dict(facecolor='white', alpha=0.8)
            )
    pl.title('Coils - side view')
    pl.xlim([zp[-1]-0.01, zn[-1]+0.01])
    # pl.legend(loc="center left")

def fitfunction(c, z, params, limit=True):
    if type(z) == type(1):
        z = np.array([z])
    ## Coil parameters
    cp = c[0:9]
    cn = c[9:18]
    I = params['I']
    R = params['R']
    cord = params['cord']
    ## Input parameters
    paramsin = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
    prange = np.linspace(0, 2*np.pi, 4001)
    res = np.array([integ.trapz(bzfield(prange, paramsin, zz), prange) for zz in z])
    pitch = getpitch(paramsin)
    plim = 4e-3
    if (pitch[0] < plim) or (pitch[1] < plim) :
        print ("Too small pitch")
        if limit:
            res *= 0
    return res

def dofit(params, z, bfield):
    data = odr.Data(z, bfield)
    fitfunc = lambda x, z: fitfunction(x, z, params, True)
    model = odr.Model(fitfunc)
    c0 = np.append(params['cp'], params['cn'])
    fit = odr.ODR(data, model, c0)
    fit.set_job(fit_type=2)
    out = fit.run()
    return out

def getpitch(params):
    cp = params['cp']
    cn = params['cn']
    cord = params['cord']
    prange = np.array([2*np.pi, 2*np.pi-1e-3])
    prange = np.array([2*np.pi, 2*np.pi-1e-3])
    # print prange
    # zp = np.polyval(cp[(cord+1):(cord+3)], prange)
    thetap = np.polyval(cp[0:(cord+1)], prange)
    zp = np.polyval(cp[(cord+1):(cord+3)], prange)
    pitchp = -np.diff(zp) / (np.diff(thetap) / (2*np.pi))

    thetan = np.polyval(cn[0:(cord+1)], prange)
    zn = np.polyval(cn[(cord+1):(cord+3)], prange)
    pitchn = -np.diff(zn) / (np.diff(thetan) / (2*np.pi))
    # print (zp, prange)
    return (pitchp[0], pitchn[0])

def coillength(params):
    cp = params['cp']
    cn = params['cn']
    coord = params['cord']
    R = params['R']
    prange = np.linspace(0, 2*np.pi, 40001)
    thetap = np.polyval(cp[0:(cord+1)], prange)
    xp = R * np.cos(thetap)
    yp = R * np.sin(thetap)
    zp = np.polyval(cp[(cord+1):(cord+3)], prange)

    dx, dy, dz = np.diff(xp), np.diff(yp), np.diff(zp)
    lenp = np.sum(np.sqrt(dx**2 + dy**2 + dz**2))

    thetan = np.polyval(cn[0:(cord+1)], prange)
    zn = np.polyval(cn[(cord+1):(cord+3)], prange)
    xn = R * np.cos(thetan)
    yn = R * np.sin(thetan)

    dx, dy, dz = np.diff(xn), np.diff(yn), np.diff(zn)
    lenn = np.sum(np.sqrt(dx**2 + dy**2 + dz**2))
    
    return (lenp, lenn)

if __name__ == "__main__":
    cord = 6
    R = 0.0383

    # ### Fixing up the numbers for Bell2010
    # cp = [0, 1.58e-3, 0, -0.30231, 8.69812, 2.99942, 0, -0.10626, 0.58149]
    # cn = [-0.0012, 0.00347, 0.02038, -0.23045, -1.25728, -1.9995, np.pi, 0.04308, 0.61819]

    # Make the pitch >4mm by adjusting cn[7]
    # cp = [0, 1.58e-3, 0, -0.30231, 8.69812, 2.99942, 0, -0.10626, 0.58149]
    # cn = [-0.0012, 0.00347, 0.02038, -0.23045, -1.25728, -1.9995, np.pi, 0.04508, 0.61819]

    # I, eta, v0, vf, detu = 110, 0.7, 360, 20, 252
    # params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}



    I, eta, v0, vf, detu = 150, 0.5, 257, 30, 160
    #### From program
    c = [0.00016638, -0.00432291, 0.0540347, -0.522765, 8.47316, 1.22958, -0.00141238, -0.101068, 0.568002, -0.00398553, 0.0506522, -0.257939, 0.540945, -2.33171, -0.575496, -0.0234875, 0.0440871, 0.620623]
    #### From fit:
    c = [0.00016641, -0.0043198, 0.0541086, -0.52173, 8.51442, 1.22973, -0.00141238, -0.100587, 0.570773, -0.0040054, 0.0501239, -0.260157, 0.539409, -2.33601, -0.575532, -0.0234875, 0.0438557, 0.613775]
    ####

    I, eta, v0, vf, detu = 110, 0.5, 257, 30, 160
    ### From program
    c = [8.67375e-05, -0.0020817, 0.0235791, -0.199735, 2.71148, 0.377943, -0.000588783, -0.0442239, 0.234247, -0.00146283, 0.0151475, -0.0644685, 0.104234, -0.587468, -0.186212, -0.00303369, 0.0180984, 0.248207]
    ### From fit
    c = [8.68791e-05, -0.00206882, 0.0238349, -0.196968, 2.78631, 0.378165, -0.000588783, -0.0434768, 0.237525, -0.00146947, 0.0150472, -0.064683, 0.104205, -0.587, -0.186184, -0.00303369, 0.0181244, 0.251633]
    ###
    c = [0.000150134, -0.00374181, 0.0436853, -0.376179, 5.09255, 1.59556, -0.00141396, -0.0749602, 0.413539, -0.000867299, 0.00955619, -0.0449806, 0.0697088, -0.758922, -0.0343091, -0.00292828, 0.026442, 0.461673]
    c = [0.000186918, -0.00232913, 0.0422248, -0.467452, 5.95027, 1.54941, -0.00141396, -0.0811131, 0.433969, -0.00140557, 0.0083009, -0.0207248, 0.0802116, -0.990341, -0.0346022, -0.00292828, 0.0263718, 0.451932]
    ###
    cp = c[0:9]
    cn = c[9:18]
    params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}

    # #### Fit:
    # # c = [  1.79648369e-02,   4.97082071e-01,   4.54901216e-04,   4.58564916e+02,
    # #        6.03699716e+04,   1.14109895e+03,   0.00000000e+00,   8.83821139e+01,
    # #        7.74853503e+02,   1.71772275e+00,   2.29233050e+00,   1.25960930e+01,
    # #        2.55986134e+02,   1.21485029e+03,   4.89355700e+02,   3.14159265e+00,
    # #        -2.99694277e-02,   8.65410378e-04]
    # # c = [  1.73266815e-03,  -1.41892193e-02,   1.97767779e-03,   9.19562958e-02,
    # #        5.61661434e+00,   5.07707741e+00,   0.00000000e+00,  -9.47365984e-02,
    # #        5.28805144e-01,  -1.04767200e-02,   9.19535619e-02,  -5.47788716e-02,
    # #        -1.55237738e+00,   3.09374681e+00,  -6.91959775e+00,   3.14159265e+00,
    # #        5.59338290e-02,   4.93951310e-01]

    # # c = [  1.72623824e-03,  -1.39936875e-02,  -5.26044608e-05,   1.01619072e-01,
    # #        5.59016985e+00,   5.11963808e+00,   0.00000000e+00,  -9.47001492e-02,
    # #        5.28614669e-01,  -8.99131445e-03,   6.43459938e-02,   1.38871640e-01,
    # #        -2.18940123e+00,   4.07914091e+00,  -7.57512654e+00,   3.14159265e+00,
    # #        4.18215913e-02,   5.94986034e-01]

    # # Fit with detuning I, v0, vf, detu = 110, 360, 20, 260
    # I, eta, v0, vf, detu = 110, 0.7, 360, 20, 260

    # # c = [ -2.64813662e-08,   1.57917280e-03,  -7.45867070e-10,  -3.03593810e-01,
    # #        8.42031791e+00,   2.98830152e+00,   0.00000000e+00,  -1.07145396e-01,
    # #        5.41570863e-01,  -1.21911992e-03,   3.44918522e-03,   2.02963812e-02,
    # #        -2.31548043e-01,  -1.26169437e+00,  -2.00422354e+00,   3.14159265e+00,
    # #        4.44841435e-02,   5.82735817e-01]

    # # Startting parameter from me
    # c = [8.43598e-05, -0.00240179, 0.0332447, -0.362301, 6.74412, 0.492459, -0.00134721, -0.091169, 0.496232, -0.00122632, 0.0186904, -0.11525, 0.288536, -1.68802, -0.5132, -0.026239, 0.0387468, 0.663752]

    # # # Fit with detuning I, v0, vf, detu = , 360, 20, 100
    # I, eta, v0, vf, detu = 80, 0.7, 230, 20, 150
    # c = [7.3407e-05, -0.00200823, 0.0259163, -0.249601, 3.83613, 0.389881, -0.00137507, -0.0387378, 0.166797, -0.000276305, 0.00407555, -0.0249578, 0.056904, -0.513007, -0.21706, -0.00583082, 0.0126123, 0.325657]
    # c = [0.000146814, -0.00401645, 0.0518326, -0.499202, 7.67226, 0.779761, -0.00275013, -0.0387378, 0.166797, -0.000552609, 0.00815109, -0.0499155, 0.113808, -1.02601, -0.43412, -0.0116616, 0.0126123, 0.325657]
    # c = [7.3407e-05, -0.00200823, 0.0259163, -0.249601, 3.83613, 0.389881, -0.00137507, -0.0387378, 0.185947, -0.000276305, 0.00407555, -0.0249578, 0.056904, -0.513007, -0.21706, -0.00583082, 0.0126123, 0.287357]

    # c = [7.27145e-05, -0.00199057, 0.0257103, -0.248054, 3.82638, 0.486438, -0.00123624, -0.0387378, 0.185947, -0.00057977, 0.00858154, -0.0524891, 0.120679, -1.03494, -0.452585, -0.0114505, 0.0126123, 0.267357]


    # # c = [ -2.64813662e-08,   0.77917280e-03,  -7.45867070e-10,  -6.03593810e-01,
    # #        8.42031791e+00,   2.98830152e+00,   0.00000000e+00,  -0.47145396e-01,
    # #        2.01570863e-01,  -1.61911992e-03,   6.44918522e-03,   7.02963812e-02,
    # #        -2.31548043e-01,  -2.66169437e+00,  -3.00422354e+00,   3.14159265e+00,
    # #        2.74841435e-02,   2.42735817e-01]

    # # c = [ -2.64813649e-08,   7.79280559e-04,  -7.45867070e-10,  -6.06985286e-01,
    # #        8.20684733e+00,   2.98924947e+00,   0.00000000e+00,  -4.87250744e-02,
    # #        2.11183230e-01,  -1.63554003e-03,   6.39867407e-03,   6.91223878e-02,
    # #        -2.33956884e-01,  -2.71102473e+00,  -3.00302504e+00,   3.14159265e+00,
    # #        2.70037662e-02,   2.62210776e-01]


    # # c = [-4.02676e-20, 1.19646e-07, -2.84578e-06, 3.1362e-05, -0.000255343, 0.00320356, 0.000283923, -0.0903782, 0.529563, 7.67097e-21, -2.15342e-09, -1.16645e-07, -2.91219e-06, -5.2993e-05, -0.00145625, -1.87751e-05, 0.0404679, 0.614639]

    # c = [7.27762e-05, -0.00198243, 0.0259533, -0.243921, 4.01206, 0.487037, -0.00123624, -0.0364344, 0.227223, -0.000582171, 0.00849813, -0.0529762, 0.120291, -1.03882, -0.452648, -0.0114505, 0.0126061, 0.270278]

    # c = [0.000145429, -0.00398115, 0.0514206, -0.496107, 7.65275, 0.972877, -0.00247248, -0.0387378, 0.224247, -0.00057977, 0.00858154, -0.0524891, 0.120679, -1.03494, -0.452585, -0.0114505, 0.0126123, 0.268207]

    # c = [0.000109072, -0.00298586, 0.0385655, -0.37208, 5.73957, 0.729657, -0.00185436, -0.0387378, 0.224247, -0.000434827, 0.00643615, -0.0393668, 0.0905094, -0.776202, -0.339439, -0.00858791, 0.0126123, 0.268207]

    # c = [0.000109072, -0.00298586, 0.0385655, -0.37208, 5.73957, 0.729657, -0.00185436, -0.0387378, 0.224247, -0.000724712, 0.0107269, -0.0656114, 0.150849, -1.29367, -0.565732, -0.0143132, 0.0126123, 0.268207]

    # c = [8.18038e-05, -0.00223939, 0.0289241, -0.27906, 4.30467, 0.547243, -0.00139077, -0.0387378, 0.224247, -0.000434828, 0.00643616, -0.0393669, 0.0905095, -0.776202, -0.339439, -0.00858792, 0.0126123, 0.268207]

    # # I, eta, v0, vf, detu = 60, 0.7, 170, 20, 150

    # # c= [1.89457e-06, -6.4773e-05, 0.00123203, -0.0224159, 0.826359, 0.354265, -1.57742e-05, -0.0146715, 0.0730335, -0.000775362, 0.0114645, -0.0701438, 0.161384, -1.38732, 0.0592676, -0.0143273, 0.0132043, 0.114406]
 
    # cp = c[0:9]
    # cn = c[9:18]
    # params = {'R': R, 'cp': cp, 'cn': cn, 'cord': cord, 'I': I}
    # testrun = False
    testrun = True

    ##### Parameters
    nz = 33 # number of points, including the +2 extra points in the ens outside
    # eta, v0, vf, detu = 0.7, 230, 20, 150
    #####

    atom = zs.Rb85()
    sl = zs.slowerlength(atom.aslow, eta, v0, vf)
    print sl
    z = np.append([-5.5*params['R']], np.append(np.linspace(0, sl, nz-2), sl+5.5*params['R']))
    bfield = zs.bideal(atom, z, eta, v0, vf, detu)


    try:
        testrun
    except NameError:
        testrun = False

    if testrun:
        c = np.append(params['cp'], params['cn'])
    else:
        out = dofit(params, z, bfield)
        out.pprint()
        copyprint(out.beta)
        cfit = out.beta
        cp = cfit[0:9]
        cn = cfit[9:18]
        params['cp'] = cp
        params['cn'] = cn
        c = cfit

    plotcoils(params)
    if savefigs:
        pl.savefig("%s_coils.png" %(figname))

    zdetail = np.linspace(z[0], z[-1], 1001)
    ff = fitfunction(c,
                     zdetail,
                     params,
                     False)

    pl.figure(figsize=(fw, fh))
    # pl.subplot(2,1,1)
    ax = pl.plot(zdetail, ff*1e4, 'g-', linewidth=3, label='Generated field')
    pl.plot(z, bfield*1e4, 'ko', markersize=7, label='Target field')
    pl.legend(loc='best')
    pl.xlim([z[0], z[-1]])
    pl.xlabel('z position (m)', fontsize=16)
    pl.ylabel('Magnetic field (G)', fontsize=16)
    pl.title("R=%.3fm, I=%.1fA, v0=%.1fm/s, vf=%0.1fm/s D=%.1fMHz" %(R, I, v0, vf, detu))
    # pl.subplot(2,1,2)
    # pl.plot(z, ff-bfield)
    if savefigs:
        pl.savefig("%s_field.png" %(figname))

    ## Analyze the coil
    print "Minimum pitch: ", getpitch(params)
    clen = coillength(params)
    print "Coil lengths: ", clen
    resistance = 0.334 / 100
    rp = clen[0] * resistance
    rn = clen[1] * resistance
    print "Resistance: ", (rp, rn)
    # powerdiss = rp*110**2 + rn*90**2
    powerdiss = rp*110**2 + rn*110**2
    voltage = (rp+rn)*I
    print "Power discipation:", powerdiss 
    print "Voltage drop:", voltage


    ### Draw pitch as function of position. Pitch is dz/d(theta/2pi)
    prange = np.linspace(0, 2*np.pi, 1001)
    thetap = np.polyval(cp[0:(cord+1)], prange)
    zp = np.polyval(cp[(cord+1):(cord+3)], prange)
    pitchp = -np.diff(zp) / (np.diff(thetap) / (2*np.pi))

    thetan = np.polyval(cn[0:(cord+1)], prange)
    zn = np.polyval(cn[(cord+1):(cord+3)], prange)
    pitchn = -np.diff(zn) / (np.diff(thetan) / (2*np.pi))

    pl.figure(figsize=(fw, fh))
    pl.semilogy(zp[1:], pitchp, linewidth=3, label='positive coil')
    pl.semilogy(zn[1:], pitchn, linewidth=3, label='negative coil')
    pl.xlim([zp[-1]-0.01, zn[-1]+0.01])
    # pl.ylim([0, 0.01])
    pl.xlabel('z position (m)', fontsize=16)
    pl.ylabel('Pitch (m)', fontsize=16)
    pl.legend(loc='upper left')
    if savefigs:
        pl.savefig("%s_pitch.png" %(figname))
    ####


    # Show everything
    pl.show()
