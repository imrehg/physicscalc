#!/usr/bin/python2
import pylab as pl
import numpy as np
import ourgui

def runanalysis():
    filename = ourgui.openFile(type='npz')

    data = np.load(filename, mmap_mode='r')
    X = data['X']
    Y = data['Y']
    BX = data['BX']
    BY = data['BY']
    BZ = data['BZ']
    pitch = data['pitch']
    windnum = data['windnum']

    phi = np.linspace(0, np.pi*2, 100)

    pl.figure(figsize=(11.69, 8.27), dpi=100)

    BXBZ = np.zeros(BZ.shape)
    BYBZ = np.zeros(BZ.shape)
    for i in xrange(BZ.shape[0]):
        for j in xrange(BZ.shape[1]):
            if abs(BZ[i, j]) > 1e-14:
                BXBZ[i, j] = BX[i, j] / BZ[i, j]
                BYBZ[i, j] = BY[i, j] / BZ[i, j]

    
    pl.subplot(2, 2, 1, aspect='equal')
    pl.title("Pitch = %g, wind number = %d" %(pitch, windnum))
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BXBZ*100, colors='k')
    pl.clabel(CS, inline=1, fontsize=10, fmt='%.1f%%')
    pl.xlabel("Bx / Bz")

    pl.subplot(2, 2, 2, aspect='equal')
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BYBZ*100, colors='k')
    pl.clabel(CS, inline=1, fontsize=10, fmt='%.1f%%')
    pl.xlabel("By / Bz")

    pl.subplot(2, 2, 3, aspect='equal')
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BZ / BZ.max()*100, 10, colors='k')
    pl.clabel(CS, inline=1, fontsize=10, fmt='%.1f%%')
    pl.xlabel("Bz (normalized)")

    pl.subplot(2, 2, 4, aspect='equal')
    FIELD = np.sqrt(BX**2 + BY**2 + BZ**2)
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, FIELD, colors='k')
    pl.clabel(CS, inline=1, fontsize=10)
    pl.xlabel("Total field strength")

    # import matplotlib.cm as cm
    # im = pl.imshow(FIELD, interpolation='bilinear', origin='lower',
    #                 cmap=cm.gray, extent=(-1,1,-1,1))
    # pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    # CS = pl.contour(X, Y, np.sqrt(BX**2 + BY**2 + BZ**2))
    pl.show()

if __name__ == "__main__":
    runanalysis()
