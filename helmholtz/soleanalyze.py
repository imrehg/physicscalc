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
    
    pl.subplot(2, 2, 1, aspect='equal')
    pl.title("Pitch = %g, wind number = %d" %(pitch, windnum))
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BX)
    pl.clabel(CS, inline=1, fontsize=10)
    pl.xlabel("Bx")

    pl.subplot(2, 2, 2, aspect='equal')
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BY)
    pl.clabel(CS, inline=1, fontsize=10)
    pl.xlabel("By")

    pl.subplot(2, 2, 3, aspect='equal')
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, BZ)
    pl.clabel(CS, inline=1, fontsize=10)
    pl.xlabel("Bz")

    pl.subplot(2, 2, 4, aspect='equal')
    pl.plot(np.cos(phi), np.sin(phi), linewidth=3)
    CS = pl.contour(X, Y, np.sqrt(BX**2 + BY**2 + BZ**2))
    pl.clabel(CS, inline=1, fontsize=10)
    pl.xlabel("Total field strength")

    pl.show()

if __name__ == "__main__":
    runanalysis()
