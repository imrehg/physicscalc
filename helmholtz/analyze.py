#!/usr/bin/python2
import pylab as pl
import numpy as np
import ourgui

def getdata(data, field):
    try:
        this = data[field]
    except:
        this = None
    return this

def runanalysis():
    filename = ourgui.openFile(type='npz')

    data = np.load(filename, mmap_mode='r')
    Y = data['Y']
    Z = data['Z']
    BY = data['BY']
    BZ = data['BZ']
    coilpos = getdata(data, 'coilpos')
    coilwidth = getdata(data, 'coilwidth')
    nturns = getdata(data, 'nturns')
    coilsize = getdata(data, 'coilsize')

    infopiece = []
    if coilpos:
        infopiece += ['Cpos: %g"' % coilpos]
    if coilwidth:
        infopiece += ['Cwidth: %g"' % coilwidth]
    if coilsize:
        infopiece += ['Csize: %g"' % coilsize]
    if nturns:
        infopiece += ["Turns: %d" % nturns]
    infotitle = ", ".join(infopiece)


    fig = pl.figure(figsize=(11.69, 8.27), dpi=100)
    fig.text(.4, .95, infotitle)

    pl.subplot(2,2,1)
    pl.quiver(Z, Y, BZ, BY)
    pl.xlabel('Z')
    pl.ylabel('Y')
    pl.title('Magnetic field direction')

    pl.subplot(2,2,2)
    CS = pl.contour(Z, Y, BY/BZ)
    pl.xlabel('Z')
    pl.ylabel('Y')
    pl.title('Y-strength/Z-strength')
    pl.clabel(CS, inline=1, fontsize=10)

    pl.subplot(2,2,3)
    zpos = Z[:, 0]
    zfield = BZ[:,0]/BZ[0,0]
    pl.plot(zpos, zfield)
    pl.xlabel('Z position')
    pl.ylabel('On axis field strength')

    pl.subplot(2,2,4)
    fieldstrength = np.sqrt(BY**2 + BZ**2)
    CS = pl.contour(Z, Y, fieldstrength)
    pl.xlabel('Z')
    pl.ylabel('Y')
    pl.title('Field strength', fontsize=10)
    pl.clabel(CS, inline=1, fontsize=10)

    pl.show()

if __name__ == "__main__":
    runanalysis()
