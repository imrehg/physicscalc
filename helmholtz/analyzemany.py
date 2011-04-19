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

def runanalysis(filename):

    data = np.load(filename, mmap_mode='r')
    Y = data['Y']
    Z = data['Z']
    BY = data['BY']
    BZ = data['BZ']
    coilpos = getdata(data, 'coilpos')
    coilwidth = getdata(data, 'coilwidth')
    nturns = getdata(data, 'nturns')
    coilsize = getdata(data, 'coilsize')
    vary = str(getdata(data, 'varied'))

    infopiece = []
    if coilpos:
        infopiece += ['Cpos: %g"' % coilpos]
    if coilwidth:
        infopiece += ['Cwidth: %g"' % coilwidth]
    if coilsize:
        infopiece += ['Csize: %g"' % coilsize]
        Y *= coilsize
        Z *= coilsize
    if nturns:
        infopiece += ["Turns: %d" % nturns]
    infotitle = ", ".join(infopiece)

    fig = pl.figure(1, figsize=(11.69, 8.27), dpi=100)
    # pl.subplot(2,2,3)
    pl.title("Bz-field as function of position and %s" %(vary))
    zpos = Z[:, 0]
    zfield = BZ[:,0]/BZ[0,0]
    pl.plot(zpos, zfield, label="%g" %(getdata(data, vary)))
    pl.xlabel('Z position')
    pl.ylabel('On axis field strength')
    pl.legend(loc='best')
    return (max(abs(zfield - 1)), getdata(data, vary))

if __name__ == "__main__":
    logname = ourgui.openFile(type="log")
    logged = np.loadtxt(logname,
                        delimiter=",",
                        comments="#",
                        dtype={'names': ('filename', 'coilpos', 'coilwidth', 'coilsize', 'nturns'),
                               'formats': ('a40', 'f4', 'f4', 'f4', 'i4')}
                        )
    
    best = 1e14
    bestval = None
    for log in logged:
        maxdev, maxval = runanalysis("%s.npz" %log['filename'])
        if maxdev < best:
            best = maxdev
            bestval = maxval
            bestsetting = log['filename']
    print "Best setting: %s (%f for %f)" %(bestsetting, best, bestval)
    pl.show()
