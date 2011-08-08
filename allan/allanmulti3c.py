from __future__ import division
import numpy as np
import pylab as pl
import scipy as sp
import ourgui
import time as t
import sys
import allanmulti3



def main(files):
    # names = ["Counter", "Non-overlapping", "Overlapping"]
    # names = ['Rubidium timebase', 'Loran-C timebase']
    # names = ['Tze-wei, Rubidium', 'Greg, Rubidium', 'Greg, Loran-C', 'Repetition rate']
    names = files
    for i, filename in enumerate(files):
        gate, allan, avgf, minf, maxf = np.loadtxt(filename, comments="#", delimiter=",", unpack=True)
        pl.loglog(gate, allan/avgf, '.-', label=names[i], linewidth=2)
    pl.xlabel('Gating time (s)')
    pl.ylabel('Allan deviation (fractional)')
    pl.legend(loc="best")
    pl.show()

if __name__ == "__main__":

    files = []
    filename  = ourgui.openFile(type="log")
    while filename:
        files += [filename]
        filename  = ourgui.openFile(type="log")
    if len(files) > 0:
        main(files)
