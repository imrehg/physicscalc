#!/usr/bin/env python2
"""
Plotting the magnetic field after the Zeeman slower
based on previous simulation
"""
from __future__ import division
import numpy as np
import pylab as pl

# Data saved by comparefield2.py
fielddata = np.loadtxt('outside.csv')
pos = fielddata[:, 0]
field = fielddata[:, 1]
pos = pos-pos[0]

pl.figure(figsize=(11.69, 8.27))
pl.plot(pos, field, 'k-', label='Magnetic field simulation', linewidth=3)
pl.plot([0.075, 0.075], [min(field), 0], 'r--', label='Approximate cube centre', linewidth=3)
pl.xlabel('Distance after slower [m]', fontsize=15)
pl.ylabel('Magnetic field [G]', fontsize=15)
pl.grid()
pl.xlim([pos[0], pos[-1]])
pl.ylim([min(field), 0])
pl.legend(loc='lower right')

pl.savefig('fieldoutside.pdf')
pl.savefig('fieldoutside.png')

pl.show()
