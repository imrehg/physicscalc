#!/usr/bin/env python
"""
Fitting data from 2011-05-26

Aiming to make it more general fitting program later
"""
import numpy as np
import pylab as pl

import twophotonfit as tw

filename = "gigahertz_110525_235920.log.npz"
data = np.load(filename)

# pl.plot(data['favgs'], data['rchan'])
# pl.show()

x = data['favgs']
y = data['rchan']
yerr = np.sqrt(data['xchanerr']**2 + data['ychanerr']**2)
p0t = [-113, 0.01, -64, 0.01, 0, 0.05, 82, 0.09, 1, 1, 1, 0.001]

# # yh = tw.multipletripple(p0t, x)

# # pl.subplot(211)
# # pl.plot(x, y, '.')
# # pl.plot(x, yh, '-')

# # pl.subplot(212)
# # pl.plot(x, y-yh, '.')

yh, outthree = tw.completefit(tw.multipletripple, x, y, p0t, filename, 'tripple', toplot=True, tosave=True)

pl.subplot(211)
pl.plot(x, y, '.')
pl.plot(x, yh, '-')

pl.subplot(212)
pl.plot(x, y-yh, '.')

pl.savetxt(filename+".csv", zip(x, y, yerr, yh), delimiter=',')
pl.show()
