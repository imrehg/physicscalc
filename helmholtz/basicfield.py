import numpy as np
import pylab as pl

def magfield(R, x):
    return R**2/(2*(R**2 + x**2)**1.5)


def plotfield(radius, offset, pos, colour="k", details=True, name=None):
    fieldright = magfield(radius, pos - offset)
    fieldleft = magfield(radius, pos + offset)
    totalfield = fieldleft+fieldright
    middle = int((len(pos)-1)/2.0 + 1)
    scale = totalfield[middle]
    if  details:
        pl.plot(pos, fieldright/scale, '%s' %colour)
        pl.plot(pos, fieldleft/scale, '%s' %colour)
    pl.plot(pos, totalfield/scale, '%s' %colour, linewidth=3, label=name)

rad = 1
zpos = 1
pos = np.linspace(-zpos, zpos, 201)
details = True
# details = False

pl.figure(figsize=(8.27, 11.69), dpi=100)

pl.subplot(2, 1, 1)
series = [(1, "r"), (2.5, "k"), (10, "y")]
for this in series:
    plotfield(rad*this[0], zpos, pos, this[1], details, "R=%g" %this[0])
pl.legend(loc='best')

pl.subplot(2, 1, 2)
series = [(1.7, "y--"), (1.8, "k:"), (1.9, "r-"), (2.0, "b-."), (2.1, "g.")]
for this in series:
    plotfield(this[0], zpos, pos, this[1], False, "R=%g" %this[0])
pl.legend(loc='best')

# pl.savefig("test.pdf")
pl.show()

