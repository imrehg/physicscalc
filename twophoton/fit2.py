from twophotonfit import *

filename = '885 3872 external_101012_174245.csv'
data = loadtxt(filename, delimiter=',')
x = data[:,1]/1e6
y = -data[:,2]
p0t = [48, 0.1, 73, 0.1, 105, 0.1, 145, 0.1, 1, 1, 1, 0.009]
# yhd = multipletripple(p0t, x)
# pl.plot(x, y, '.')
# pl.plot(x, yhd)
yhthree, outthree = completefit(multipletripple, x, y, p0t, filename, 'tripple', toplot=True, tosave=True)

filename = '885 3991 spectrum external_101012_124354.csv'
data = loadtxt(filename, delimiter=',')
x = data[:,1]/1e6
y = -data[:,2]
p0t = [106, 0.1, 145, 0.1, 182, 0.1, 203, 0.1, 1, 1, 1, 0.009]
# yhd = multipletripple(p0t, x)
# pl.plot(x, y, '.')
# pl.plot(x, yhd)
yhthree, outthree = completefit(multipletripple, x, y, p0t, filename, 'tripple', toplot=True, tosave=True)

pl.show()
