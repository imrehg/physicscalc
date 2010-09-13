from twophotonfit import *
from sys import exit


def gruntwork(filename, p0l=None, p0d=None, p0df=None, p0t=None, p0tf=None, justplot=0):

    ##### Start
    data = loadtxt(filename, delimiter=',')
    x = array(range(len(data)))
    y = data[:]
    x = x[:]
    if justplot == 1:
        yhl = multilorentz(p0l, x)
        pl.plot(x, y, '.')
        pl.plot(x, yhl)
        return
    if justplot == 2:
        yhd = multipledouble(p0d, x)
        pl.plot(x, y, '.')
        pl.plot(x, yhd)
        return
    if justplot == 3:
        yhd = multipledouble_freeg(p0df, x)
        pl.plot(x, y, '.')
        pl.plot(x, yhd)
        return
    if justplot == 4:
        yhd = multipletripple(p0t, x)
        pl.plot(x, y, '.')
        pl.plot(x, yhd)
        return
    if justplot == 5:
        yhd = multipletripple_freeg(p0tf, x)
        pl.plot(x, y, '.')
        pl.plot(x, yhd)
        return

    out = []
    if (p0l is not None):
        #### Lorentzian
        yhl, outl = completefit(multilorentz, x, y, p0l, filename, 'lorentzian', toplot=True, tosave=True)
        print "Residual variance:", outl.res_var
        print "##############################"
        out.append(yhl)
    if (p0d is not None):
        #### Transit
        yhtrans, outtrans = completefit(multipledouble, x, y, p0d, filename, 'transit', toplot=True, tosave=True)
        print "Residual variance:", outtrans.res_var
        print "##############################"
        out.append(yhtrans)
    if (p0df is not None):
        #### Transit with free Gamma
        yhtrans, outtrans = completefit(multipledouble_freeg, x, y, p0df, filename, 'transitfree', toplot=True, tosave=True)
        print "Residual variance:", outtrans.res_var
        print "##############################"
    if (p0t is not None):
        #### Three-convolution
        yhthree, outthree = completefit(multipletripple, x, y, p0t, filename, 'tripple', toplot=True, tosave=True)
        print "Residual variance:", outthree.res_var
        print "##############################"
    if (p0tf is not None):
        #### Transit - free Gamma
        yhtrans, outtrans = completefit(multipletripple_freeg, x, y, p0tf, filename, 'transitfree', toplot=True, tosave=True)
        print "Residual variance:", outtrans.res_var
        print "##############################"
    
    saved = zeros((len(x),len(out)+2))
    saved[:, 0] = x[:]
    saved[:, 1] = y[:]
    for i, yh in enumerate(out):
        saved[:, i+2] = yh[:]
    savetxt("%s_fitting.csv"%filename, saved)


#### Params
filename = '883_7222nm.csv'
p0l = [800, 20, 1, 950, 20, 2, 1120, 20, 5, 1308, 20, 8, 1550, 20, 10, 0.2]
p0d = [800, 0.03, 950, 0.06, 1120, 0.15, 1308, 0.24, 1550, 0.40, 16, 14, 0.02]
p0df = [800, 0.03, 16, 950, 0.06, 16, 1120, 0.15, 16, 1308, 0.24, 16, 1550, 0.40, 16, 14, 0.02]
# p0t = [800, 0.03, 950, 0.06, 1120, 0.15, 1308, 0.24, 1550, 0.40, 10, 10, 10, 0.02]
###
# End parameters
###
gruntwork(filename, p0l, p0d, p0df)

#### Params
filename = '883_7102nm.csv'
p0l = [685, 20, 3, 777, 20, 5, 888, 20, 5, 1070, 20, 4.5, 1281, 20, 4, 0.02]
p0d = [685, 0.2, 777, 0.2, 888, 0.2, 1070, 0.15, 1281, 0.10, 16, 14, 0.02]
p0df = [685, 0.2, 16, 777, 0.2, 16, 888, 0.2, 16, 1070, 0.15, 16, 1281, 0.10, 16, 14, 0.02]
# p0t = [685, 0.2, 777, 0.2, 888, 0.2, 1070, 0.15, 1281, 0.10, 10, 10, 10, 0.02]
###
# End parameters
###
gruntwork(filename, p0l, p0d, p0df)


#### Params
filename = '885_3872.csv'
p0l = [854, 20, 5, 1233, 30, 20, 1714, 20, 40, 2300, 20, 90, 0.06]
p0d = [854, 0.2, 1233, 0.8, 1714, 1.7, 2300, 4.4, 16, 14, 0.065]
p0df = [854, 0.2, 16, 1233, 0.8, 16, 1714, 1.7, 16, 2300, 4.4, 16, 14, 0.065]
p0t = [854, 0.2, 1233, 0.8, 1714, 1.7, 2300, 4.4, 10, 10, 14, 0.065]
p0tf = [854, 0.2, 10, 1233, 0.8, 10, 1714, 1.7, 10, 2300, 4.4, 8, 8, 8, 0.065]
###
# End parameters
###
gruntwork(filename, p0l=p0l, p0d=p0d, p0df=p0df, p0t=p0t, p0tf=p0tf, justplot=5)

#### Params
filename = '885_3991.csv'
p0l = [1060, 20, 70, 1244, 20, 60, 1452, 20, 50, 1725, 20, 45, 0.09]
p0d = [1060, 3, 1244, 3, 1452, 2.4, 1725, 2, 9, 10, 0.09]
p0df = [1060, 3, 10, 1244, 3, 10, 1452, 2.4, 10, 1725, 2, 10, 10, 0.09]
p0t = [1060, 3, 1244, 3, 1452, 2.4, 1725, 2, 10, 10, 14, 0.09]
p0tf = [1060, 3, 10, 1244, 3, 10, 1452, 2.4, 10, 1725, 2, 8, 8, 8, 0.09]
###
# End parameters
###
# gruntwork(filename, p0l=p0l)
gruntwork(filename, p0l=p0l, p0d=p0d, p0df=p0df, p0t=p0t, p0tf=p0tf, justplot=0)

# pl.show()
