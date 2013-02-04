import numpy as np
import pylab as pl

import layeroptimize as lo
import zeemanslower as zs
import sys

from scipy.integrate import odeint
import time


fw, fh = (11.69, 8.27)

def dataplot(filename):
    saveddata = np.load("%s.npz" %filename)
    v, y = saveddata['arr_0'], saveddata['arr_1']

    

filename = "trajectory01"

dataplot(filename)
