import numpy as np
from time import strftime
import logging
import fields

nturns = 10
numval = 5
sizemin, sizemax = (0.75, 2.0)
widthmin, widthmax = (0.5, 1.5)
edgemin, edgemax = (0.5, 1)
maxy, numy = (0.5, 10)
maxz, numz = (0.5, 10)


sizes = np.linspace(sizemin, sizemax, numval)
widths = np.linspace(widthmin, widthmax, numval)


# Setup output file
logger = logging.getLogger()
logfile = "helmholtz_%s.log" %(strftime("%y%m%d_%H%M%S"))
hdlr = logging.FileHandler(logfile)
formatter = logging.Formatter('%(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO) 
logger.info("#filename, CoilPosition(inch), CoilWidth(inch), CoilSize(inch), Number_of_turns")


for csize in sizes:
    for cwidth in widths:
        positions = np.linspace(edgemin, edgemax, numval) + cwidth/2.0
        for cpos in positions:
            print (csize, cwidth, cpos)
            hhcoil = fields.HelmholtzCoils(cpos/csize, cwidth/2.0/csize, nturns)
            Y, Z, BY, BZ =  hhcoil.getfield((maxy/csize, numy), (maxz/csize, numz))

            filename = "helmholtz_%s" %(strftime("%y%m%d_%H%M%S"))
            np.savez(filename, Y=Y, Z=Z, BY=BY, BZ=BZ, coilpos=cpos, coilwidth=cwidth, nturns=nturns, coilsize=csize)
            logger.info("%s, %g, %g, %g, %d" %(filename, cpos, cwidth, csize, nturns))


