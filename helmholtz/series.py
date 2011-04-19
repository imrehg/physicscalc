import numpy as np
from time import strftime
import logging
import fields

choice = ['', 'coilsize', 'coilwidth', 'coilpos', 'nturns']
series = None
while series not in range(1, len(choice)):
    try:
        series = int(raw_input("What to vary: 1) radius, 2) width, 3) position, 4) number of turns "))
    except:
        pass

vmin = float(raw_input("Minimum: "))
vmax = float(raw_input("Maximum: "))
steps = int(raw_input("Number of steps: "))
varied = np.linspace(vmin, vmax, steps)

params = {}

for ch in range(1, len(choice)):
    if ch <> series:
        params[choice[ch]] = float(raw_input("%s: " %(choice[ch])))

### Keep this the same for the moment:
maxz, numz, maxy, numy = (0.5, 10, 0.5, 10)
# maxz = float(raw_input("Z-max: "))
# numz = int(raw_input("Z-steps: "))
# maxy = float(raw_input("Y-max: "))
# numy = int(raw_input("Y-steps: "))


# Setup output file
logger = logging.getLogger()
logfile = "helmholtz_%s.log" %(strftime("%y%m%d_%H%M%S"))
hdlr = logging.FileHandler(logfile)
formatter = logging.Formatter('%(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO) 
logger.info("#filename, CoilPosition(inch), CoilWidth(inch), CoilSize(inch), Number_of_turns")

for val in range(steps):
    params[choice[series]] = varied[val]
    params['nturns'] = int(params['nturns'])
    print params
    hhcoil = fields.HelmholtzCoils(params['coilpos']/params['coilsize'],
                                   params['coilwidth']/(2.0*params['coilsize']),
                                   params['nturns'],
                                   )
    Y, Z, BY, BZ =  hhcoil.getfield((maxy/params['coilsize'], numy), (maxz/params['coilsize'], numz))
    filename = "helmholtz_%s" %(strftime("%y%m%d_%H%M%S"))
    np.savez(filename, Y=Y, Z=Z, BY=BY, BZ=BZ,
             coilpos=params['coilpos'],
             coilwidth=params['coilwidth'],
             nturns=params['nturns'],
             coilsize=params['coilsize'],
             varied=choice[series],
             )
    logger.info("%s, %g, %g, %g, %d" %(filename, params['coilpos'], params['coilwidth'], params['coilsize'], params['nturns']))
