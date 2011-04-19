### Information about the functionality in these files

This should be included in the .py files but not always sure how

#### General:

Every internal value is scaled by the coil radius. The output values
are in this format, this if in the output there is y = 1.0, that is
the distance is equal to the radius (the real value saved in the 
'coilsize' parameter).

#### fields.py

General field simulation of one set of parameters

#### series.py

Run a series of simulation varying either of:

 * coil size (radius)
 * coil width
 * coil centre position
 * number of turns

Output: a single log file wich holds the filenames and parameters of
this sequence, plus all the saved files. The saved files hold the mesh
locations (Z & Y) and fields (BZ(z,y), BY(z,y)) as well as the coil
parameters and the name of the parameter being varied in the series.

#### analyzemany.py

Analyize this series for Bz(z, 0) by checking the .log file of the series.
Outputs plot and parameter values of the solution with smallest deviation from
the centre value

#### analyze.py

Analyze a single file, getting field direction, By/Bz, Bz(z,0)/Bz(0,0), and Bz.
