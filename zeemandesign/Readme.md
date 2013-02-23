# Zeemandesign

Zeeman slower design software

## Files

The list of files and their roles

* `getlayers.py`: display and list the calculated layers and fields (run after layerstudy)
* `layeroptimize.py`: the main program for calculation
* `layerstudy.py`: in its current form, the interface to set up optimization for the best Zeeman slower arrangement, using layeroptimize
* `wires.py`: collection of different wires with their properties (radius, resitance)
* `zeemanslower.py`: general slower calculation for slower length, ideal field, the different atomic species

## Required libraries

* [Numpy](http://www.numpy.org/)
* [Scipy](http://www.scipy.org/)
* [Matplotlib (pylab)](http://matplotlib.org/)
* [WxPython](http://www.wxpython.org/) (for ourgui)
