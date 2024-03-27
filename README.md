# normPlot
A tool for normalizing stellar spectra, with an interactive GUI.


## Overview

This program is designed for continuum normalizing reduced 1D stellar spectra, typically from an echelle spectrograph.  

The program is intended to work on one echelle spectral order at a time.  The algorithm was built for spectra with heavy blending, so the goal is to identify a relatively small number of 'good' continuum points, then fit a low degree polynomial through them.  The algorithm breaks each spectral order up into a set of bins (of user defined size in velocity units), and looks for the best continuum point to use in each bin.  The 'best' continuum point is found by calculating a running/sliding average of the spectrum (with a user defined size in pixels), then taking the highest averaged point in the bin as the most likely continuum.  These points are used for fitting the continuum polynomial and are shown in the plot as big black points.  Once a satisfactory set of continuum polynomials have been found, the original spectrum is divided by those polynomials to produce the normalized spectrum.

When using the program, typically one would set the running average large enough to reduce the impact of noise.  Then set the bin size large enough that there is generally some adequate continuum inside the bin.  Then emission lines or regions with particularly broad or blended lines can be dealt with by by excluding them from the fit by hand.  A polynomial degree of 5 or 6 is often sufficient, but for short or very noisy spectral orders a lower degree may be better, while large spectral orders with more complex behaviour may benefit from a higher degree.

The program has an optional feature for dealing with spectral orders that start/end in a broad absorption feature (e.g. a Balmer line).  If the spectrum around the feature is excluded from the fit, and the spectral order ends in this excluded region, the algorithm can try to take a point from the previous/next valid included region (in the previous/next order) and use that to constrain the polynomial.  While this is not particularly accurate, it can help keep the polynomial well behaved.

## Installing

normPlot can be downloaded and installed from PyPI with:
```
pip install normPlot
```
Or you can download the source files, and some example input files, from Github at [https://github.com/folsomcp/normPlot](https://github.com/folsomcp/normPlot)


## Using

If you have installed the program using pip you can run normPlot with
```
normplot [observation_file]
```
You can also use `normplot -h` for information about additional command line parameters.  You can also use normPlot inside other scripts like
```
import normPlot
normPlot.normplot('[observation_filename]')
```
see `help(normPlot.normplot)` for a full list of function arguments.

If you have downloaded the source files from Github, the main executable code is `normPlot2.py` in the`normPlot/` folder.  The program is in Python 3 and uses tkinter for the GUI and matplotlib for visualizing spectra.  You can run the program on an observation with
```
python3 normPlot/normPlot2.py [observation_file]
```
and you can run `python3 normPlot/normPlot2.py -h` for some other command line parameters.

When run, program will try to read additional information from exclude.dat, poly-deg.dat, and params.dat, if they exist (if not, defaults will be assumed).  Alternate files can be specified using command line parameters.  When the program is running these parameters can be modified interactively, and the chosen parameters can be saved to those files using "save params" button.  The program can be ran without the interactive window, just reading saved parameter from files, with the -b option. 

The observation file should be in a text format with columns of wavelength, intensity, and optionally errors or polarimetric information and errors.  Spectral orders are identified by looking for overlap or gaps in wavelength.  Spectra with gaps in wavelength for other reasons (e.g. a large number of omitted bad pixels) may cause issues for the order identification.

The program has an interactive UI for defining regions to include/exclude from the fit, simply by clicking on a plot of the observation.  It also plots the selected 'good' continuum points and the fit continuum polynomial.  When the interactive window is closed, the normalized spectrum is saved as [observation_file].norm

The "set poly. degree..." button can set the degree of the polynomial used for each spectral order.  The "set output params..." button provides a few options. You can chose to merge spectral orders in a simple fashion (this just cuts the orders off at the midpoint of the overlapping region and then splices them together).  You can convert the wavelength with "Scale output wavelength by" (e.g. convert from nm to Angstroms by setting it to 10).  You can also convert wavelengths between wavelength in vacuum and in air.  The calibration requires wavelength in A, and is applied after the scaling from "Scale output wavelength by", so set that scaling appropriately to convert to A. 

You should be able to pan the plotted spectra with arrow keys and zoom to a selected region with the 'zoom' button or the z key.  The "include range" and "exclude range" allow interactively selecting regions of the spectrum to include in the fit by clicking on the plot.  The "fit cont." button updates the fit of the continuum polynomial.  The "fill order edge gaps" check box enables the feature where, if one spectral order ends in an exclude region, the algorithm can take a point in the adjacent spectra order and use that to constrain the continuum polynomial (as mentioned above).  Closing the window will then save the normalized spectrum to [observation_file].norm
