#!/usr/bin/python3
#
# Program for normalizing and plotting an observed spectrum.
# Can work interactivly with a graphical interface, through matplotlib.
__version__ = "2.3.1"

#User modifiable parameters
observationName = 'observed.dat'
excludeRegionName = 'exclude.dat'
polynomialsName = 'poly-deg.dat'
paramsName = 'params.dat'

import scipy.constants
c = scipy.constants.c*1e-3 #km/s

import numpy as np
from matplotlib.figure import Figure

#import guiControls2 as guic
import guiMainWin2 as mainWin
import fittingFuncs2 as ff
import time

#Optional input files as command line arguments
import argparse
parser = argparse.ArgumentParser(description='Normalize an observed spectrum, working order by order.  Plots the spectrum and normalizing polynomials, with interactive controls.')
parser.add_argument("observation", nargs='?', default='', help='Observed spectrum file. If none is given, defaults to {:}'.format(observationName))
parser.add_argument("-e", "--exclude",dest='exclude_file',  default='', help='File of regions to exclude from the fit. If no file is given, attempts to use {:}. If the file does not exist the whole spectrum is used.'.format(excludeRegionName))
parser.add_argument("-d", "--degree", dest='degree_file', default='', help='File of degrees of the fitting polynomials. If no file is given, attempts to use {:}. If the file does not exist a default is assumed.'.format(polynomialsName))
parser.add_argument("-c", "--control", dest='control_file', help='File of additional control parameters.  If no file is given, attempts to use {:}. If the file does not exist defaults are used.'.format(paramsName))
parser.add_argument("-b", "--batch", action="store_true", help='Skip interactive plotting and just output a normalized spectrum based on read input parameters.')
args = parser.parse_args()
if args.observation != '':
    observationName = args.observation
if args.exclude_file != '':
    excludeRegionName = args.exclude_file
if args.degree_file != '':
    polynomialsName = args.degree_file
if args.control_file != None:
    paramsName = args.control_file
batchMode = args.batch

#read an observation
nObsCol, inSpec, ords = ff.readObs(observationName, trimMax=100.)
if nObsCol == 2:
    obsWl, obsI, obsSig = inSpec
elif nObsCol == 3:
    obsWl, obsI, obsSig = inSpec
elif nObsCol == 6:
    obsWl, obsI, obsV, obsN1, obsN2, obsSig = inSpec
elif nObsCol == 10:
    obsWl, obsI, obsV, obsN1, obsN2, obsSig = inSpec
elif nObsCol == 30:
    obsWl, obsI, obsSig, obsTel = inSpec

#Read in the exclude from fit wavelength regions, if the file exists
try:
    fExclude = open(excludeRegionName, 'r')
except IOError:
    print('did not find {:}, using defaults'.format(excludeRegionName))
    excludeWls = []
else:
    print('reading exclude ranges from {:}'.format(excludeRegionName))
    excludeWls = []
    i=0
    for line in fExclude:
        if len(line.split()) > 0:
            if line.split()[0][0] != '#':
                excludeWls += [[float(line.split()[0]), float(line.split()[1])]]
                i += 1
    fExclude.close()

#Read in the polynomial degrees from a file, if it exists
try:
    fPolyDeg = open(polynomialsName, 'r')
except IOError:
    print('did not find {:}, using defaults'.format(polynomialsName))
    polyDegs = []
else:
    print('reading polynomial degrees from {:}'.format(polynomialsName))
    polyDegs = []
    i=0
    for line in fPolyDeg:
        if len(line.split()) > 0:
            if line.split()[0][0] != '#':        
                polyDegs += [int(line.split()[0])]
                i += 1
    fPolyDeg.close()

#Initialize a few control parameters for the routine
par = ff.controlPars(averageLen=11, velBin=500., bMergeOrd=False, bFillEdgeGaps=True)
#Read in values for those parameters, if the file exists
par.readParams(paramsName)


#Match up the (optionally) read polynomial degrees with the found spectral orders
polyDegs = ff.fixExtraPolyDegs(polyDegs, ords, defaultDegree=5)

#Moving average (running average) for the intensity spectrum
obsIavg = ff.runningAvg(obsI, ords, par.averageLen)

#Get flag for points that are in fittable regions of the spectrum
bFittable = ff.getIndFittable2(obsWl, obsSig, excludeWls)

#Bin in velocity (km/s) units to search for the best continuum point
fittingWl, fittingI, fittingSig, fittingOrder = ff.getBestInBin(obsWl, obsIavg, obsSig, ords.obsOrder, bFittable, par)

#Fit a polynomial to the selected best continuum points, geting  
fitIvals = ff.fitPoly(obsWl, ords, fittingOrder, fittingWl, fittingI, fittingSig, polyDegs)


#Plotting
if not batchMode:
    #Generate the figure iwth matplotlib, then pass it to the UI generator
    #fig = Figure(figsize=(5, 4), dpi=100) #this seems to cause some odd flikering on updates.
    #Setting figsize here seems to be a minor problem, or at least inefficent.
    fig = Figure()
    ax = fig.add_subplot(111)
    
    #Plot the spectra
    plObs = ax.plot(obsWl, obsI, 'k', linewidth=1.0) #plot observation
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')
    
    #Plot individual orders
    setPlObsO=[]
    setPlPoly=[]
    for i in range(ords.numOrders):
        #Plot fittable regions of individual orders
        obsWlOr = obsWl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1]
        obsIOr = obsI[ords.iOrderStart[i]:ords.iOrderEnd[i]+1]
        obsSigOr = obsSig[ords.iOrderStart[i]:ords.iOrderEnd[i]+1]
        bFittableO = ff.getIndFittable2(obsWlOr, obsSigOr, excludeWls)
        setPlObsO += [ax.plot(obsWlOr[bFittableO], obsIOr[bFittableO])]
        setPlPoly += [ax.plot(obsWl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1], fitIvals[i])]
    #Plot points used for fitting
    plFitting = ax.plot(fittingWl, fittingI, 'ko')

    #Make the main GUI window, with Tkinter, embedding the matplotlib figure.
    #This runs a main loop untill the window is closed.
    mainWin.makeWin(fig, ax, par, polyDegs, ords, obsWl, obsI, obsSig, obsIavg,
                    bFittable, plObs, setPlObsO, setPlPoly, plFitting)
    
    #Update the final continuum polynomial from the y data
    # of the single line plotted, for each order.
    for i in range(ords.numOrders):
        fitIvals[i] =  setPlPoly[i][0].get_data()[1]
    
    print('Done plotting')

#Apply normalization and save output

#Reorder the fit polynomial values to match the 1D observation
fitIvalsFlat = np.zeros(obsI.shape)
i=0
for fitIval in fitIvals:
    ilen = fitIval.shape[0]
    fitIvalsFlat[i:i+ilen] = fitIval
    i += ilen

#Apply the normalization,
#optionally merge spectral orders
#then write the normalized spectra.
#With specific logic depending on the format of the spectrum
merWl = obsWl
merI = obsI = obsI/fitIvalsFlat

if nObsCol == 2:
    if par.bMergeOrd:
        merWl, merI = ff.mergeOrders(ords, obsWl, obsI)
    merWl = ff.scaleWavelength(merWl, par.outputWavelengthScale)
    merWl = ff.convertAirVacuum(merWl, par.outputConvertAirVac)
    ff.writeSpec(observationName+'.norm', nObsCol, merWl, merI)    
if nObsCol == 3:
    merSig = obsSig = np.where(obsSig > 1e9, 10., np.abs(obsSig/fitIvalsFlat))
    if par.bMergeOrd:
        merWl, merI, merSig = ff.mergeOrders(ords, obsWl, obsI, obsSig)
    merWl = ff.scaleWavelength(merWl, par.outputWavelengthScale)
    merWl = ff.convertAirVacuum(merWl, par.outputConvertAirVac)
    ff.writeSpec(observationName+'.norm', nObsCol, merWl, merI, merSig)
if nObsCol == 6:
    merV   = obsV   = obsV/fitIvalsFlat
    merN1  = obsN1  = obsN1/fitIvalsFlat
    merN2  = obsN2  = obsN2/fitIvalsFlat
    merSig = obsSig = np.where(obsSig > 1e9, 10., np.abs(obsSig/fitIvalsFlat))
    if par.bMergeOrd:
        merWl, merI, merV, merN1, merN2, merSig = ff.mergeOrders(ords, obsWl, obsI, obsV, obsN1, obsN2, obsSig)
    merWl = ff.scaleWavelength(merWl, par.outputWavelengthScale)
    merWl = ff.convertAirVacuum(merWl, par.outputConvertAirVac)
    ff.writeSpec(observationName+'.norm', nObsCol, merWl, merI, merV, merN1, merN2, merSig)
if nObsCol == 10:
    merV   = obsV   = obsV/fitIvalsFlat
    merN1  = obsN1  = obsN1/fitIvalsFlat
    merN2  = obsN2  = obsN2/fitIvalsFlat
    merSig = obsSig = np.where(obsSig > 1e9, 10., np.abs(obsSig/fitIvalsFlat))
    if par.bMergeOrd:
        merWl, merI, merV, merN1, merN2, merSig = ff.mergeOrders(ords, obsWl, obsI, obsV, obsN1, obsN2, obsSig)
    merWl = ff.scaleWavelength(merWl, par.outputWavelengthScale)
    merWl = ff.convertAirVacuum(merWl, par.outputConvertAirVac)
    ff.writeSpec(observationName+'.norm', 6, merWl, merI, merV, merN1, merN2, merSig)
#Experimental mode including a telluric spectrum.  The telluric spectrum should already be normalized.
if nObsCol == 30:
    if par.bMergeOrd:
        merWl, merI, merTel = ff.mergeOrders(ords, obsWl, obsI, obsTel)
    merWl = ff.scaleWavelength(merWl, par.outputWavelengthScale)
    merWl = ff.convertAirVacuum(merWl, par.outputConvertAirVac)
    ff.writeSpec(observationName+'.norm', 3, merWl, merI, merTel)

