#Helper functions for continuum normalization fitting

import numpy as np
import time
import scipy.constants
c = scipy.constants.c*1e-3 #km/s

#Save a few basic control parameters for the routine,
#Use a little object for this, to make passing to other functions easier.
class controlPars():
    def __init__(self, averageLen=1, velBin=1, bMergeOrd=False, bFillEdgeGaps=False, outWaveScale=1.0, outConvertAirVac=0):
        #length of the running average for smoothing the continuum
        self.averageLen = int(averageLen) 
        #velocity bin (km/s) to search for a continuum point
        self.velBin = velBin  
        #Flag for merging spectral orders
        self.bMergeOrd = bMergeOrd  
        #Flag for, if a order ends in an exclude region,
        #take an extra good point from the next/last order
        self.lookToNextOrderForGaps = bFillEdgeGaps 
        #Scale the output wavelength by some value, for converting units
        self.outputWavelengthScale = outWaveScale
        #Flag for converting wavelength between air and vacuum
        #(-1 = air to vacuum, 0 = no conversion, 1 = vacuum to air)
        self.outputConvertAirVac = outConvertAirVac

    #Set control parameters from an input file (if the file exists)
    def readParams(self, fileName):
        try:
            fParams = open(fileName, 'r')
        except IOError:
            #if the file dosen't exist just use the defaults above
            print('did not find {:}, using defaults'.format(fileName))
            return
        else:
            print('reading control parameters from {:}'.format(fileName))
            i = 0
            for line in fParams:
                if len(line.split()) > 0:
                    if line.split()[0][0] != '#':
                        if i == 0:
                            self.averageLen = int(line.split()[0])
                        elif i == 1:
                            self.velBin = float(line.split()[0])
                        elif i == 2:
                            if int(line.split()[0]) == 0:
                                self.bMergeOrd = False
                            else:
                                self.bMergeOrd = True
                        elif i == 3:
                            if int(line.split()[0]) == 0:
                                self.lookToNextOrderForGaps = False
                            else:
                                self.lookToNextOrderForGaps = True
                        elif i == 4:
                            waveScale = float(line.split()[0])
                            if waveScale > 0.0:
                                self.outputWavelengthScale = waveScale
                        elif i == 5:
                            convAirVacuum = int(line.split()[0])
                            if (-1 <= convAirVacuum) & (convAirVacuum <= 1):
                                self.outputConvertAirVac = convAirVacuum
                        i += 1

    #Helper function for using a button, link this object to some usefull data
    def linkPolyDegExclude(self, polyDegs, bFittable, obsWl, ords):
        self.polyDegs = polyDegs
        self.bFittable = bFittable
        self.obsWl = obsWl
        self.ords = ords

    #Save the parameters used for the fit & normalization, on a button press 
    def saveValues(self, **event):
        fname = 'params.dat'
        fnamePoly = 'poly-deg.dat'
        fnameExclude = 'exclude.dat'

        print('saving general parameters to {:}'.format(fname))
        fout = open(fname, 'w')
        fout.write('#length of average (pixels)\n')
        fout.write('{:n}\n'.format(self.averageLen))
        fout.write('#size of search bin (km/s)\n')
        fout.write('{:.2f}\n'.format(self.velBin))
        fout.write('#merge spectral orders (1=yes, 0=no)\n')
        fout.write('{:n}\n'.format(self.bMergeOrd))
        fout.write('#fill gaps at order edges with points from nearby orders\n')
        fout.write('{:n}\n'.format(self.lookToNextOrderForGaps))
        fout.write('#scale output wavelength by a factor (e.g. for unit conversion)\n')
        fout.write('{:.3f}\n'.format(self.outputWavelengthScale))
        fout.write('#convert wavelengths in air to vacuum (-1), no conversion (0), or vacuum to air (+1)\n')
        fout.write('#(this requires units in A, so scale wavelength appropriately first)\n')
        fout.write('{:n}\n'.format(self.outputConvertAirVac))
        fout.close()

        print('saving polynomial degrees to {:}'.format(fnamePoly))
        fpoly = open(fnamePoly, 'w')
        fpoly.write('#Polynomial degree, for each spectral order\n')
        i = 1
        for polyDeg in self.polyDegs:
            fpoly.write('{:n}  {:n}\n'.format(polyDeg, i))
            i += 1
        fpoly.close()

        #Infer exclude regions (wavelength ranges) from the array of fittable points bFittable
        #First strip out order overlap, since that would just add extraneous exclude regions
        merWl, merFit = mergeOrders(self.ords, self.obsWl, self.bFittable)
        merFit = merFit.astype(int)
        
        #Find where bFittable flags change, for exclude range edges
        edges =  np.nonzero(merFit[1:-1] != merFit[0:-2])[0] + 1
        #Use edges where the bFittable flag has gone false for a exclude region start,
        #and back to true for ends.
        starts =  edges[merFit[edges] == 0] 
        ends =  edges[merFit[edges] == 1]
        #Protect against exclude regions at the start or end of the spectrum
        if ends[0] < starts[0]:
            starts = np.insert(starts, 0, 0)
        if starts[-1] > ends[-1]:
            ends = np.append(ends, merFit.shape[0]-1)
        if starts.shape[0] != ends.shape[0]:
            print('ERROR calculating exclude regions, not writting to file')
            return

        print('saving exclude regions to {:}'.format(fnameExclude))
        fexclude = open(fnameExclude, 'w')
        for i in range(starts.shape[0]):
            fexclude.write('{:.2f} {:.2f}\n'.format(merWl[starts[i]], merWl[ends[i]]))
        fexclude.close()
        return


#read an observation
def readObs(observationName, trimMax=-1.):
    #check number of columns
    fObs = open(observationName, 'r')
    fObs.readline() #skip possible headder
    fObs.readline()
    checkLine = fObs.readline()
    nObsCol = len(checkLine.split())
    fObs.close()
    #assume observation format based on number of columns
    if nObsCol == 2:
        print('reading {:}, 2 column spectrum, assuming constant uncertainty'.format(observationName))
        obsWl, obsI = np.loadtxt(observationName, usecols = (0,1), skiprows=2, unpack=True)
        obsSig = np.std(obsI)*np.ones(obsI.shape)
        inSpec = [obsWl, obsI, obsSig]
    elif nObsCol == 3:
        print('reading {:}, 3 column spectrum, assuming input uncertainties'.format(observationName))
        obsWl, obsI, obsSig = np.loadtxt(observationName, usecols = (0,1,2), skiprows=2, unpack=True)
        inSpec = [obsWl, obsI, obsSig]
        #if np.mean(obsSig) > 0.1*np.mean(obsI):
        #    print('experimental: treating 3rd colmun as SPIRou telluric spectrum')
        #    obsWl, obsI, obsTel = np.loadtxt(observationName, usecols = (0,1,2), skiprows=2, unpack=True)
        #    obsSig = np.std(obsI)*np.ones(obsI.shape)
        #    inSpec = [obsWl, obsI, obsSig, obsTel]
        #    nObsCol = 30
    elif nObsCol == 6:
        print('reading {:}, 6 column spectrum, assuming ESPaDOnS format'.format(observationName))
        obsWl, obsI, obsV, obsN1, obsN2, obsSig = np.loadtxt(observationName, usecols = (0,1,2,3,4,5), skiprows=2, unpack=True)
        inSpec = [obsWl, obsI, obsV, obsN1, obsN2, obsSig]
    elif nObsCol == 10:
        print('experimental: reading {:}, 10 column spectrum, assuming Donati SPIRou format'.format(observationName))
        obsWl, obsI, obsV, obsN1, obsN2, obsSig = np.loadtxt(observationName, usecols = (0,2,3,4,5,7), skiprows=2, unpack=True)
        inSpec = [obsWl, obsI, obsV, obsN1, obsN2, obsSig]
    elif nObsCol == 7:
        print('experimental: reading {:}, 7 column spectrum, from SPIRou DRS "p" (skipping second error column)'.format(observationName))
        obsWl, obsI, obsV, obsN1, obsN2, obsSig = np.loadtxt(observationName, usecols = (0,1,2,3,4,5), skiprows=2, unpack=True)
        nObsCol = 6
        inSpec = [obsWl, obsI, obsV, obsN1, obsN2, obsSig]
    else:
        print('ERROR: found an unexpected number of columns ({:}) in observation: {:}'.format(nObsCol,observationName))
        import sys
        sys.exit()

    orders = orderEdges(obsWl)
    
    #Weak protection against small/bad values in I or sigma
    meanIapprox = np.mean(obsI)
    indUse = (obsI > 1e-10*meanIapprox) & (obsSig > 1e-10*meanIapprox)
    for i in range(1,len(inSpec)):
        inSpec[i][np.logical_not(indUse)] = 0.0
    obsSig[np.logical_not(indUse)] = 1e10
    
    #Optionally protect against very large values too
    if trimMax > 0.:
        meanIapprox = np.mean(obsI)
        indUse = (inSpec[1] < trimMax*meanIapprox)
        for i in range(1,len(inSpec)):
            inSpec[i][np.logical_not(indUse)] = 0.0
        obsSig[np.logical_not(indUse)] = 1e10
    

    return nObsCol, inSpec, orders



#Find spectral order edges
#Assume order edges occur where there is overlap in the 1D spectrum
#(i.e. wavelength goes backwards), or where there are unusualy large gaps
# in the obseration (relative to the mean pixel size in velocity).
class orderEdges:
    def __init__(self, obsWl):
        #Define spectral order edges by a step backwards in wavelength (velocity)
        #or a step forward in velocity more than 10x the average velocity step size.
        #(use velocity rather than wavelength since pixel size in velocity is
        # more consistent across a spectrum)
        velSteps = (obsWl[1:-1] - obsWl[0:-2])/obsWl[1:-1]*c
        meanVelStep = np.mean(velSteps)
        orderEdges = np.logical_or(velSteps < 0., velSteps > 20.*meanVelStep)
        indOrderEdges =  np.nonzero(orderEdges) #last point in a spectral order
        self.numOrders = indOrderEdges[0].shape[0]+1
        print('found {:} spectral orders:'.format(self.numOrders))
        #Get order start and end array indexes for convenience
        self.iOrderStart = np.zeros(self.numOrders, dtype=int)
        self.iOrderEnd = np.zeros(self.numOrders, dtype=int)
        self.obsOrder = np.zeros(obsWl.shape, dtype=int)
        for i in range(self.numOrders):
            if i == 0:
                self.iOrderStart[i]=0
            else:
                self.iOrderStart[i] = indOrderEdges[0][i-1]+1
            if i == self.numOrders-1:
                self.iOrderEnd[i] = obsWl.shape[0]-1
            else:
                self.iOrderEnd[i] = indOrderEdges[0][i]
            self.obsOrder[self.iOrderStart[i]:self.iOrderEnd[i]+1] = i
            print('{:.4f}  {:.4f}'.format(obsWl[self.iOrderStart[i]], obsWl[self.iOrderEnd[i]]) )
        self.wlOrderStart = obsWl[self.iOrderStart]
        self.wlOrderEnd = obsWl[self.iOrderEnd]

    #When trim out points from the spectrum (e.g. bad pixels), update the spectral orders
    #Takes an array of booleans for usable points.
    #This needs to update the array of orders for each pixel,
    #and the index positions of order edges
    def trimPoints(self, indUse):
        for i in range(self.numOrders):
            #Get the position of the usable pixel closest to the order start
            #More accurately find the last useable pixel before the order start,
            #then take the next useable pixel after that.
            #(since [0:start] dosen't include start)
            #(and using shape of an array as an index gives the ending index +1)
            istart = np.nonzero(indUse[:self.iOrderStart[i]])[0].shape[0]
            #get the position of last usable pixel before/at the end of an order
            iend = np.nonzero(indUse[:self.iOrderEnd[i]+1])[0].shape[0]-1
            self.iOrderStart[i] = istart
            self.iOrderEnd[i] = iend
            
            ##Alternate logic, conceptualy simpler, but more complicated code
            ##get an array of positions for the old array
            #inddic = np.nonzero(indUse)[0]
            ##get the last usable pixel before/at the end of an order
            #indEndUse = np.nonzero(indUse[:self.iOrderEnd[i]+1])[0][-1]
            ##get the first usable pixel after/at the start of an order
            #indStartUse = np.nonzero(indUse[self.iOrderStart[i]:])[0][0] + self.iOrderStart[i]
            ##Get the position in the new array corrisponding to the index in the old array
            #istart2 = np.nonzero(inddic == indStartUse)[0][0]
            #iend2 = np.nonzero(inddic == indEndUse)[0][0]
            
        self.obsOrder = self.obsOrder[indUse]
        return

#Match up the (optionally) read polynomial degrees with the found spectral orders
def fixExtraPolyDegs(polyDegs, ords, defaultDegree=5):
    #If there were input polynomial degrees
    if len(polyDegs) > 0:
        #cut extra polynomial degree values, or padd out with the default value
        if len(polyDegs) > ords.numOrders:
            print('Read too many polynomial degrees ({:}) for the number of spectral orders ({:}), ignoring the extras'.format(len(polyDegs), ords.numOrders))
            polyDegs = polyDegs[0:ords.numOrders]
        elif len(polyDegs) < ords.numOrders:
            print('Read too few polynomial degrees ({:}) for the number of spectral orders ({:}), padding with default {:}'.format(len(polyDegs), ords.numOrders, defaultDegree))
            polyDegs += [defaultDegree]*(ords.numOrders - len(polyDegs))
    else:
        #If we have no input information just set all polynomial degrees to some default.
        polyDegs = [defaultDegree]*ords.numOrders  #make a list ords.numOrders long
    return polyDegs

#Get a boolean array of pixels that are not in exclude regions
def getIndFittable2(obsWl, obsSig, excludeWls):
    bFittable = np.ones(obsWl.shape, dtype=bool)
    for exclude in excludeWls:
        bFittable &= (obsWl < exclude[0]) | (obsWl > exclude[1])
    indBadPix = obsSig > 1e9
    bFittable[indBadPix] = False
    return bFittable

#Modify a boolean array of pixels that are not in exclude regions
def cutIndFittable(bFittable, obsWl, excludeWls):
    for exclude in excludeWls:
        bFittable &= (obsWl < exclude[0]) | (obsWl > exclude[1])
    return bFittable

#Modify a boolean array of pixels that are not in exclude regions
def addIndFittable(bFittable, obsWl, includeWls):
    for include in includeWls:
        bFittable |= (obsWl > include[0]) & (obsWl < include[1])
    return bFittable


#Moving average (Running average) for the intensity spectrum
def runningAvg(obsI, ords, averageLen):
    obsIavg = np.zeros(obsI.shape)
    #Compute the averages only within one spectral order.
    for i in range(ords.numOrders):
        i1 = ords.iOrderStart[i]
        i2 = ords.iOrderEnd[i]+1
        
        #Fast but more confusing moving average: use a cumulative sum to do the summation.
        # Then take the difference in cumsum between start and end points of
        # this bin of the running average to get the sum across this bin.
        #Padd using the starting and ending values
        #Note: even values of averageLen will produce a 1/2 pixel wavelengh error
        obsIrange = np.insert(obsI[i1:i2], 0, obsI[i1]*np.ones(int((averageLen/2))+1))
        obsIrange = np.append(obsIrange, obsI[i2-1]*np.ones(int(np.ceil(averageLen/2.))-1))
        cumsum = np.cumsum(obsIrange)
        obsIavg[i1:i2] =  (cumsum[averageLen:] - cumsum[:-averageLen])/float(averageLen)
        #cumsum = np.cumsum(np.insert(obsI[i1:i2], 0, 0))  #for no padding
        #obsIavgT[i1:i2] =  (cumsum[averageLen:] - cumsum[:-averageLen])/float(averageLen)

    return obsIavg


#Look ahead to the next order, or back to the previous order, to find a good point for fitting.
#A good point is the maximum point in a bin, and not in an exclude region.
#The direction of the search is inferred based on whether the order for this bin
#matches the order for the previous/next bins.
def lookToNextOrderForGaps(obsWl, obsIavg, obsSig, obsOrder, bFittable, par, 
                           orderReplace, iBinStart, iBinNext):
    fittingWl = -1.
    fittingI = 0.
    fittingSig = 0.
    #if we are at the start of an order
    #(the order number does not match the previous bin's order number)
    if (orderReplace != obsOrder[iBinStart-1]) & (iBinNext > 0):
        #Run backwards looking for a good bin with a good point
        binWlEndB = obsWl[iBinStart]
        binWlStartB = binWlEndB - par.velBin/c*binWlEndB
        iBinEndB = iBinStart
        for j in range(iBinStart-1,0,-1):
            #If we moved to the wrong part of an overlap,
            #update the end index for the search bin, so we search
            #a continuous (in wl) set of pixels from binWlStartB to binWlEndB
            if (obsWl[j] > binWlEndB):
                iBinEndB = j
            #use points from any spectral order (hopefully an adjacent one)
            if (obsWl[j] < binWlStartB):
                #protect against too small bins
                if iBinEndB > j+1:
                    indMax = j+1 + np.argmax(obsIavg[j+1:iBinEndB])
                else:
                    indMax = j+1
                #If this point is not in an exclude region,
                # save it and exit the loop
                if bFittable[indMax] == True:
                    fittingWl = obsWl[indMax]
                    fittingI = obsIavg[indMax]
                    fittingSig = obsSig[indMax]
                    return fittingWl, fittingI, fittingSig
                binWlEndB = obsWl[j]
                binWlStartB = binWlEndB - par.velBin/c*binWlEndB
                iBinEndB = j
                
    #if we are at the end of an order
    #(the order number does not match the next bin's order number)
    elif orderReplace != obsOrder[iBinNext]:
        #Run ahead looking for the next good bin with a good point
        binWlStartB = obsWl[iBinNext-1]
        binWlEndB = binWlStartB + par.velBin/c*binWlStartB
        iBinStartB = iBinNext
        for j in range(iBinNext, obsWl.shape[0]):
            #If we are in the wrong part of an overlap,
            #update the start index for the search bin, so we search
            #a continuous (in wl) set of pixels from binWlStartB to binWlEndB
            if (obsWl[j] < binWlStartB):
                iBinStartB = j
            #If we hit a bin end for any spectral order
            if (obsWl[j] > binWlEndB):
                #protect against too small bins
                if iBinStartB < j:
                    indMax = np.argmax(obsIavg[iBinStartB:j])+iBinStartB
                else:
                    indMax = j-1
                #If this point is not in an exclude region,
                # save it and exit the loop
                if bFittable[indMax] == True:
                    fittingWl = obsWl[indMax]
                    fittingI = obsIavg[indMax]
                    fittingSig = obsSig[indMax]
                    return fittingWl, fittingI, fittingSig
                binWlStartB = obsWl[j]
                binWlEndB = binWlStartB + par.velBin/c*binWlStartB
                iBinStartB = j

    return fittingWl, fittingI, fittingSig


#Get the highest ('most continuum') point in each bin
#The exclude points that are inside an exclude region
def getBestInBin(obsWl, obsIavg, obsSig, obsOrder, bFittable, par):
    fittingWl = np.zeros(obsWl.shape)
    fittingI = np.zeros(obsWl.shape)
    fittingSig = np.zeros(obsWl.shape)
    fittingOrder = np.zeros(obsWl.shape, dtype=int)
    binWlStart = obsWl[0]
    binWlEnd = binWlStart + par.velBin/c*binWlStart
    iBinStart=0
    nFitPts = 0
    for i in range(obsWl.shape[0]):
        #If we hit a bin end or a spectral order end
        if (obsWl[i] > binWlEnd) | (obsOrder[i] != obsOrder[iBinStart]):
            #protect against too small bins
            if iBinStart < i:
                indMax = np.argmax(obsIavg[iBinStart:i])+iBinStart
            else:
                indMax = i-1
            
            #If this point is not in an exclude region
            if bFittable[indMax] == True:
                fittingWl[nFitPts] = obsWl[indMax]
                fittingI[nFitPts] = obsIavg[indMax]
                fittingSig[nFitPts] = obsSig[indMax]
                fittingOrder[nFitPts] = obsOrder[indMax]
                nFitPts += 1
            else:
                #if excluded a point, add an filler point at the start/end of a spectral order
                if par.lookToNextOrderForGaps:
                    iBinNext = i
                    orderReplace = obsOrder[indMax]
                    #run the function for all excluded points,
                    #since the function needs logic for figuring out order ends anyway,
                    #in order to look forwards/backwards properly. 
                    tmpWl, tmpI, tmpSig = lookToNextOrderForGaps(
                        obsWl, obsIavg, obsSig, obsOrder, bFittable, par, 
                        orderReplace, iBinStart, iBinNext)
                    if tmpWl > 0.0:
                        fittingWl[nFitPts] = tmpWl
                        fittingI[nFitPts] = tmpI
                        fittingSig[nFitPts] = tmpSig
                        fittingOrder[nFitPts] = orderReplace
                        nFitPts += 1
                    
            binWlStart = obsWl[i]
            binWlEnd = binWlStart + par.velBin/c*binWlStart
            iBinStart = i
    fittingWl = fittingWl[:nFitPts]
    fittingI = fittingI[:nFitPts]
    fittingSig = fittingSig[:nFitPts]
    fittingOrder = fittingOrder[:nFitPts]
    return fittingWl, fittingI, fittingSig, fittingOrder
    

#Fit a polynomial to each spectral order, using the best points chosen.
#Allow choice of polynomial type p=regular polynomial, c=Chebyshev l=Legendre polynomial
#Returns a list of fits to the continuum points for each spectral order.
def fitPoly(obsWl, ords, fittingOrder, fittingWl, fittingI, fittingSig, polyDegs, polyType='p'):
    #Run fit on each spectral order
    if polyType == 'p':
        from numpy.polynomial import polynomial as Poly
    elif polyType == 'c':
        from numpy.polynomial import chebyshev as Cheb
    elif polyType == 'l':
        from numpy.polynomial import legendre as Leg
    else:
        print('ERROR: trying to use an unknown polynomial type!')
        return 
    fitIvals=[]
    polyfitvals=[]
    for i in range(ords.numOrders):
        #Get points in this spectral order, a brute force solution
        iOrd = np.where(fittingOrder == i)
        if fittingWl[iOrd].shape[0] < polyDegs[i]+1:
            polyDegO = fittingWl[iOrd].shape[0]-1
            print('WARNING: in spectral order {:} too few points to constrain polynomial fit ({:}), assume degree {:}'.format(i+1, fittingWl[iOrd].shape[0], polyDegO))
            if polyDegO < 0:
                print('ERROR: not enough points to fit spectral order {:} ({:}) !'.format(i+1, fittingWl[iOrd].shape[0]))
                fitIvals += [np.ones(obsWl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1].shape)]
                continue
        else:
            polyDegO = polyDegs[i]

        #Centering the polynomial x values on zero makes them more well conditioned for the fit,
        #so shift everything by the mean of the fitting wavelengths
        meanFitWl = np.mean(fittingWl[iOrd])
        fittingWlShift = fittingWl[iOrd] - meanFitWl
        obsWlShift = obsWl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1] - meanFitWl
        
        #Run the actual fit with numpy's routines
        if polyType == 'p':
            polyfitval = Poly.polyfit(fittingWlShift, fittingI[iOrd], polyDegO, w=1./fittingSig[iOrd])
            polyfitvals += [polyfitval]
            fitIvals += [Poly.polyval(obsWlShift, polyfitval)]
        elif polyType == 'c':
            chebfitVal = Cheb.chebfit(fittingWlShift, fittingI[iOrd], polyDegO, w=1./fittingSig[iOrd])
            polyfitvals += [chebfitVal]
            fitIvals += [Cheb.chebval(obsWlShift, chebfitVal)]
        elif polyType == 'l':  
            legfitVal = Leg.legfit(fittingWlShift, fittingI[iOrd], polyDegO, w=1./fittingSig[iOrd])
            polyfitvals += [legfitVal]
            fitIvals += [Leg.legval(obsWlShift, legfitVal)]
            
    return fitIvals


#Merge spectral orders
#Can take a variable number of input spectra for Stokes parameters 
#or errorbars in a multi-column spectrum
def mergeOrders(ords, wl, *specList):
    #Simplistic merging of spectral orders, by just truncating them
    #at the midpoint of order overlap.
    merSpecList = []
    for spec in specList:
        merSpec = np.zeros(0)
        merWl = np.zeros(0)
        wlStartMid = ords.wlOrderStart[0]
        for i in range(ords.numOrders):
            #Get the midpoint of any overlap
            if i >= ords.numOrders-1:
                wlEndMid = ords.wlOrderEnd[-1]
            else:
                wlEndMid = (ords.wlOrderEnd[i]+ords.wlOrderStart[i+1])*0.5
                
            indRange = (wl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1] >= wlStartMid) \
                & (wl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1] <= wlEndMid)
            merSpec = np.append(merSpec, spec[ords.iOrderStart[i]:ords.iOrderEnd[i]+1][indRange])
            if len(merSpecList) == 0:  #only merge the wavelength array once
                merWl = np.append(merWl, wl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1][indRange])
            wlStartMid = wlEndMid
        if len(merSpecList) == 0:
            merSpecList += [merWl]
        merSpecList += [merSpec]
    
    return tuple(merSpecList) #return a tuple of variable length


#Scale input wavelngths array by a value (if it is positive and not 1)
def scaleWavelength(wl, scale):
    wlScaled = wl
    if scale > 0. and scale != 1.0:
        wlScaled = wl*scale
    return wlScaled

#Convert wavelengths between air and vacuum assuming they are Angstroms
#flag = -1: air-to-vacuum, flag = +1: vacuum-to-air, flag = 0 (or other) do nothing
def convertAirVacuum(wl, flag):
    # using the formula from VALD3's website:
    # http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    # From their documentation:
    # 
    # For the vacuum to air conversion the formula from Donald Morton (2000, ApJ. Suppl., 130, 403)
    # is used for the refraction index, which is also the IAU standard:
    # n = 1 + 0.0000834254 + 0.02406147 / (130 - s^2) + 0.00015998 / (38.9 - s^2), 
    # where s = 10^4 / lambda_vac and lambda_vac is in Angstroms.
    # The conversion is then: lambda_air = lambda_vac / n. 
    # This formula comes from Birch and Downs (1994, Metrologia, 31, 315) and applies to 
    # dry air at 1 atm pressure and 15 C with 0.045% CO2 by volume. The corrections to 
    # Edlen (1953, J. Opt. Soc. Am., 43, 339) are less than 0.0001 A at 2000 A and less 
    # than 0.001 A at 30000 A.
    # 
    # The opposite conversion (air-to-vacuum) is less trivial because n depends on lambda_vac 
    # and conversion equations with sufficient precision are not readily available. 
    # VALD3 tools use the following solution derived by N. Piskunov:
    # n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s^2) + 0.0001599740894897 / (38.92568793293 - s^2), 
    # where s = 10^4 / lambda_air and the conversion is: lambda_vac = lambda_air * n.
    
    wlConv = wl
    if flag == -1:  #air to vacuum
        swl = 1e4/wl
        refraction = 1.0 + 0.00008336624212083 + 0.02408926869968/(130.1065924522 - swl**2) + 0.0001599740894897/(38.92568793293 - swl**2);

        wlConv = wl*refraction
    elif flag == 1:  #vacuum to air
        swl = 1e4/wl
        refraction = 1.0 + 0.0000834254 + 0.02406147/(130.0 - swl**2) + 0.00015998/(38.9 - swl**2)
        wlConv = wl/refraction
    return wlConv
      

#Write a spectrum to an input file
#Can use multiple input colmns of spectra, for Stokes parameters and and uncertainty
#Uses the input nObsCol value to specify the format that is written
#Protects against some erronious values in Stokes I.
def writeSpec(fname, nObsCol, *cols):
    print('Saving to {:}'.format(fname))
    if nObsCol != len(cols):
        print('ERROR: got an incorrect number of columns for writting')
        print('got {:} cols for {:} arrays.  Nothing saved.'.format(nObsCol, len(cols)))
        return
    wl = cols[0]
    obsI =  np.where(cols[1] > 0., cols[1], 0.0)
    
    fOut = open(fname, 'w')    
    if nObsCol == 2:
        for i in range(wl.shape[0]):
            fOut.write('{:10.4f} {:11.4e}\n'.format(wl[i], obsI[i]))
            
    if nObsCol == 3:
        for i in range(wl.shape[0]):
            fOut.write('{:10.4f} {:11.4e} {:11.4e}\n'.format(wl[i], obsI[i], cols[2][i]))
            
    if nObsCol == 6:
        for i in range(wl.shape[0]):
            fOut.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(
                wl[i], obsI[i], cols[2][i], cols[3][i], cols[4][i], cols[5][i]))

    #if nObsCol == 10:
    #    for i in range(wl.shape[0]):
    #        fOut.write('{:10.4f} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n'.format(
    #            obsI, cols[1][i], cols[2][i], cols[3][i], cols[4][i], cols[5][i]))

    fOut.close()
