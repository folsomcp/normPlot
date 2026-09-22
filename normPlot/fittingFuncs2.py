#Helper functions for continuum normalization fitting

import numpy as np
import scipy.constants
c = scipy.constants.c*1e-3 #km/s

supportedPolyTypes = ('Chebyshev', 'Legendre', 'Geometric', 'Spline', 'SmSpline', 'SplRep')

#Save polynomial degrees and type into one object
class polySet():
    def __init__(self, ords, obsWl, obsI, velBin=500.,
                 polyType=supportedPolyTypes[0]):
        self.type = polyType
        self.degs = [5]*ords.numOrders
        self.nknots = [4]*ords.numOrders
        self.splam = [1.0]*ords.numOrders

        # For the smoothed spline fit:
        # Use a heuristic to roughly estimate the regularization parameter lambda
        # Typically the chi^2 and regularization terms should be of the same
        # order of magnitude in regularized fitting.  Here we use chi^2 with no
        # error bars (no weights), since we are unlikely to get statistically 'good'
        # fits, and that makes lambda a bit easier to estimate.
        # The regularization term is the second derivative of the spline squared,
        # integrated over the wavelength range used.
        # As a rough order-of-magnitude heuristic for this:
        # chi^2 ~~ npts * (how_close * specI)^2
        # regularization ~~ (how_flat * specI / wl^2)^2 * widthWl
        #                ~~ how_flat^2 * specI^2 * widthWl / wl^4
        # With tuned parameters now_close and how_flat
        # The lambda = chi^2/regularization
        # lambda ~~ (how_close/how_flat)^2 * (npts * wl^4 / widthWl)
        # This needs to be evaluated for each order
        for i in range(ords.numOrders):
            #Get points in this spectral order
            iOrd = (ords.obsOrder == i)
            obsWlOrd = obsWl[iOrd]
            deltaWl = obsWlOrd[-1] - obsWlOrd[0]
            midWl = 0.5*(obsWlOrd[-1] + obsWlOrd[0])
            # Get the number of fitting points likely to be used
            # (as the window size in velocity / velocity search bin size)
            velSize = c*deltaWl/midWl
            numFitPts = np.rint(velSize/velBin)

            lamb = ((1e-5)**2)*numFitPts*(midWl**4)/deltaWl
            # Now round lambda to just 1 significant figure,
            # to make it cleaner to display and easier to change
            # First get the number of digits in the decimal number lambda
            nDigits = np.floor(np.log10(lamb)).astype('int')
            lambRound = np.around(lamb, -nDigits) #then round to that many digits
            self.splam[i] = lambRound


    def readPolyParams(self, polynomialsName):
        """
        Read in the type of fitting function, and polynomial degree for each
        spectral order, from a file.  Verify that the number of orders read
        matches the number of orders needed, and correct if necessary.
        """
        try:
            fPolyDeg = open(polynomialsName, 'r')
        except IOError:
            print('did not find {:}, using defaults'.format(polynomialsName))
            self.type = supportedPolyTypes[0]
        else:
            print('reading polynomial degrees from {:}'.format(polynomialsName))
            polyDegs = []
            i=0
            for line in fPolyDeg:
                if len(line.split()) > 0:
                    if line.split()[0][0] != '#':
                        if i == 0: # check if this is a header line (support older versions with no header)
                            try: # if this is an int (no header)
                                float(line.split()[0])
                            except ValueError: # if this is a string (a polynomial type)
                                tstPolyType = line.strip()
                                if tstPolyType in supportedPolyTypes:
                                    self.type = tstPolyType
                                else:
                                    raise ValueError('Unsupported polynomial type read from '
                                                     'input file! {:}'.format(tstPolyType))
                                continue
                        if self.type == 'SmSpline':
                            polyDegs += [float(line.split()[0])]                                
                        else:
                            polyDegs += [int(line.split()[0])]                                
                        i += 1
            fPolyDeg.close()
            
            # If we have read the file successfully, check the number of orders
            # read against the number needed, and correct if necessary
            if len(polyDegs) > len(self.degs):
                print('Read too many polynomial degrees ({:}) for the number of '
                      'spectral orders ({:}), ignoring the extras'.format(
                          len(polyDegs), len(self.degs)))
                polyDegs = polyDegs[:len(self.degs)]
            elif len(polyDegs) < len(self.degs):
                if self.type == 'Spline':
                    print('Read too few numbers of knots ({:}) for the number of '
                          'spectral orders ({:}), padding with default {:}'.format(
                              len(polyDegs), len(self.degs), self.nknots[-1]))
                elif self.type == 'SmSpline':
                    print('Read too few spline smoothing values ({:}) for the number of '
                          'spectral orders ({:}), padding with default {:}'.format(
                              len(polyDegs), len(self.degs), self.splam[-1]))
                else:
                    print('Read too few polynomial degrees ({:}) for the number of '
                          'spectral orders ({:}), padding with default {:}'.format(
                              len(polyDegs), len(self.degs), self.degs[-1]))
            # Store the read values into the appropriate place based on type
            if self.type == 'Spline':
                self.nknots[0:len(polyDegs)] = polyDegs
            elif self.type == 'SmSpline':
                self.splam[0:len(polyDegs)] = polyDegs
            else:
                self.degs[0:len(polyDegs)] = polyDegs
        return


#Save a few basic control parameters for the routine,
#Use a little object for this, to make passing to other functions easier.
class controlPars():
    def __init__(self, averageLen=1, velBin=1, bMergeOrd=False, bFillEdgeGaps=False,
                 outWaveScale=1.0, outConvertAirVac=0):
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
    def linkPolyDegExclude(self, polys, bFittable, obsWl, ords):
        self.polys = polys
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
        fpoly.write('{:}\n'.format(self.polys.type))
        i = 1
        for polyDeg in self.polys.degs:
            fpoly.write('{:n}  {:n}\n'.format(polyDeg, i))
            i += 1
        fpoly.close()

        #Infer exclude regions (wavelength ranges) from the array of fittable points bFittable
        print('saving exclude regions to {:}'.format(fnameExclude))
        #First strip out order overlap, since that would just add extraneous exclude regions
        merWl, merFit = mergeOrders(self.ords, self.obsWl, self.bFittable)
        merFit = merFit.astype(int)
        
        #Find where bFittable flags change, for exclude range edges
        edges =  np.nonzero(merFit[1:-1] != merFit[0:-2])[0] + 1
        # If there are no exclude regions found
        if len(edges) < 1:
            #print('no exclude regions to save')
            fexclude = open(fnameExclude, 'w')
            fexclude.close()
            return
        
        #Use edges where the bFittable flag has gone false for a exclude region start,
        #and back to true for ends.
        starts = edges[merFit[edges] == 0] 
        ends = edges[merFit[edges] == 1]
        #Protect against exclude regions at the start or end of the spectrum
        if ends[0] < starts[0]:
            starts = np.insert(starts, 0, 0)
        if starts[-1] > ends[-1]:
            ends = np.append(ends, merFit.shape[0]-1)
        if starts.shape[0] != ends.shape[0]:
            print('ERROR calculating exclude regions, not writting to file')
            return

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
        raise ValueError('Found an unexpected number of columns ({:}) in observation: {:}'.format(nObsCol,observationName))

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
#The direction of the search is controlled by iDer (+1 = forward -1 = backward).
def lookToNextOrderForGaps(obsWl, obsIavg, obsSig, bFittable, par,
                           iBinLast, iDer, npts=1):
    if(iDer != 1 and iDer != -1): raise ValueError
    fittingWl = []
    fittingI = []
    fittingSig = []
    nptsFound = 0
    #Run forward or backward (using iDer) looking for the next good bin with a good point
    indEnd = obsWl.shape[0]
    if iDer < 0: indEnd = 0
    binWlStart = obsWl[iBinLast]
    binWlEnd = binWlStart + iDer*par.velBin/c*binWlStart
    iBinStart = iBinLast
    for j in range(iBinLast, indEnd, iDer):
        #If we are in the wrong part of an overlap,
        #just update the start index for the search bin, so we search
        #a continuous (in wl) set of pixels from binWlStart to binWlEnd
        if (iDer > 0 and obsWl[j] < binWlStart) or (iDer < 0 and obsWl[j] > binWlStart):
            iBinStart = j
        #If we hit a bin end for any spectral order
        if (iDer > 0 and obsWl[j] > binWlEnd) or (iDer < 0 and obsWl[j] < binWlEnd):
            #(protect against zero size bins)
            if iBinStart != j:
                indMax = np.argmax(obsIavg[iBinStart:j:iDer])*iDer + iBinStart
            else:
                indMax = j - iDer
            #If this point is not in an exclude region,
            # save it and exit the loop
            if bFittable[indMax] == True:
                inpos = len(fittingWl) # insert points at the end of the list,
                if iDer < 0: inpos = 0 # or the start, depending on search direction
                fittingWl.insert(inpos, obsWl[indMax])
                fittingI.insert(inpos, obsIavg[indMax])
                fittingSig += [obsSig[indMax]]
                nptsFound += 1
                if nptsFound >= npts: # if we have enough points, quit
                    return fittingWl, fittingI, fittingSig
            # Lastly, update for the next velocity bin
            binWlStart = obsWl[j]
            binWlEnd = binWlStart + iDer*par.velBin/c*binWlStart
            iBinStart = j

    return fittingWl, fittingI, fittingSig



#Get the highest ('most likely continuum') point in each bin
#The exclude points that are inside an exclude region
def getBestInBin(obsWl, obsIavg, obsSig, obsOrder, bFittable, par, polyType):
    fittingWl = np.zeros(obsWl.size + 6) #initialize arrays with a bit of extra space
    fittingI = np.zeros(obsWl.size + 6)
    fittingSig = np.zeros(obsWl.size + 6)
    fittingOrder = np.zeros(obsWl.size + 6, dtype=int)
    if polyType == 'Spline':
        nFillPts = 2
    elif polyType == 'SmSpline':
        nFillPts = 3
    else:
        nFillPts = 1
    # Set up the first velocity bin
    iBinStart = 0
    binWlStart = obsWl[0]
    binWlEnd = binWlStart + par.velBin/c*binWlStart
    nFitPts = 0
    # loop through pixels, defining velocity bins as we go,
    # and resetting the bin when changing orders
    for i in range(obsWl.shape[0]):
        # If we have completed a velocity bin, or hit a spectral order end
        # (otherwise do nothing, just iterate to the next pixel)
        if (obsWl[i] > binWlEnd) | (obsOrder[i] != obsOrder[iBinStart]):
            #protect against too small bins (with no points)
            if iBinStart < i:
                indMax = np.argmax(obsIavg[iBinStart:i]) + iBinStart
            else:
                indMax = i-1
            
            #If the best point in this bin is not in an exclude region
            if bFittable[indMax] == True:
                fittingWl[nFitPts] = obsWl[indMax]
                fittingI[nFitPts] = obsIavg[indMax]
                fittingSig[nFitPts] = obsSig[indMax]
                fittingOrder[nFitPts] = obsOrder[indMax]
                nFitPts += 1
            else:
                #if the best point is excluded, add a filler point at
                #the start/end of a spectral order
                if par.lookToNextOrderForGaps:
                    tmpWl = []
                    # Check if this is an order end or an order start
                    # then run the function to get replacement point(s)
                    if obsOrder[i] > obsOrder[i-1]: # at an order end
                        tmpWl, tmpI, tmpSig = lookToNextOrderForGaps(
                            obsWl, obsIavg, obsSig, bFittable, par,
                            i-1, 1, nFillPts)
                    elif obsOrder[iBinStart] > obsOrder[iBinStart-1]: # at an order start
                        tmpWl, tmpI, tmpSig = lookToNextOrderForGaps(
                            obsWl, obsIavg, obsSig, bFittable, par,
                            iBinStart, -1, nFillPts)
                    # If there are replacement points to use, add them
                    if len(tmpWl) > 0:
                        nFillPtsUsed = len(tmpWl)
                        fittingWl[nFitPts:nFitPts + nFillPtsUsed] = tmpWl
                        fittingI[nFitPts:nFitPts + nFillPtsUsed] = tmpI
                        fittingSig[nFitPts:nFitPts + nFillPtsUsed] = tmpSig
                        fittingOrder[nFitPts:nFitPts + nFillPtsUsed] = obsOrder[indMax]
                        nFitPts += nFillPtsUsed
            # And lastly, set up the next velocity bin
            binWlStart = obsWl[i]
            binWlEnd = binWlStart + par.velBin/c*binWlStart
            iBinStart = i
            # If the bin after next would hit an order end, and be too small,
            # then extend this bin's size to avoid very small bins
            indOrderEnds = np.nonzero(obsOrder > obsOrder[i])[0] - 1
            if len(indOrderEnds) < 1: indOrderEnds = [-1]
            wlOrderEnd = obsWl[indOrderEnds[0]]
            if wlOrderEnd - binWlEnd < 0.33*(par.velBin/c*binWlEnd):
                binWlEnd = wlOrderEnd

    fittingWl = fittingWl[:nFitPts]
    fittingI = fittingI[:nFitPts]
    fittingSig = fittingSig[:nFitPts]
    fittingOrder = fittingOrder[:nFitPts]
    return fittingWl, fittingI, fittingSig, fittingOrder


#Fit a polynomial to each spectral order, using the best points chosen.
#Allow choice of polynomial type p=regular polynomial, c=Chebyshev l=Legendre polynomial
#Returns a list of fits to the continuum points for each spectral order.
def fitPoly(obsWl, ords, fittingOrder, fittingWl, fittingI, fittingSig, polys):
    #Run fit on each spectral order
    if polys.type == 'Geometric':
        from numpy.polynomial import polynomial as Poly
    elif polys.type == 'Chebyshev':
        from numpy.polynomial import chebyshev as Cheb
    elif polys.type == 'Legendre':
        from numpy.polynomial import legendre as Leg
    elif polys.type == 'Spline':
        from scipy.interpolate import make_lsq_spline
    elif polys.type == 'SmSpline':
        from scipy.interpolate import make_smoothing_spline, splrep, splev
    elif polys.type == 'SplRep':
        from scipy.interpolate import make_smoothing_spline, splrep, splev
    else:
        raise ValueError('Trying to fit with an unknown polynomial type! {:}'.format(
            polys.type))
    fitIvals=[]
    polyfitvals=[]
    for i in range(ords.numOrders):
        #Get points in this spectral order, a brute force solution
        iOrd = np.where(fittingOrder == i)
        numObsPts = fittingWl[iOrd].shape[0]
        #Check that this order has points to fit
        if polys.type == 'SmSpline' and numObsPts < 5:
            fitIvals += [np.ones(ords.iOrderEnd[i] - ords.iOrderStart[i] + 1)]
            print('ERROR: not enough points to fit spectral order '
                  '{:} ({:} points) need 5!'.format(i+1, numObsPts))
            continue
        elif polys.type == 'Spline' and numObsPts < 4:
            fitIvals += [np.ones(ords.iOrderEnd[i] - ords.iOrderStart[i] + 1)]
            print('ERROR: not enough points to fit spectral order '
                  '{:} ({:} points) need 4!'.format(i+1, numObsPts))
            continue
        elif numObsPts < 2:
            fitIvals += [np.ones(ords.iOrderEnd[i] - ords.iOrderStart[i] + 1)]
            print('ERROR: not enough points to fit spectral order '
                  '{:} ({:} points) need 2!'.format(i+1, numObsPts))
            continue

        if polys.type == 'Geometric' or polys.type == 'Chebyshev' or polys.type == 'Legendre':
            polyDegO = polys.degs[i]
            if numObsPts < polys.degs[i] + 1:
                polyDegO = numObsPts - 1
                print('WARNING: in spectral order {:} too few points to constrain '
                      'polynomial fit ({:}), assume degree {:}'.format(
                          i+1, numObsPts, polyDegO))
        elif polys.type == 'Spline':
            polyDegO = polys.nknots[i]
            # There are typically (4 + number of interior knots) free parameters in
            # the spline fit for make_lsq_spline, so make sure there are at least that 
            # many datapoints for constraints.  (2 + knots, counting exterior end knots)
            if numObsPts < polyDegO + 4:
                newPolyDegO = max(numObsPts - 4, 0)
                print('WARNING: in spectral order {:} too few points ({:}) to constrain '
                      'spline with {:} interior knots, reducing to {:}'.format(
                          i+1, numObsPts, polyDegO, newPolyDegO))
                polyDegO = newPolyDegO
        elif polys.type == 'SmSpline':
            polyDegO = polys.splam[i]
        else:
            polyDegO = polys.degs[i]

        #Centering the polynomial x values on zero makes them more well conditioned for the fit,
        #so shift everything by the mean of the fitting wavelengths
        meanFitWl = np.mean(fittingWl[iOrd])
        fittingWlShift = fittingWl[iOrd] - meanFitWl
        obsWlShift = obsWl[ords.iOrderStart[i]:ords.iOrderEnd[i]+1] - meanFitWl
        
        #Run the actual fit with numpy's routines
        if polys.type == 'Geometric':
            polyfitval = Poly.polyfit(fittingWlShift, fittingI[iOrd], polyDegO,
                                      w=1./fittingSig[iOrd])
            polyfitvals += [polyfitval]
            fitIvals += [Poly.polyval(obsWlShift, polyfitval)]
        elif polys.type == 'Chebyshev':
            chebfitVal = Cheb.chebfit(fittingWlShift, fittingI[iOrd], polyDegO,
                                      w=1./fittingSig[iOrd])
            polyfitvals += [chebfitVal]
            fitIvals += [Cheb.chebval(obsWlShift, chebfitVal)]
        elif polys.type == 'Legendre':
            legfitVal = Leg.legfit(fittingWlShift, fittingI[iOrd], polyDegO,
                                   w=1./fittingSig[iOrd])
            polyfitvals += [legfitVal]
            fitIvals += [Leg.legval(obsWlShift, legfitVal)]
        elif polys.type == 'Spline':
            splDegree = 3 # spline degree (cubic spline = 3)
            # make knot positions for the spline.
            ## To just evenly space 'polyDegO' knots across the order.
            ## There should be degree+1 boundary knots padding the ends of the x range
            #knots = np.linspace(fittingWlShift[0], fittingWlShift[-1], polyDegO + 2)
            #knots = np.concatenate((np.ones(splDegree)*fittingWlShift[0],
            #                       knots,
            #                       np.ones(splDegree)*fittingWlShift[-1]))
            # Place knots at datapoints, spaced by datapoint number
            # Get polyDegO fractions of the way through the array of positions,
            # then convert those fractions to array indices.
            fracKnots = np.linspace(0.0, 1.0, polyDegO + 2)
            iknots = np.rint(fracKnots*(numObsPts - 1)).astype(int)
            knots = fittingWlShift[iknots]
            # There should be degree+1 boundary knots padding the ends of the x range
            knots = np.concatenate((np.ones(splDegree)*fittingWlShift[0],
                                   knots,
                                   np.ones(splDegree)*fittingWlShift[-1]))
            spl = make_lsq_spline(fittingWlShift, fittingI[iOrd],
                                  knots, k=splDegree, w=1./fittingSig[iOrd])
            polyfitvals += [spl]
            fitIvals += [spl(obsWlShift)]
        elif polys.type == 'SmSpline':
            splDegree = 3 # spline degree (cubic spline = 3)
            lam = polyDegO
            spl = make_smoothing_spline(fittingWlShift, fittingI[iOrd], lam=lam)
            # To include error bars (but that makes setting lambda harder) use:
            #spl = make_smoothing_spline(fittingWlShift, fittingI[iOrd],
            #                            lam=lam, w=1./fittingSig[iOrd]**2)
            polyfitvals += [spl]
            fitIvals += [spl(obsWlShift)]
            
        elif polys.type == 'SplRep':
            splDegree = 3 # spline degree (cubic spline = 3)
            #knots = np.linspace(fittingWlShift[0], fittingWlShift[-1], 5, endpoint=False)
            tck = splrep(fittingWlShift, fittingI[iOrd],
                         w=1./fittingSig[iOrd], k=splDegree,
                         s=(fittingWlShift.size)*polyDegO)  #, t=knots[1:])
            # Adding manual knots seems to override smoothing
            # (Or maybe there just aren't enough degrees of freedom for smoothing algorithm in that case?
            # But smoothed version seems to have a variable number of knots,
            # so knot choice seems to be part of that algorithm).
            #print(tck[0])
            print(fittingWlShift.size, tck[0].size, tck[1].size)
            #print(fittingWlShift.size - np.sqrt(2*fittingWlShift.size),
            #      fittingWlShift.size + np.sqrt(2*fittingWlShift.size),
            #      fittingWlShift.size - splDegree,
            #      fittingWlShift.size - tck[1].size)
            chi2 = np.sum(((fittingI[iOrd] - splev(fittingWlShift, tck))/fittingSig[iOrd])**2)
            print('chi2:', chi2,
                  'target:', fittingWlShift.size*polyDegO,
                  'reduced chi2', chi2/(fittingWlShift.size - (tck[1].size - polyDegO - 1)) )
            #print(tck)
            polyfitvals += [tck]
            fitIvals += [splev(obsWlShift, tck)]
    return fitIvals


#Merge spectral orders
#Can take a variable number of input spectra for Stokes parameters 
#or errorbars in a multi-column spectrum
def mergeOrders(ords, wl, *specList):
    #Simplistic merging of spectral orders, by just truncating them
    #at the midpoint of order overlap.
    overlapFracMerge = 0.5
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
                #if ords.wlOrderEnd[i] >= 840.: #wavelength dependence could look like
                #    overlapFracMerge = 0.2
                #else: 
                #    overlapFracMerge = 0.5
                wlEndMid = (ords.wlOrderEnd[i]*overlapFracMerge
                            +ords.wlOrderStart[i+1]*(1-overlapFracMerge))
                
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
