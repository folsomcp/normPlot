#Classes/Funcs for interactive GUI

import numpy as np
from matplotlib import colors
import tkinter.messagebox as messagebox
import tkinter as tk
import tkinter.ttk as ttk
import tkinter.font as tkfont
import fittingFuncs2 as ff
import time


#Functions for panning and Zooming
#A little event manager thing for dealing with callbacks
class bViewFuncs:
    def __init__(self, ax, plObs, canvas):
        self.ax = ax
        self.plObs = plObs
        self.canvas = canvas
        self.zoomRecActive = False
        return
    def autoScale(self, *event):
        fracSpace = 0.05
        xmin = np.min(self.plObs[0].get_xdata())
        xmax = np.max(self.plObs[0].get_xdata())
        ymin = np.min(self.plObs[0].get_ydata())
        ymax = np.max(self.plObs[0].get_ydata())
        xspan = xmax-xmin
        yspan = ymax-ymin
        self.ax.set_xlim(xmin-xspan*fracSpace, xmax+xspan*fracSpace)
        self.ax.set_ylim(ymin-yspan*fracSpace, ymax+yspan*fracSpace)
        self.canvas.draw()
        return
    def autoScaleY(self, *event):
        fracSpace = 0.05
        xmin, xmax = self.ax.get_xlim()
        iInViewX = (self.plObs[0].get_xdata() > xmin) & (self.plObs[0].get_xdata() < xmax)
        ymin = np.min(self.plObs[0].get_ydata()[iInViewX])
        ymax = np.max(self.plObs[0].get_ydata()[iInViewX])
        yspan = ymax-ymin
        self.ax.set_ylim(ymin-yspan*fracSpace, ymax+yspan*fracSpace)        
        self.canvas.draw()
        return
    def zoomIn(self, *event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        modFracZoom = fracZoom/(1.+2.*fracZoom) #modify so zoom in reverses zoom out exactly
        self.ax.set_xlim(xmin+xspan*modFracZoom, xmax-xspan*modFracZoom)
        self.ax.set_ylim(ymin+yspan*modFracZoom, ymax-yspan*modFracZoom)
        self.canvas.draw()
        return
    def zoomOut(self, *event):
        fracZoom = 0.1
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        xspan = xmax-xmin
        yspan = ymax-ymin
        self.ax.set_xlim(xmin-xspan*fracZoom, xmax+xspan*fracZoom)
        self.ax.set_ylim(ymin-yspan*fracZoom, ymax+yspan*fracZoom)
        self.canvas.draw()
        return
    def panLeft(self, *event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax-xmin
        self.ax.set_xlim(xmin-xspan*fracPan, xmax-xspan*fracPan)
        self.canvas.draw()
        return
    def panRight(self, *event):
        fracPan = 0.2
        xmin, xmax = self.ax.get_xlim()
        xspan = xmax-xmin
        self.ax.set_xlim(xmin+xspan*fracPan, xmax+xspan*fracPan)
        self.canvas.draw()
        return
    def panUp(self, *event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax-ymin
        self.ax.set_ylim(ymin+yspan*fracPan, ymax+yspan*fracPan)
        self.canvas.draw()
        return
    def panDown(self, *event):
        fracPan = 0.2
        ymin, ymax = self.ax.get_ylim()
        yspan = ymax-ymin
        self.ax.set_ylim(ymin-yspan*fracPan, ymax-yspan*fracPan)
        self.canvas.draw()
        return
        #Optionally check if the event is not comming from a widget we don't want to respond to
        #if event:
        #    print(event[0].widget, event[0].type)
        #    #parse .widget to make sure it's not an Entry ?
        
    def zoomRec(self, *event):
        #Zoom to a rectangle selected by the user with the mouse
        if not self.zoomRecActive:
            self.zoomRecActive = True
            self.canvasWidget = self.canvas.get_tk_widget()
            self.numClick = 0
            self.x0rec, self.y0rec, self.x1rec, self.y1rec = 0., 0., 0., 0.
            self.oldBindClick = self.canvasWidget.bind('<Button-1>')
            self.canvasWidget.bind('<Button-1>', self.onClickRec)
            self.oldBindClick3 = self.canvasWidget.bind('<Button-3>')
            self.canvasWidget.bind('<Button-3>', self.onClickCancel)
            self.oldCursor = self.canvasWidget.cget('cursor')
            self.canvasWidget.config(cursor='crosshair')
        return
    def zoomRecDeactivate(self):
        self.canvasWidget.bind('<Button-1>', self.oldBindClick)
        self.canvasWidget.bind('<Button-3>', self.oldBindClick3)
        self.canvasWidget.config(cursor=self.oldCursor)
        self.zoomRecActive = False
        #when deactivating also set the plot view to the selected range
        if self.x1rec > 0:
            #get the canvas position in window pixel coordinates
            x0full, x1full = self.canvas.figure.gca().bbox.intervalx
            y0inv, y1inv = self.canvas.figure.gca().bbox.intervaly
            #and correct for matplotlib using y in the opposide direction from tk
            height = self.canvas.figure.bbox.height
            y0full = height - y1inv
            y1full = height - y0inv
            axesxrange = x1full - x0full
            axesyrange = y1full - y0full
            #get the current view range in data coordinates
            axdatax0, axdatax1, axdatay0, axdatay1 = self.canvas.figure.gca().axis()
            dataxrange = axdatax1 - axdatax0
            datayrange = axdatay1 - axdatay0
            #convert the selected window pixel range into data coordinates
            datax0 = (self.x0rec - x0full)/(axesxrange)*(dataxrange) + axdatax0
            datax1 = (self.x1rec - x0full)/(axesxrange)*(dataxrange) + axdatax0
            datay0 = (self.y0rec - y0full)/(axesyrange)*(-datayrange) + axdatay1
            datay1 = (self.y1rec - y0full)/(axesyrange)*(-datayrange) + axdatay1
            #set the plotted range into the selected range in data coordiantes
            if abs(datax1 - datax0) > 0. and abs(datay1 - datay0) > 0.:
                self.ax.set_xlim(min(datax0, datax1), max(datax0, datax1))
                self.ax.set_ylim(min(datay0, datay1), max(datay0, datay1))
                self.canvas.draw()
        return
    def onClickRec(self, event):
        #On the first click start drawing a selection rectangle
        if self.numClick == 0:
            self.numClick += 1
            self.oldBindMove = self.canvasWidget.bind('<Motion>')
            self.canvasWidget.bind('<Motion>', self.onMoveRec)
            self.x0rec = event.x
            self.y0rec = event.y
            self.lastrect = self.canvasWidget.create_rectangle(
                self.x0rec, self.y0rec, self.x0rec+1, self.y0rec+1, width=1)
        #On the second click deactivate and finish the zoom
        else:
            self.x1rec = event.x
            self.y1rec = event.y
            self.canvasWidget.delete(self.lastrect)
            self.canvasWidget.bind('<Motion>', self.oldBindMove)
            self.zoomRecDeactivate()
        return
    def onClickCancel(self, event):
        #cancel the zoom to rectangle, usually in right click
        self.canvasWidget.delete(self.lastrect)
        if self.numClick > 0:
            self.canvasWidget.bind('<Motion>', self.oldBindMove)
        self.zoomRecDeactivate()
        return
    def onMoveRec(self, event):
        #delete and redraw the rectangle when the mouse moves
        x1 = event.x
        y1 = event.y
        self.canvasWidget.delete(self.lastrect)
        self.lastrect = self.canvasWidget.create_rectangle(
            self.x0rec, self.y0rec, x1, y1, width=1)
        return
              

#Recalculate the running average, using the usere's input values
class changeRunningAvg:
    def __init__(self, txt_runningAvg, obsI, obsIavg, ords, par):
        self.txt_runningAvg = txt_runningAvg
        self.obsI = obsI
        self.obsIavg = obsIavg
        self.ords = ords
        self.par = par
    def redoAverage(self, event):
        text = self.txt_runningAvg.get()
        try:
            int(text)
        except ValueError:
            print('got a non-integer value for the running average length: {:}'.format(text))
            messagebox.showwarning(
                "Bad value",
                "Need an integer for the running average length.\n"
                +"You entered:\n{:}".format(text) )
            return False
        averageLen = int(text)
        if not(averageLen > 0 and averageLen < 1000):
            print('extrem value for the running average length: {:}'.format(averageLen))
            messagebox.showwarning(
                "Extreme value",
                "Got a value outside the range for the running average length:\n"
                +"{:}".format(text) )
            return False
        self.par.averageLen = averageLen
        if event is not None:
            self.obsIavg[:] = ff.runningAvg(self.obsI, self.ords, averageLen)
            #No need to draw unless you also re-calculate/re-fit something
            #for dline in self.plObs:
            #    dline.set_ydata(obsIavg)
            #    self.canvas.draw()
        return True


#Read and rerun the selecton of best continuum points in the velocity sized spectral bins
class changeBinSize:
    def __init__(self, canvas, txt_velBin, obsWl, obsIavg, obsSig, 
                 obsOrder, par, bFittable, plFitting):
        self.canvas = canvas
        self.txt_velBin = txt_velBin
        self.obsWl = obsWl
        self.bFittable = bFittable
        self.obsIavg = obsIavg
        self.obsSig = obsSig
        self.obsOrder = obsOrder
        self.par = par
        self.plFitting = plFitting
    def redoBestInBin(self, event):
        #tStart = time.clock()
        text = self.txt_velBin.get()
        try:
            float(text)
        except ValueError:
            print('got a non-number for the velocity bin size: {:}'.format(text))
            messagebox.showwarning(
                "Bad value",
                "Need a real number for the velocity bin size.\n"
                +"You entered:\n{:}".format(text) )
            return False
        velBin = float(text)
        if not(velBin > 0. and velBin < 100000.):
            print('extrem value for the running average length: {:}'.format(velBin))
            messagebox.showwarning(
                "Extreme value",
                "Got a value outside the range for the velocity bin size:\n"
                +"{:}".format(text) )
            return False
        self.par.velBin = velBin

        if event is not None:
            fittingWl, fittingI, fittingSig, fittingOrder = ff.getBestInBin(self.obsWl, self.obsIavg, self.obsSig, self.obsOrder, self.bFittable, self.par)
            #tBinned = time.clock()
            #print('time Binned {:}'.format(tBinned-tStart))
            
            for dline in self.plFitting:
                dline.set_data(fittingWl, fittingI)
            self.canvas.draw()
        return True
        

class uiIncludeRange:
    def __init__(self, canvas, obsWl, obsI, bFittable, ords, setPlObsO):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.obsWl = obsWl
        self.obsI = obsI
        self.bFittable = bFittable
        self.ords = ords
        self.setPlObsO = setPlObsO
        self.active = False
        self.rangeSelect = rangeSelect(canvas, self)

        #Set up a style for the 'actively selecting a range' button state
        #Seems like I need to set the base style changes first,
        #then use map to set state specific differences,
        #since map doesn't provide easy access to being in
        #not any of the special states (only not a specific state)
        self.style = ttk.Style()
        self.style.configure('ActRange.TButton', relief='sunken',
                        background='#fbfbfb')
        self.style.map('ActRange.TButton',
                  relief=[('active','sunken')],
                  background=[('pressed','#eeeeee'),('active','#fbfbfb')])
        #note: order matters for map here, 
        #when the button is 'active' (mouse over) and also being pressed
        #The scope of the style names seems to be pretty wide...

    def linkButton(self, button, otherRange):
        self.button = button
        self.otherRange = otherRange
    def runSpanSelect(self):
        #First turn off the other button (only include or exclude not both!)
        bOtherActive = self.otherRange.active
        if (bOtherActive == True): 
            self.otherRange.deactivate()
        self.active = not self.active
        if self.active:
            self.button.configure(style='ActRange.TButton')
            #get a string of tk info that we can use later to restore the current bindings
            self.oldBindClick = self.canvasWidget.bind('<Button-1>')
            #override the current binding
            self.fidClick = self.canvasWidget.bind('<Button-1>',
                                                   self.rangeSelect.onClick)
            self.oldBindClick3 = self.canvasWidget.bind('<Button-3>')
            self.canvasWidget.bind('<Button-3>', self.rangeSelect.cancel)
        else:
            self.deactivate()
    def deactivate(self):
        #stop span selector...
        self.rangeSelect.deactivate()
        self.button.configure(style='TButton')
        #reset to the previous bindings once we are done
        #(Safer to bind to something different, since unbind seems to have poorly defined behaviour)
        self.canvasWidget.bind('<Button-1>', self.oldBindClick)
        self.canvasWidget.bind('<Button-3>', self.oldBindClick3)
        self.active = False
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the plotted 'not excluded' spectrum
        wlRange = [[minval, maxval],]
        ff.addIndFittable(self.bFittable, self.obsWl, wlRange)
        #update the data in the plot for each spectal order
        i = 0
        for plObsO in self.setPlObsO:
            obsWlOr    = self.obsWl[    self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            obsIOr     = self.obsI[     self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            bFittableO = self.bFittable[self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            for dline in plObsO:
                dline.set_data(obsWlOr[bFittableO], obsIOr[bFittableO])
            i += 1
        self.canvas.draw()

        
class uiExcludeRange:
    def __init__(self, canvas, obsWl, obsI, bFittable, ords, setPlObsO):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.obsWl = obsWl
        self.obsI = obsI
        self.bFittable = bFittable
        self.ords = ords
        self.setPlObsO = setPlObsO
        self.active = False
        self.rangeSelect = rangeSelect(canvas, self)
        
    def linkButton(self, button, otherRange):
        self.button = button
        self.otherRange = otherRange
    def runSpanSelect(self):
        #First turn off the other button (only include or exclude not both!)
        bOtherActive = self.otherRange.active
        if (bOtherActive == True):
            self.otherRange.deactivate()
        self.active = not self.active
        if self.active:
            self.button.configure(style='ActRange.TButton')
            #get a string of tk info that we can use later to restore the current bindings
            self.oldBindClick = self.canvasWidget.bind('<Button-1>')
            #override the current binding
            self.fidClick = self.canvasWidget.bind('<Button-1>',
                                                   self.rangeSelect.onClick)
            self.oldBindClick3 = self.canvasWidget.bind('<Button-3>')
            self.canvasWidget.bind('<Button-3>', self.rangeSelect.cancel)
        else:
            self.deactivate()
    def deactivate(self):
        #stop span selector...
        self.rangeSelect.deactivate()
        self.button.configure(style='TButton')
        #reset to the previous bindings once we are done
        #(Safer to bind to something different, since unbind seems to have poorly defined behaviour)
        self.canvasWidget.bind('<Button-1>', self.oldBindClick)
        self.canvasWidget.bind('<Button-3>', self.oldBindClick3)
        self.active = False
    def selectedWl(self, minval, maxval):
        #After a region has been selected, update the plotted 'not excluded' spectrum
        wlRange = [[minval, maxval],]
        ff.cutIndFittable(self.bFittable, self.obsWl, wlRange)
        #update the data in the plot for each spectal order
        i = 0
        for plObsO in self.setPlObsO:
            obsWlOr    = self.obsWl[    self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            obsIOr     = self.obsI[     self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            bFittableO = self.bFittable[self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1]
            for dline in plObsO:
                dline.set_data(obsWlOr[bFittableO], obsIOr[bFittableO])
            i += 1
        self.canvas.draw()

#uiIncludeRange and uiExcludeRange are almost identical,
#so you could probably do something where they are both derived
#from the same base class.  But for now it's easier to keep them separate


#Container to make sure all parameters are read from the UI
#before calling the par saveValues method to write out the parameters
class saveParams:
    def __init__(self, par, changeRunningAvg, changeBinSize):
        self.par = par
        self.chRunningAvg = changeRunningAvg
        self.chBinSize = changeBinSize
    def doSave(self, *event):
        #Check for an updated average length
        #(with the 'event' None this just saves the values)
        okChk = self.chRunningAvg.redoAverage(None)
        if okChk is False: return
        #Check for an updated bin size
        #(with the 'event' None this just saves the values)
        okChk = self.chBinSize.redoBestInBin(None)
        if okChk is False: return

        self.par.saveValues()
        return


class runFitCont:
    def __init__(self, canvas, changeRunningAvg, changeBinSize,
                 obsWl, obsI, obsSig, ords, bFittable, obsIavg,
                 par, polyDegs, setPlPoly, plFitting):
        self.canvas = canvas
        self.chRunningAvg = changeRunningAvg
        self.chBinSize = changeBinSize
        self.obsWl = obsWl
        self.obsI = obsI
        self.obsSig = obsSig
        self.ords = ords
        self.bFittable = bFittable
        self.obsIavg = obsIavg
        self.par = par
        self.polyDegs = polyDegs
        self.setPlPoly = setPlPoly
        self.plFitting = plFitting
    def refitCont(self):
        #Refit the polynomial, redo any previous steps to ensure this is up to date.
        #Check for an updated average length (with the 'event' None)
        okChk = self.chRunningAvg.redoAverage(None)
        if okChk is False: return
        #We could re-do the running average with an event other than none,
        #but to be safe it's done explicity here
        self.obsIavg[:] = ff.runningAvg(self.obsI, self.ords, self.par.averageLen)
        #Check for an updated bin size (with the 'event' None)
        okChk = self.chBinSize.redoBestInBin(None)
        if okChk is False: return
        #We may need to redo the search for best points in each bin,
        #since exclude regions or the average may have changed.
        fittingWl, fittingI, fittingSig, fittingOrder = ff.getBestInBin(
            self.obsWl, self.obsIavg, self.obsSig, self.ords.obsOrder,
            self.bFittable, self.par)
        fitIvals = ff.fitPoly(self.obsWl, self.ords, fittingOrder, fittingWl,
                              fittingI, fittingSig, self.polyDegs)

        #re-draw the fitting polynomial
        i=0
        for plPoly in self.setPlPoly:
            for dline in plPoly:
                dline.set_data(self.obsWl[self.ords.iOrderStart[i]:self.ords.iOrderEnd[i]+1],
                               fitIvals[i])
            i += 1

        #re-draw points used in the fit
        for dline in self.plFitting:
            dline.set_data(fittingWl, fittingI)
        self.canvas.draw()
        return


#Save the results of the check button for whether to merge spectral orders
#Actual merging is done just before the normalized spectrum is written.
class mergeOrders:
    def __init__(self, varMergeOrders, par):
        self.varMergeOrders = varMergeOrders
        self.par = par
    def setFlag(self):
        val = self.varMergeOrders.get()
        self.par.bMergeOrd = bool(val)
        return

    
#Container for the check button for whether to use points from adjacent spectral orders,
#when an order starts/ends in an exclude region. 
class fillEdgeGaps:
    def __init__(self, varFillEdgeGaps, par):
        self.varFillEdgeGaps = varFillEdgeGaps
        self.par = par
    def setFlag(self, *event):
        val = self.varFillEdgeGaps.get()
        self.par.lookToNextOrderForGaps = bool(val)
        return

    
#Open a new matplotlib window with a set of textboxes for setting polynomial degrees
#This is pretty inefficient, but lets you see all the polynomial orders in one place.
class newWindowDeg:
    def __init__(self, parent, polyDegs, ords, setPlObsO, setPlPoly):
        self.parent = parent
        self.polyDegs = polyDegs
        self.ords = ords
        self.setPlObsO = setPlObsO
        self.setPlPoly = setPlPoly
        self.win = None
        self.active = False
        #some params mostly for printing errors
        self.wframe = None
        self.lblErrMessage = None
        self.nrow = 0
        self.ncol = 0
        self.styleErr = ttk.Style()
        self.styleErr.configure('error.TLabel', foreground='#ff0000')
        self.errFont = tkfont.Font(font='TkDefaultFont')
        self.errFont.configure(weight='bold')
        
    def closeWindow(self, *event):
        self.active = False
        try:
            self.win.destroy()
        except:
            pass
        return

    def openWindow(self, *event):
        #If the window is already open do nothing
        if self.active:
            return
        self.active = True
        
        self.win = tk.Toplevel(self.parent)
        self.win.title("Set fitting polynomial degrees")
        #Associate this window with its parent,
        #and keep this window out of the window manager.
        self.win.transient(self.parent)
        #Set what happens when the window manager asks to close the window.
        self.win.protocol("WM_DELETE_WINDOW", self.closeWindow)

        
        #styleWframe = ttk.Style()
        #styleWframe.configure('wframe.TLabelframe', relief='flat')
        wframe = ttk.Labelframe(self.win, text='Set polynomial degree for each observation order',
                                padding="5 5 5 5")
        wframe.pack(fill=tk.BOTH, expand=1, padx=4, pady=4)
        self.wframe = wframe

        self.setChangePolyDeg = []
        setTxtDegrees = []
        nrows = 20
        for iorder in range(self.ords.numOrders):
            #Get the matplotlib colors used on the plot and conver them to Tk format strings
            colorObs  = self.setPlObsO[iorder][0].get_color()
            colorObs  = colors.to_hex(colorObs, keep_alpha=False)
            colorLine = self.setPlPoly[iorder][0].get_color()
            colorLine = colors.to_hex(colorLine, keep_alpha=False)
            styleOrd = ttk.Style()
            styleOrd.configure('order{:n}.TLabel'.format(iorder), foreground=colorObs)
            
            #Make a label and entry box for each spectral order
            lblOrderNum = ttk.Label(wframe, text='order {:n}:'.format(iorder+1),
                                    style='order{:n}.TLabel'.format(iorder),
                                    justify='right', padding=(4,0,0,0))
            txt_degree = tk.StringVar()
            txt_degree.set('{:n}'.format(self.polyDegs[iorder]))
            setTxtDegrees += [txt_degree]
            changePolyDeg = changePolyDegree(self, wframe, setTxtDegrees[iorder],
                                             self.polyDegs, iorder)
            self.setChangePolyDeg += [changePolyDeg]
            inputDegree = ttk.Entry(master=wframe, textvariable=setTxtDegrees[iorder], width=6)
            inputDegree.bind('<Key-Return>', changePolyDeg.updateDegree)
            #ToolTip(inputDegree, 'set polynomial degree')
            
            #Set up buttons in columns nrows long, then however many columns necessary
            ypos = (iorder//nrows)
            xpos = (iorder%nrows)
            lblOrderNum.grid(row=xpos, column=ypos*2, sticky=tk.E, pady=2)
            inputDegree.grid(row=xpos, column=ypos*2+1, sticky=tk.W, pady=2)

        self.nrow = nrow = xpos+1
        self.ncol = ncol = ypos*2+2
        #Let the columns expand
        for i in range(ncol):
            wframe.columnconfigure(i, weight=1)

        butApply = ttk.Button(master=wframe, text="Apply",
                              command=self.updateAllDegrees)
        butClose = ttk.Button(master=wframe, text="Close",
                              command=self.closeWindow)
        #leave space for an error message at row=nrow
        butApply.grid(row=nrow+1, column=ncol-2, sticky='S', pady=2)
        butClose.grid(row=nrow+1, column=ncol-1, sticky='S', pady=2)
        wframe.rowconfigure(nrow, weight=1)

    def updateAllDegrees(self, *event):
        #Run update degree (on the changePolyDegree holder calss),
        #for each spectral order.
        errs = 0
        self.printClear()
        for i in range(self.ords.numOrders):
            err = self.setChangePolyDeg[i].updateDegree(batch=True)
            errs += err
        return
    
    def printError(self, message, batch=False):
        #Print an error, called from changePolyDegree.updateDegree
        #Use batch=True for running this multiple times and printing multiple messages
        wraplen=self.wframe.winfo_width()
        if wraplen < 10:
            wraplen = 400
        #Generate a new message label if necessary
        if not self.lblErrMessage:
            self.lblErrMessage = ttk.Label(self.wframe, text=message,
                style='error.TLabel', font=self.errFont, wraplength=wraplen)
            self.lblErrMessage.grid(row=self.nrow, column=0, columnspan=10)
        else:
            if not batch:
                self.lblErrMessage.configure(text=message, wraplength=wraplen)
            else:  #if in batch mode, append messages
                oldtext = self.lblErrMessage['text'] + '\n'
                self.lblErrMessage.configure(text=oldtext+message, wraplength=wraplen)                
        return

    def printClear(self, batch=False):
        #Clear any error messages, called from changePolyDegree.updateDegree
        if self.lblErrMessage and not batch:
            #self.lblErrMessage.configure(text='') #alternative
            tmp = self.lblErrMessage
            self.lblErrMessage = None
            tmp.destroy()
            #destroy() is supposed to release all resources,
            #but this may still leak a small amount of memory
            #since Tkinter seems to at least save names.
            #(Hopefully the garbage collector will deallocate the
            # tkinter object when it has no reference.)
        return
     
           
class changePolyDegree:
    def __init__(self, newWindowDeg, wframe, txt_degree, polyDegs, iorder):
        self.newWindowDeg = newWindowDeg
        self.wframe = wframe
        self.txt_degree = txt_degree
        self.polyDegs = polyDegs
        self.iorder = iorder #this is just an interger, so not mutable
        
    def updateDegree(self, *event, batch=False):
        text = self.txt_degree.get()
        try:
            int(text)
        except ValueError:
            print('got a non-integer value for a polynomial degree: {:}'.format(text))
            self.newWindowDeg.printError(
                'Need an integer value for polynomial {:n} degree: {:}'.format(
                    self.iorder+1, text), batch)
            return 1
        polyDeg = int(text)
        if not(polyDeg >= 0 and polyDeg < 1000):
            print('extreme value for a polynomial degree: {:}'.format(polyDeg))
            if polyDeg < 0:
                self.newWindowDeg.printError(
                    'Need a positive value for polynomial {:n} degree: {:}'.format(
                        self.iorder+1, polyDeg), batch)
            else:
                self.newWindowDeg.printError(
                    'Need a less extreme value for polynomial {:n} degree: {:}'.format(
                        self.iorder+1, polyDeg), batch)
            return 1
        if not batch:
            self.newWindowDeg.printClear()
        self.polyDegs[self.iorder] = polyDeg
        return 0
    
    
#Open a new window for inputing parameters controling the writing of the spectrum
#e.g. scale the wavelength to change units or convert wavelength from vacuum to air
class newWindowOutPar:
    def __init__(self, parent, par):
        self.parent = parent
        self.par = par
        self.win = None
        self.active = False
        
    def cancelWindow(self, *event):
        self.active = False
        try:
            self.win.destroy()
        except:
            pass
        return

    def doneWindow(self, *event):
        self.active = False
        ok = self.setWavelengthScale()
        self.airVacConv()
        if ok:
            try:
                self.win.destroy()
            except:
                pass
        return

    def openWindow(self, *event):
        #Try making this window modal (deactivating the parent window)
        #If the window is already open do nothing
        if self.active:
            return
        self.active = True

        self.win = tk.Toplevel(self.parent)
        self.win.title("Set output parameters")
        #Associate this window with its parent,
        #and keep this window out of the window manager.
        self.win.transient(self.parent)
        #Send keyboard and mouse events only to this window,
        #mostly deactivating other windows.
        self.win.grab_set()
        #Set what happens when the window manager asks to close the window.
        self.win.protocol("WM_DELETE_WINDOW", self.cancelWindow)
        #win.bind('<Destroy>', self.closedWindow) #seems to run on all child widget's destuction too

        wframe = ttk.Labelframe(self.win, text='Set parameters used when writing out the spectrum', padding="10 10 10 10")
        wframe.pack(fill=tk.BOTH, expand=1, padx=4, pady=4)

        #Toggle merging spectral orders for output spectrum
        varMergeOrders = tk.IntVar()
        varMergeOrders.set(self.par.bMergeOrd)
        mergeOrdersF = mergeOrders(varMergeOrders, self.par)
        chkbutMergeOrders = ttk.Checkbutton(wframe, text='Merge spectral\norders',
                        variable=varMergeOrders, command=mergeOrdersF.setFlag)
        ToolTip(chkbutMergeOrders, "Merge spectral orders when writing to file.")
        chkbutMergeOrders.grid(row=0, column=0, columnspan=10, padx=10, pady=4)
        wframe.columnconfigure(0, weight=1)
        wframe.rowconfigure(0, weight=1)
        
        #Set a scaling factor for the output wavelength
        frameScale = ttk.Frame(wframe, padding="0 0 0 0")
        frameScale.grid(row=1, column=0, columnspan=10, sticky=(tk.N, tk.S, tk.E, tk.W))
        wframe.columnconfigure(1, weight=1)
        wframe.rowconfigure(1, weight=1)
        lblScaleWl = ttk.Label(frameScale, text='Scale output wavelength by', justify='right')
        lblScaleWl.grid(row=0, column=0, sticky=tk.E, pady=4, padx=2)
        self.txt_ScaleWl = tk.StringVar()
        self.txt_ScaleWl.set('{:.1f}'.format(self.par.outputWavelengthScale))
        entryScaleWl = ttk.Entry(master=frameScale, textvariable=self.txt_ScaleWl, width=10)
        entryScaleWl.bind('<Key-Return>', self.setWavelengthScale)
        ToolTip(entryScaleWl, 'Scale the output wavelength by some factor.  E.g. 0.1 to convert A to nm', wraplength = 300)
        entryScaleWl.grid(row=0, column=1, sticky=tk.W, pady=4)
        frameScale.columnconfigure(0, weight=1)
        frameScale.columnconfigure(1, weight=1)
        frameScale.rowconfigure(0, weight=1)

        #Optionaly select a vacuum-air conversion for wavelength
        frameAirVac = ttk.Frame(wframe, padding="0 0 0 0")
        frameAirVac.grid(row=2, column=0, columnspan=10, sticky=(tk.N, tk.S, tk.E, tk.W))
        wframe.rowconfigure(2, weight=1)
        lblAirVac = ttk.Label(frameAirVac, text='Convert wavelength between air and vacuum \n'
                                          +'(requires units in A, use the scaling above \n'
                                          +'to convert units)')
        lblAirVac.grid(row=0, column=0, rowspan=3, sticky=tk.E)
        #Set an inital air-vac value
        self.txt_AirVac = tk.StringVar()
        if self.par.outputConvertAirVac < 0:
            self.txt_AirVac.set('air-to-vac')
        elif self.par.outputConvertAirVac > 0:
            self.txt_AirVac.set('vac-to-air')
        else:
            self.txt_AirVac.set('none')
        #Build radio buttons
        convertAirVac = ttk.Radiobutton(frameAirVac, text='Air-to-Vacuum',
            variable=self.txt_AirVac, value='air-to-vac', command=self.airVacConv)
        convertNone = ttk.Radiobutton(frameAirVac, text='None',
            variable=self.txt_AirVac, value='none', command=self.airVacConv)
        convertVacAir = ttk.Radiobutton(frameAirVac, text='Vacuum-to-Air',
            variable=self.txt_AirVac, value='vac-to-air', command=self.airVacConv)
        convertAirVac.grid(row=0, column=1, sticky=tk.W)
        convertNone.grid(row=1, column=1, sticky=tk.W)
        convertVacAir.grid(row=2, column=1, sticky=tk.W)
        frameAirVac.columnconfigure(0, weight=1)
        frameAirVac.columnconfigure(1, weight=1)

        #Buttons to accept or cancle, closing the window either way
        butDone = ttk.Button(master=wframe, text="Done", command=self.doneWindow)
        butCancle = ttk.Button(master=wframe, text="Cancel", command=self.cancelWindow)
        butDone.grid(row=5, column=0, pady=4, sticky='S')
        butCancle.grid(row=5, column=1, pady=4, sticky='S')
        wframe.columnconfigure(1, weight=1)
        wframe.rowconfigure(5, weight=1)
        
        #Optionally set the keyboard focus to a widget
        #self.win.initial_focus = wframe
        #self.win.initial_focus.focus_set()

        #Enter a local event loop, returning only when the given window is destroyed.
        #(Maybe not absolutely necessary with grab_set() above)
        self.win.wait_window(self.win)
        
    def setWavelengthScale(self, *event):
        text = self.txt_ScaleWl.get()
        try:
            waveScale = float(text)
        except ValueError:
            print('got a non-number value for output wavelength scaling: {:}'.format(text))
            messagebox.showwarning("Bad value",
            "Need a number for the wavelength scaling.\n"
                +"You entered:\n{:}".format(text) )
            return False
        if not (waveScale > 0.):
            print('got negative value for output wavelength scaling: {:}'.format(waveScale))
            messagebox.showwarning("Error negative value",
            "The wavelength scaling needs to be positive.\n"
                +"You entered: {:}".format(waveScale) )
            return False
        self.par.outputWavelengthScale = waveScale
        return True
    
    def airVacConv(self, *event):
        airVacdict = {'air-to-vac': -1, 'none': 0, 'vac-to-air': +1}
        self.par.outputConvertAirVac = airVacdict[self.txt_AirVac.get()]
        return
        

#Handle drawing a range selection on the screen, and getting the range selected
class rangeSelect:
    def __init__(self, canvas, parentUIRange):
        self.canvas = canvas
        self.canvasWidget = canvas.get_tk_widget()
        self.parentUIRange = parentUIRange
        self.bStartSelect = True
        self.lastrect = None
        
    def onClick(self, event):
        #Tkinter uses units from the top left, while matplotlib provides units from the bottom left.
        #So we need to flip matplotlilb's provided values (even in screen space)
        x0full, x1full = self.canvas.figure.gca().bbox.intervalx
        y0full, y1full = self.canvas.figure.gca().bbox.intervaly
        height = self.canvas.figure.bbox.height
        y0mod = height - y1full
        y1mod = height - y0full
        
        if self.bStartSelect:  #For the 1st click
            self.x0 = x0 = event.x
            self.y0 = y0 = event.y
            
            #drawn then update a rectangle as the mouse moves
            self.bStartSelect = False
            self.lastrect = self.canvasWidget.create_rectangle(x0, y0mod, x0, y1mod)
            self.oldBindMove = self.canvasWidget.bind('<Motion>')
            self.funcid = self.canvasWidget.bind('<Motion>', self.onMovePointer)
        else:
            self.x1 = event.x
            self.y1 = event.y

            #Get the cursor position in plotted data coordinates
            #(Testing this against matplotlib's toolbar it seems accurate to a pixel level.)
            axesdatax0, axesdatax1, axesdatay0, axesdatay1 = self.canvas.figure.gca().axis()
            axesxrange = x1full - x0full
            dataxrange = axesdatax1 - axesdatax0
            self.datax0 = (self.x0 - x0full)/(axesxrange)*(dataxrange) + axesdatax0
            self.datax1 = (self.x1 - x0full)/(axesxrange)*(dataxrange) + axesdatax0
            
            #Call the creating object's function for processing this range.
            #This is a bit recursive, so there may be a small possibility
            #to leak memory or odd behavior?
            self.parentUIRange.selectedWl(min(self.datax0, self.datax1),
                                          max(self.datax0, self.datax1))
            self.deactivate()

    def cancel(self, *event):
            self.deactivate()
            
    def deactivate(self):
        #if there is an active selection region (not starting to select) undo it
        if not self.bStartSelect:
            self.bStartSelect = True
            self.canvasWidget.delete(self.lastrect)
            #Bind can take a string of tk binding info,
            #and bind retuns such a string if called with just an event name,
            #so we can use that to restore a previous binding.
            self.funcid = self.canvasWidget.bind('<Motion>', self.oldBindMove)

    def onMovePointer(self, event):
        #Delete and redraw the rectangle when the mouse moves
        x1 = event.x
        y1 = event.y
        self.canvasWidget.delete(self.lastrect)
        #Get the full vertical range of the plot for the rectangle
        y0full, y1full = self.canvas.figure.gca().bbox.intervaly
        #Change from matplotlib to tk's y-axis direction
        height = self.canvas.figure.bbox.height
        y0mod = height - y1full
        y1mod = height - y0full
        self.lastrect = self.canvasWidget.create_rectangle(self.x0, y0mod, x1, y1mod, dash=[5,5], width=2)
        #self.lastrect = self.canvasWidget.create_rectangle(self.x0, y0mod, x1, y1mod, dash=[3,3])


#A fairly simple, fairly general tooltip
class ToolTip(object):
    """
    Create a tooltip for a given widget, using Tkinter
    This is a fairly simple class, but it wraps text at a specified length,
    puts the tooltip below the widget, and moves the tooltip to keep it 
    within the screen limits.  
    Based on the recipe from:
    http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml#e387
    and some discussion from: 
    https://stackoverflow.com/questions/3221956/how-do-i-display-tooltips-in-tkinter
    """
    def __init__(self, widget, text='widget info',
                 waittime = 750, wraplength = 180):
        self.waittime = waittime       #miliseconds
        self.wraplength = wraplength   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.waitid = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.waitid = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        waitid = self.waitid
        self.waitid = None
        if waitid:
            self.widget.after_cancel(waitid)

    def showtip(self, event=None):
        #Create a window showing the tooltip
        #shouldn't be necessary, but just in case remove any existing window
        self.hidetip() 
        #Make the tooltip appear outside the mouse position (below the button)
        #helps prevent accidental loops of showing and hiding the tooltip.  
        x = self.widget.winfo_pointerx() + 1
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 1
        #Creates a toplevel window to hold the tooltip
        #and remove all the standard window parts.
        self.tw = tk.Toplevel(self.widget)
        self.tw.wm_overrideredirect(True)
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       wraplength = self.wraplength)
        #extra 1 pixel spacing around the text
        label.pack(ipadx=1) 
        #Check that the tooltip will fit on the screen,
        #and if it won't move the window containing it.
        if (x + label.winfo_reqwidth() + 2) > self.widget.winfo_screenwidth():
            x = self.widget.winfo_screenwidth() - (label.winfo_reqwidth() + 2)
        if (y + label.winfo_reqheight()) > self.widget.winfo_screenheight():
            y = self.widget.winfo_rooty() - label.winfo_reqheight() - 2
        #The new window seems to not draw until update_idletasks is called
        #or when this function exits and control returns to the main loop.
        #So we can change the after everything else.
        self.tw.wm_geometry("+%d+%d" % (x, y))

    def hidetip(self):
        #set self.to to None first, in case the destroy method throws an exception?
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()
