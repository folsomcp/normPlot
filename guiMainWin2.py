#Make the main window for the GUI with Tkinter

import tkinter as tk
import tkinter.ttk as ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
try: #for matplotlib 3.x
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
    mplToolbar=3
except ImportError:  #Probably not in matplotlib 3.x
    try: #for matplotlib 2.x
        from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
        mplToolbar=2
    except ImportError: #not sure what failed, try just not using this extra toolbar
        mplToolbar=0

#For implementing the default Matplotlib key bindings.
#from matplotlib.backend_bases import key_press_handler
import guiControls2 as guic

def makeWin(fig, ax, par, polyDegs, ords, obsWl, obsI, obsSig, obsIavg,
            bFittable, plObs, setPlObsO, setPlPoly, plFitting):
    #Build GUI with tkinter
    root = tk.Tk(className='norm')
    #the className seems to set an icon title
    root.title("Normalize spectra")
    root.iconname("normSpec") #not sure what this does!
    #(iconname only used by some window managers)
    try:
        imgIcon = tk.Image("photo", file="./iconNorm.png")
        root.iconphoto(True, imgIcon)
    except:
        pass

    #For a frame around the window (change padding)
    #mainframe = ttk.Frame(root, padding="3 3 12 12")
    #mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
    #mainframe.pack(fill=tk.BOTH, expand=1)
    #Then widgets in the frame should use it as their parent.

    #Generate a tk.DrawingArea using matplotlib's Tk backend,
    #from the given figure
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw() #this draw call may not be strictly necessary? 
    canvasWidget = canvas.get_tk_widget()
    canvasWidget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas.figure.tight_layout()  #maybe use matplotlib's tighter boarders?
    
    #Add a matplotlib's standard toolbar.
    if mplToolbar > 0:
        #Ideally I might modify this to remove some less usefull buttons?
        if mplToolbar == 2:
            toolbarMpl = NavigationToolbar2TkAgg(canvas, root)
        if mplToolbar == 3:
            toolbarMpl = NavigationToolbar2Tk(canvas, root)
        #Note, making the toolbar calls the parent container's .pack()
        toolbarMpl.update() #maybe not strictly necessary?
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        #this seems to be equivelent, but probably better to call the intended public method
        #canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    #Pass key presses on to matplotlib (may not want this?)
    # uses matplotlib's slower(?) event handling
    # seems to only work when that canvas has the keyboard focus
    #def on_key_press(event):
    #    key_press_handler(event, canvas, toolbarMpl)
    #canvas.mpl_connect("key_press_event", on_key_press)
    
    #Make a frame to hold our tool buttons,
    #and allow us to use use the tk grid geometry manager for them
    tools = ttk.Frame(root, padding="3 3 3 3")
    #tools.pack(fill=tk.BOTH, expand=1)
    tools.pack(side=tk.BOTTOM, fill=tk.X)

    
    #Generate UI for controls
    
    #Create a little event manager-ish object for processing
    #callbacks for view changing buttons
    bView = guic.bViewFuncs(ax, plObs, canvas)
    #and put view changing buttons in their own frame
    viewtools = ttk.Frame(tools)
    viewtools.grid(row=4, column=0, columnspan=8)
    #
    butAuto = ttk.Button(master=viewtools, text="auto",
                         command=bView.autoScale, width=7)
    guic.ToolTip(butAuto, "Auto-scale axes ranges [a]")
    butAuto.grid(row=0, column=0, sticky=(tk.E,tk.W))
    butAutoY = ttk.Button(master=viewtools, text="auto-y",
                          command=bView.autoScaleY, width=7)
    guic.ToolTip(butAutoY, "Auto-scale y-axis range [A]")
    butAutoY.grid(row=0, column=1, sticky=(tk.E,tk.W))
    butZin = ttk.Button(master=viewtools, text="zoom",
                        command=bView.zoomRec, width=5)
    guic.ToolTip(butZin, "Zoom to selected region [z]")
    butZin.grid(row=0, column=2, sticky=(tk.E,tk.W))
    butZin = ttk.Button(master=viewtools, text="Zin",
                        command=bView.zoomIn, width=5)
    guic.ToolTip(butZin, "Zoom in [i]")
    butZin.grid(row=0, column=3, sticky=(tk.E,tk.W))
    butZout = ttk.Button(master=viewtools, text="Zout",
                         command=bView.zoomOut, width=5)
    guic.ToolTip(butZout, "Zoom out [o]")
    butZout.grid(row=0, column=4, sticky=(tk.E,tk.W))
    butLeft = ttk.Button(master=viewtools, text="\u2190",
                         command=bView.panLeft, width=3)
    guic.ToolTip(butLeft, "Pan left [arrow key]")
    butLeft.grid(row=0, column=5, sticky=(tk.E,tk.W))
    butRight = ttk.Button(master=viewtools, text="\u2192",
                          command=bView.panRight, width=3)
    guic.ToolTip(butRight, "Pan right [arrow key]")
    butRight.grid(row=0, column=6, sticky=(tk.E,tk.W))
    butUp = ttk.Button(master=viewtools, text="\u2191",
                       command=bView.panUp, width=3)
    guic.ToolTip(butUp, "Pan up [arrow key]")
    butUp.grid(row=0, column=7, sticky=(tk.E,tk.W))
    butDown = ttk.Button(master=viewtools, text="\u2193",
                         command=bView.panDown, width=3)
    guic.ToolTip(butDown, "Pan down [arrow key]")
    butDown.grid(row=0, column=8, sticky=(tk.E,tk.W))
    #Allow buttons to expand if necessary
    viewtools.rowconfigure(0, weight=1)
    viewtools.columnconfigure(0, weight=1)
    viewtools.columnconfigure(1, weight=1)
    viewtools.columnconfigure(2, weight=1)
    viewtools.columnconfigure(3, weight=1)
    viewtools.columnconfigure(4, weight=1)
    viewtools.columnconfigure(5, weight=1)
    viewtools.columnconfigure(6, weight=1)
    viewtools.columnconfigure(7, weight=1)
    viewtools.columnconfigure(8, weight=1)
    

    #Optional parameters for the spectrum when writen,
    #set in a new window to reduce the amount of clutter.
    winOutPars = guic.newWindowOutPar(root, par)
    butSetParams = ttk.Button(master=tools, text="set output\nparams...",
                              command=winOutPars.openWindow)
    guic.ToolTip(butSetParams, "Set additional parameters for the output spectrum.")
    butSetParams.grid(row=0, column=0)
    
    #Open a new window for seting order polynomial degrees
    winSetPoly = guic.newWindowDeg(root, polyDegs, ords, setPlObsO, setPlPoly)
    butSetPoly = ttk.Button(master=tools, text="set poly.\ndegree...",
                            command=winSetPoly.openWindow)
    guic.ToolTip(butSetPoly, "Set the degrees of the fitting polynomials.")
    butSetPoly.grid(row=0, column=1)

    #Toggle looking ahead to fill in gaps, when there is an exclude region
    # at the end of a spectra order
    varFillEdgeGaps = tk.IntVar()
    varFillEdgeGaps.set(par.lookToNextOrderForGaps)
    fillEdgeGaps = guic.fillEdgeGaps(varFillEdgeGaps, par)
    chkbutFillEdgeGaps = ttk.Checkbutton(tools, text='fill order\nedge gaps',
                         variable=varFillEdgeGaps, command=fillEdgeGaps.setFlag)
    guic.ToolTip(chkbutFillEdgeGaps, 'Look for extra pixels from the next/previous spectral order, when ending an order in an exclude region.', wraplength = 300)
    chkbutFillEdgeGaps.grid(row=0, column=3, padx=2)
    
    #Set the moving/running average length
    lblRunningAvg = ttk.Label(tools, text='average\nlength', justify='right', padding=(2,0,0,0))
    lblRunningAvg.grid(row=0, column=4, sticky=tk.E)
    txt_runningAvg = tk.StringVar()
    txt_runningAvg.set('{:n}'.format(par.averageLen))
    changeRunningAvg = guic.changeRunningAvg(txt_runningAvg, obsI, obsIavg,
                                             ords, par)
    entryRunAvg = ttk.Entry(master=tools, textvariable=txt_runningAvg, width=5)
    entryRunAvg.bind('<Key-Return>', changeRunningAvg.redoAverage)
    guic.ToolTip(entryRunAvg, 'Set the length of the moving/running average applied to the observation, in pixels.', wraplength = 300)
    entryRunAvg.grid(row=0, column=5, sticky=tk.W)
    #don't pass events from this widget to the root widget
    #(possibly overkill for dealing with only a few problematic events!)
    tmpBindTags = entryRunAvg.bindtags()
    entryRunAvg.bindtags((tmpBindTags[0], tmpBindTags[1], tmpBindTags[3]))

    #Set the bin in velocity for selecting best continuum points
    lblVelBin = ttk.Label(tools, text='srch. bin\n(km/s)', justify='right', padding=(2,0,0,0))
    lblVelBin.grid(row=0, column=6, sticky=tk.E)
    txt_velBin = tk.StringVar()
    txt_velBin.set('{:.0f}'.format(par.velBin))
    changeBinSize = guic.changeBinSize(canvas, txt_velBin, obsWl, obsIavg, 
                            obsSig, ords.obsOrder, par, bFittable, plFitting)
    entryBinSize = ttk.Entry(master=tools, textvariable=txt_velBin, width=6)
    entryBinSize.bind('<Key-Return>', changeBinSize.redoBestInBin)
    guic.ToolTip(entryBinSize, 'Set the size of the bin searched for the best continuum point, in km/s.', wraplength = 300)
    entryBinSize.grid(row=0, column=7, sticky=tk.W)
    #don't pass events from this widget to the root widget
    tmpBindTags = entryBinSize.bindtags()
    entryBinSize.bindtags((tmpBindTags[0], tmpBindTags[1], tmpBindTags[3]))    

    #Set regions to be included in fit
    uiIncludeRange = guic.uiIncludeRange(canvas, obsWl, obsI, bFittable,
                                         ords, setPlObsO)
    butIncRange = ttk.Button(master=tools, text='include\nrange',
                             command=uiIncludeRange.runSpanSelect)
    guic.ToolTip(butIncRange, 'Include selected range in the fit.')
    butIncRange.grid(row=0, column=9)

    #Set regions to be excluded in fit
    uiExcludeRange = guic.uiExcludeRange(canvas, obsWl, obsI, bFittable,
                                         ords, setPlObsO)
    butExcRange = ttk.Button(master=tools, text='exclude\nrange',
                             command=uiExcludeRange.runSpanSelect)
    guic.ToolTip(butExcRange, 'Exclude selected range from the fit.')
    butExcRange.grid(row=0, column=10)

    #Link the include and exclude buttons so they can turn eachother off
    uiIncludeRange.linkButton(butIncRange, uiExcludeRange)
    uiExcludeRange.linkButton(butExcRange, uiIncludeRange)
    
    #Save the used fitting parameters for late use
    #setup function used to pass information for writing
    par.linkPolyDegExclude(polyDegs, bFittable, obsWl, ords)
    #helper function that calls par.saveValues, gets info from the Entry boxes
    saveParams = guic.saveParams(par, changeRunningAvg, changeBinSize) 
    butSavePar = ttk.Button(master=tools, text='save\nparams',
                            command=saveParams.doSave)
    guic.ToolTip(butSavePar, 'Save the parameters currently used, including polynomial degrees and included wavelength ranges.', wraplength = 300)
    butSavePar.grid(row=4, column=9)
    
    #Run the continuum fitting 
    bFitCont = guic.runFitCont(canvas, changeRunningAvg, changeBinSize,
                               obsWl, obsI, obsSig, ords, bFittable, obsIavg,
                               par, polyDegs, setPlPoly, plFitting)
    butFitCont = ttk.Button(master=tools, text='fit cont.',
                            command=bFitCont.refitCont)
    guic.ToolTip(butFitCont, 'Fit the continuum and plot the results')
    butFitCont.grid(row=4, column=10, sticky=(tk.N, tk.S, tk.E, tk.W))

    #Empty buffer slot for the grid manager to expand
    #(this may be a bit of a hack!)
    lblBuffer = ttk.Label(tools, text='')
    lblBuffer.grid(row=0, column=8)

    #let .grid() managed cells expand
    tools.columnconfigure(0, weight=1)
    tools.columnconfigure(1, weight=1)
    tools.columnconfigure(2, weight=1)
    tools.columnconfigure(3, weight=1)
    tools.columnconfigure(4, weight=1)
    tools.columnconfigure(5, weight=1)
    tools.columnconfigure(6, weight=1)
    tools.columnconfigure(7, weight=1)
    tools.columnconfigure(8, weight=30)
    tools.columnconfigure(9, weight=1)
    tools.columnconfigure(10, weight=1)
    tools.rowconfigure(0, weight=1)
    tools.rowconfigure(1, weight=3)
    tools.rowconfigure(4, weight=1)

    #Set some general key bindings
    root.bind('<Key-Left>', bView.panLeft)
    root.bind('<Key-Right>', bView.panRight)
    root.bind('<Key-Up>', bView.panUp)
    root.bind('<Key-Down>', bView.panDown)
    root.bind('a', bView.autoScale)
    root.bind('A', bView.autoScaleY)
    root.bind('z', bView.zoomRec)
    root.bind('i', bView.zoomIn)
    root.bind('o', bView.zoomOut)

    ##Not used, but a good quit function
    #def _quit():
    #    root.quit()     # stops mainloop
    #    root.destroy()  # this is necessary on Windows

    #Run the main loop!
    root.mainloop()
