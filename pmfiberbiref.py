#!/usr/bin/python
from __future__ import print_function
import sys
import time

if sys.version_info > ( 3, 0 ):
    is_python3 = True
else:
    is_python3 = False

if is_python3:
    from tkinter import *
    import tkinter.ttk as ttk
    from tkinter import font as tkFont
    import tkinter.filedialog as tkf
    import tkinter.messagebox as tkMessageBox
else:
    from Tkinter import *
    import ttk
    import tkFont
    import tkFileDialog as tkf
    import tkMessageBox

# need these for pyinstaller
import six

import numpy as np
import pylab as plt

# import tooltip

import os
import time
import string
import sys

import version
import config
import bir
import panda_fiber
import bowtie_fiber
import graphs
import ellipses
import help

svFormats = [
	('All Files','*.*'),
	('CSV Files','*.csv'),
	]

# global
logfile = None

NoneType = type(None)

win_main = "Polarization Maintaining Fiber Birefringence/Stress Calculator (version "+version.version+")"

# -------------------------------------------------------------------------
class IORedirector(object):
    '''A general class for redirecting I/O to this Text widget.'''
    def __init__(self,text_area, root=None):
        self.text_area = text_area
        self.root = root
    def flush(self, *args):
        if logfile != None:
            logfile.flush()
        if self.root != None:
            self.root.update()
        return
    def write(self,str):
        # self.text_area.write(str,False)
        self.text_area.insert('end', str)
        self.text_area.see('end')
        if logfile != None:
            logfile.write(str)

class Controller:
    def __init__(self,debug=False, save=False):
        self.root = Tk()
        self.root.title("Controller - "+win_main)
        self.savesnapshot = save
        
        s = ttk.Style()
        s.configure('pink.TButton', background='#FFAAAA')
        s.configure('red.TButton', background='#FF0000')
        s.configure('blue.TButton', background='#AAFFFF')
        s.configure('green.TButton', background='#AAFFAA')
        s.configure('yellow.TButton', background='yellow')
        
        # fonts
        self.titlefont = titlefont = tkFont.Font( family='Helvectica', size=14, weight='bold')
        self.boldfont  = boldfont  = tkFont.Font( family='Helvectica', weight='bold')
        
        # tabs
        nt = self.nt = ttk.Notebook(self.root, padding="12 12 12 12") # main frame
        
        mn_frame = mf = ttk.Frame(nt)
        hp_frame = self.hp_frame  = ttk.Frame(nt)
        te_frame = self.te_frame  = ttk.Frame(nt)
        
        nt.add(mn_frame, text='Parameters') # 0
        nt.add(te_frame, text='Material TEC')        # 1: TEC (thermal expansion coefficient)
        nt.add(hp_frame, text='Help')       # 2
        
        marea = ttk.Frame(self.root, padding="12 3 12 12") # message area
        tarea = ttk.Frame(self.root, padding="12 3 12 12") # console area
        
        self.instructions()
        
        # =========================================
        rownum = 0
        
        # harea = ttk.Frame(mf, padding=(200,0,200,0)) # title
        harea = ttk.Frame(mf, padding=(100,0,100,0)) # title
        harea.grid( column=0, row=rownum, columnspan=2, sticky=(N,W,E,S) )
        rownum += 1
        
        # zarea stacked on top of karea, both inside column varea
        varea = ttk.Frame(mf, relief='solid') # column 0, row 1: geometry parameters
        varea.grid( column=0, row=rownum, sticky=(N,S,W), padx=15,pady=5, ipadx=5,ipady=5 )
        
        zarea = ttk.Frame(varea) # column 0, row 1: geometry parameters
        zarea.grid( column=0, row=0, sticky=(N,W), padx=15,pady=5, ipadx=5,ipady=5 )
        karea = ttk.Frame(varea) # column 0, row 2: material parameters
        karea.grid( column=0, row=1, sticky=(N,W), padx=15,pady=5 )
        
        # carea stacked on top of barea, both inside column v2area;
        v2area = ttk.Frame(mf, relief='solid') # column 0, row 1: geometry parameters
        v2area.grid( column=1, row=rownum, sticky=(N,S,W), padx=15,pady=5, ipadx=5,ipady=5 )
        
        carea = ttk.Frame(v2area) # column 1, row 1: simulation parameters
        carea.grid( column=0, row=0, sticky=(N,W), padx=15,pady=5  )
        barea = ttk.Frame(v2area) # column 1, row 2: buttons
        barea.grid( column=0, row=1, sticky=(N,W), padx=15,pady=5  )
        
        # ---------- title -----------------
        title = ttk.Label(harea, anchor='center',
                          text="Polarization Maintaining Fiber Birefringence/Stress Calculator", 
                          foreground='blue', font=titlefont) # , background='yellow' )
        title.grid( column=0, row=0, padx=0,pady=5, ipadx=0,ipady=5, sticky=(E,W))
        
        # -------------------------------------
        zrow = -1
        
        zrow += 1
        ttk.Label(zarea,text='Fiber Geometry',font=boldfont).\
        grid(column=0, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        self.sap_type_var = StringVar( value = 'panda' )
        pic = ttk.Combobox( zarea, textvariable=self.sap_type_var, state='readonly')
        pic.grid( column=1, row=zrow, columnspan=4, sticky=(W), pady=5 )
        pic.bind('<<ComboboxSelected>>', self.change_N)
        pic['values'] = ('panda','bow-tie')
        
        # headers
        zrow += 1
        ttk.Label(zarea,text='Preform').grid(column=1,row=zrow, ipadx=5,ipady=2)
        ttk.Label(zarea,text='Fiber').grid(column=3,row=zrow, ipadx=5,ipady=2)
        
        # ............................................
        zrow += 1
        ttk.Label(zarea,text='cladding diameter: ').grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.clad_var = DoubleVar( value=30 )
        xp = ttk.Entry(zarea,textvariable=self.clad_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.clad_fb_var = DoubleVar( value=500 )
        xp = ttk.Entry(zarea,textvariable=self.clad_fb_var, width=6)
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        # ............................................
        zrow += 1
        ttk.Label(zarea,text='core diameter: ').grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.core_var = DoubleVar( value=2.5 )
        xp = ttk.Entry(zarea,textvariable=self.core_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.core_fb_var = DoubleVar()
        xp = ttk.Label(zarea,textvariable=self.core_fb_var, width=6, relief='sunken', background='#CCCCCC')
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        # ............................................
        zrow += 1
        ttk.Label(zarea,text='panda rod diameter: ',background='#AAFFFF')\
        .grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.t1_var = DoubleVar( value=6.5 ) # diameter instead of radius!
        xp = ttk.Entry(zarea,textvariable=self.t1_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.t1_fb_var = DoubleVar()
        xp = ttk.Label(zarea,textvariable=self.t1_fb_var, width=6, relief='sunken', background='#CCCCCC')
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        # ............................................
        zrow += 1
        ttk.Label(zarea,text='core to panda rod center: ',background='#AAFFFF')\
        .grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.d2_var = DoubleVar( value=7.0 )
        xp = ttk.Entry(zarea,textvariable=self.d2_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.d2_fb_var = DoubleVar()
        xp = ttk.Label(zarea,textvariable=self.d2_fb_var, width=6, relief='sunken', background='#CCCCCC')
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)

        # ............................................
        zrow += 1
        ttk.Label(zarea,text='bow-tie small radius: ',background='#FFAAFF')\
        .grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.r1_var = DoubleVar( value=5.0 )
        xp = ttk.Entry(zarea,textvariable=self.r1_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.r1_fb_var = DoubleVar()
        xp = ttk.Label(zarea,textvariable=self.r1_fb_var, width=6, relief='sunken', background='#CCCCCC')
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)

        # ............................................
        zrow += 1
        ttk.Label(zarea,text='bow-tie large radius: ',background='#FFAAFF')\
        .grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.r2_var = DoubleVar( value=10.0 )
        xp = ttk.Entry(zarea,textvariable=self.r2_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' mm').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        xp.bind('<Return>', self.preform_to_fiber)
        
        self.r2_fb_var = DoubleVar()
        xp = ttk.Label(zarea,textvariable=self.r2_fb_var, width=6, relief='sunken', background='#CCCCCC')
        xp.grid(column=3, row=zrow, sticky=(E), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' um').grid(column=4, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        # ............................................
        zrow += 1
        ttk.Label(zarea,text='bow-tie half-angle: ',background='#FFAAFF')\
        .grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.theta1_var = DoubleVar( value=35.0 )
        xp = ttk.Entry(zarea,textvariable=self.theta1_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(zarea,text=' deg').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        # --------------------- material parameters -------------------------------------
        zrow = -1
        
        zrow += 1
        ttk.Label(karea,text='Material Parameters',font=boldfont).\
        grid(column=0, row=zrow, columnspan=3, sticky=(W), ipadx=5, ipady=2)

        zrow += 1
        ttk.Label(karea,text='ambient to draw/melt temperature diff.: ').grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.DT_var = DoubleVar( value=2030.0 )
        xp = ttk.Entry(karea,textvariable=self.DT_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' deg C (positive value)').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="Young's modulus: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.E_var = DoubleVar( value=7.8 )
        xp = ttk.Entry(karea,textvariable=self.E_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' x10^(10) N/m^2').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="Poisson's ratio: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.nu_var = DoubleVar( value=0.186 )
        xp = ttk.Entry(karea,textvariable=self.nu_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text='').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="photoelastic coeff. C1: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.C1_var = DoubleVar( value=0.65 )
        xp = ttk.Entry(karea,textvariable=self.C1_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' x10^(-12) m^2/N').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="photoelastic coeff. C2: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.C2_var = DoubleVar( value=4.2 )
        xp = ttk.Entry(karea,textvariable=self.C2_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' x10^(-12) m^2/N').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="core-to-cladding thermal elastic coeff.: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.alpha12_var = DoubleVar( value=0.1 )
        xp = ttk.Entry(karea,textvariable=self.alpha12_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' x10^(-6) 1/deg-C').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        zrow += 1
        ttk.Label(karea,text="SAP-to-cladding thermal elastic coeff.: ").grid(column=0, row=zrow, sticky=(E), ipadx=5, ipady=2)
        self.alpha32_var = DoubleVar( value=1.9 )
        xp = ttk.Entry(karea,textvariable=self.alpha32_var, width=6)
        xp.grid(column=1, row=zrow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(karea,text=' x10^(-6) 1/deg-C').grid(column=2, row=zrow, sticky=(W), ipadx=5, ipady=2)
        
        # -------------------------------------
        # buttons
        brow = -1
        
        brow += 1
        ttk.Label(carea,text='Simulation Parameters and Results',font=boldfont).\
        grid(column=0, row=brow, columnspan=3, sticky=(W), ipadx=5, ipady=2)
        
        brow += 1
        ttk.Label(carea,text="Grid = ").grid(column=0, row=brow, sticky=(E), ipadx=5, ipady=2)
        self.G_var = IntVar( value=101 )
        xp = ttk.Entry(carea,textvariable=self.G_var, width=6)
        xp.grid(column=1, row=brow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(carea,text=' (along x and y)').grid(column=2, row=brow, sticky=(W), ipadx=5, ipady=2)
        
        brow += 1
        ttk.Label(carea,text="N = ").grid(column=0, row=brow, sticky=(E), ipadx=5, ipady=2)
        self.N_var = IntVar( value=10 )
        xp = ttk.Entry(carea,textvariable=self.N_var, width=6)
        xp.grid(column=1, row=brow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(carea,text=' (number of terms in Airy stress function)').grid(column=2, row=brow, sticky=(W), ipadx=5, ipady=2)

        brow += 1
        ttk.Label(carea,text="sub = ").grid(column=0, row=brow, sticky=(E), ipadx=5, ipady=2)
        self.subsample_var = IntVar( value=2 )
        xp = ttk.Entry(carea,textvariable=self.subsample_var, width=6)
        xp.grid(column=1, row=brow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(carea,text=' (subsampling for tensor visualization)').grid(column=2, row=brow, sticky=(W), ipadx=5, ipady=2)

        brow += 1
        self.plot2dxy_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot 2-D X-Y stress distribution', variable=self.plot2dxy_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)

        brow += 1
        self.plot2drt_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot 2-D r-theta stress distribution', variable=self.plot2drt_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)

        brow += 1
        self.plot_el_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot stress tensor visualization', variable=self.plot_el_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)
        
        brow += 1
        self.plotbc_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot boundary stress distribution', variable=self.plotbc_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)
        
        brow += 1
        self.plot_xydn_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot X-Y delta n (unchecked: plot stress)', variable=self.plot_xydn_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)
        
        brow += 1
        self.plot_pdesign_var = BooleanVar( value=True )
        ttk.Checkbutton(carea, text='Plot panda fiber design', variable=self.plot_pdesign_var,
                        onvalue=True, offvalue=False).\
                        grid(column=1,row=brow,columnspan=2, stick=(W), ipadx=5,ipady=2)
        
        brow += 1
        self.zmap_var = StringVar( value='x' )
        pic = ttk.Combobox( carea, textvariable=self.zmap_var, state='readonly', width=8)
        pic.grid( column=1, row=brow, sticky=(W), pady=5 )
        # pic.bind('<<ComboboxSelected>>', self.export)
        pic['values'] = ('x','y','volume')
        
        ttk.Label(carea, text=' (color mapping for tensor visual)').\
        grid( column=2, row=brow, sticky=(W), ipadx=5, ipady=2)
        
        # spacer
        brow += 1
        ttk.Label(carea,text='').grid(column=0, row=brow, sticky=(E), ipadx=5, ipady=2)

        brow += 1
        ttk.Label(carea,text='Pand Fiber Design',font=boldfont).\
        grid(column=0, row=brow, columnspan=3, sticky=(W), ipadx=5, ipady=2)
        
        brow += 1
        ttk.Label(carea,text="r = ").grid(column=0, row=brow, sticky=(E), ipadx=5, ipady=2)
        self.r_var = DoubleVar( value=3.0 )
        xp = ttk.Entry(carea,textvariable=self.r_var, width=6)
        xp.grid(column=1, row=brow, sticky=(W), ipadx=5, ipady=2)
        ttk.Label(carea,text=' mm (preform center to edge of SAP)').grid(column=2, row=brow, sticky=(W), ipadx=5, ipady=2)
        
        # ------------------ buttons -------------------------------------------
        brow = -1
        
        brow += 1
        ttk.Button(barea,text='Load config',command=self.load_config, style='green.TButton').\
        grid( column=0, row=brow, sticky=(E,W),padx=5,pady=2 )
        
        ttk.Button(barea,text='Save config',command=self.save_config, style='green.TButton').\
        grid( column=1, row=brow, sticky=(E,W),padx=5,pady=2 )
                
        brow += 1
        ttk.Button(barea,text='Calculate birefringence',command=self.calc, style='green.TButton').\
        grid( column=0, row=brow, sticky=(E,W),padx=5,pady=2 )
        
        ttk.Button(barea,text='Calculate best config. (r fixed)',command=self.best_config, style='green.TButton').\
        grid( column=1, row=brow, sticky=(E,W),padx=5,pady=2 )
        
        brow += 1
        ttk.Button(barea,text='Export cross-sections',command=self.export1, style='green.TButton').\
        grid( column=0, row=brow, sticky=(E,W),padx=5,pady=2 )
        
        # -----------------------------------
        # status message area
        self.note = ttk.Label( marea, text='Ready', width=50, font=titlefont)
        self.note.grid( column=0, row=0, sticky=(W,E), ipadx=5, ipady=2, pady=5, padx=5 )
        self.note['relief'] = 'sunken'
        self.note['background'] = '#AAFFAA'
        self.note['text'] = "Ready"
        
        ttk.Button(marea,text='Clear message area',command=self.clear_msg, style='yellow.TButton').\
        grid( column=1, row=0, sticky=(E),padx=5,pady=2)
        
        ttk.Button(marea,text='Quit',command=self.quit, style='red.TButton').\
        grid( column=2, row=0, sticky=(E),padx=5,pady=2 )
        
        # =========================================
        # TEC (thermal expansion coefficient)
        self.tec_setup()
        
        # =========================================
        # text/message display
        # self.msgarea = Text(tarea)
        self.msgarea = Text(tarea, height=10)
        # self.msgarea.grid( column=0, row=0, sticky=(W,E) )
        self.msgarea.pack(side=LEFT, fill=BOTH, expand=1 )
        
        sb = ttk.Scrollbar(tarea, orient='vertical', command=self.msgarea.yview )
        # sb.grid( row=0, column=1, sticky=(N,S,W) )
        sb.pack( side=TOP, fill=BOTH, expand=1 )
        
        self.msgarea['yscrollcommand'] = sb.set
        self.msgarea.bind('<Key>', lambda e: 'break') # ignore all key presses
        
        # --------------------------- display ---------------
        nt.pack(side=TOP,expand=1,fill=BOTH)
        marea.pack(side=TOP, expand=1, fill=BOTH)
        tarea.pack(side=LEFT, expand=1, fill=BOTH)
        
        # ---------------------------------------------------
        self.display = graphs.Display(self, win_main)
        
        # redirect stdout and stderr to text area
        sys.stdout = sys.stderr = IORedirector( self.msgarea, root=self.root )
        
        # initialize
        self.disc_stress = None
        self.bc_stress = None
        self.xcut_stress = None
        self.ycut_stress = None
        
        # track first time running
        # self.first_xy = True
        # self.first_rt = True
        
        self.calc_cut_run = False
        self.preform_to_fiber()
        self.calc_tec()
        
        self.root.mainloop()
        
        return
    
    # ---------------------------------------------------------------
    def quit(self, *args):
        # restore stdout before signaling the run thread to exit!
        # it can get stuck in trying to dump to the redirected text message area
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        
        try:
            self.display.quit()
        except:
            pass
        
        self.root.destroy()
        print( "* All done!" )
        # sys.exit() # force quit!
        
        return
    
    # ----------------------------------------------------------
    def instructions(self):
        # msgarea = Text( self.hp_frame, width=140, height=18 )
        # msgarea.grid( column=0, row=0, sticky=(W,E) )

        msgarea = Text( self.hp_frame, height=18 )
        msgarea.pack(side=LEFT, fill=BOTH, expand=1 )
        
        sb = ttk.Scrollbar( self.hp_frame, orient='vertical', command=msgarea.yview )
        # sb.grid( row=0, column=1, sticky=(N,S,W) )
        sb.pack( side=TOP, fill=BOTH, expand=1 )
        
        msgarea['yscrollcommand'] = sb.set
        msgarea.bind('<Key>', lambda e: 'break') # ignore all key presses
        
        msgarea.insert('end', help.usage)
        
        return
    
    # ----------------------------------------------------------
    def clear_msg(self, *args):
        self.msgarea.delete(1.0, 'end')
        return
    
    # ----------------------------------------------------------
    def get_param(self):
        header = "get_param(): "
        
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Calculate parameters"
        sys.stdout.flush()
        
        # materials parameters
        self.param = param = {}
        param['DT'] = self.DT_var.get()
        param['E']  = self.E_var.get()*1e10
        param['nu'] = self.nu_var.get()
        # NOTE: dn = -C * stress; we are scaling to 'units' of delta n, but the
        # sign is going to be opposite!
        param['C1']  = self.C1_var.get()*1e-12
        param['C2']  = self.C2_var.get()*1e-12
        
        param['C'] = param['C2'] - param['C1']
        param['alpha12'] = self.alpha12_var.get()*1e-6
        param['alpha32'] = self.alpha32_var.get()*1e-6
        
        print( header," DT = %.1f deg C" % param['DT'] )
        print( header," E  = %.2g N/m^2" % param['E'] )
        print( header," nu = %.3f" % param['nu'] )
        print( header," C1 = %.3g m^2/N" % param['C1'] )
        print( header," C2 = %.3g m^2/N" % param['C2'] )
        print( header," C  = C2 - C1 = %.3g m^2/N" % param['C'] )
        print( header," alpha12  = %.3g 1/deg-C" % param['alpha12'] )
        print( header," alpha32  = %.3g 1/deg-C" % param['alpha32'] )
        
        # derived
        param['beta1'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha12'] # core
        param['beta3'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha32'] # SAP
        
        # geometry parameters
        b = param['b'] = self.clad_var.get()/2.0 # in um; divide by 2 to get radius
        # everything else is normalized to b!
        param['a']  = 0.5*self.core_var.get()/b # divide by 2 to get radius
        
        print( header," cladding %.2f um" % param['b'] )
        print( header," normalized core radius a = %.3f" % param['a'] )
        
        # we get *all* geometry values, but don't use all over them later
        param['sap_type'] = self.sap_type_var.get()
        param['d1'] = (0.5*self.t1_var.get())/b # diameter to radius
        param['d2'] = self.d2_var.get()/b       # but not d2
        
        param['r1'] = self.r1_var.get()/b # small radius of bow-tie
        param['r2'] = self.r2_var.get()/b # large radius of bow-tie
        param['theta1'] = self.theta1_var.get() # degrees
        
        if param['sap_type'] == 'panda':
            print( header," normalized SAP radius d1 = %.3f" % param['d1'] )
            print( header," normalized SAP center to fiber center d2 = %.3f" % param['d2'] )
        elif param['sap_type'] == 'bow-tie':
            print( header," normalized SAP radius r1 = %.3f" % param['r1'] )
            print( header," normalized SAP radius r2 = %.3f" % param['r2'] )
            print( header," normalized SAP half-angle = %.3f deg " % param['theta1'] )
            param['theta1'] *= np.pi/180.0 # degrees to radian
        
        # simulation set up
        param['N'] = self.N_var.get()
        param['gridsize'] = self.G_var.get()
        param['subsample'] = self.subsample_var.get()
        
        return
    
    # ----------------------------------------------------------
    def calc(self, *args):
        header = "calc(): "
        
        self.get_param()
        
        if self.plot_pdesign_var.get():
            self.calc_design()
        
        if self.plotbc_var.get() == True:
            self.calc_bc()
        
        self.calc_cut()
        
        if self.plot2drt_var.get() == True or \
        self.plot2dxy_var.get() == True or \
        self.plot_el_var.get():
            self.calc_disc()
        
        # ----------------------------------------------------------
        # all done!
        self.note['background'] = '#AAFFAA'
        self.note['text'] = "Ready"
        sys.stdout.flush()
        
        return
    
    # design --------------------------------------------------------------------
    def calc_design(self):
        header = "\tcalc_design(): "
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Calculate optimum stress at core"
        print( header," Calculate stress design..." )
        sys.stdout.flush()
        
        self.bir = bir.Design(self.param, ax=self.ax_design)
        
        self.bir.plot()
        self.display.tabs.select( self.display.tab_design )
        
        return
    
    # boundary condition --------------------------------------------------------
    def calc_bc(self):
        header = "\tcalc_bc(): "
        th = np.arange(-np.pi, np.pi, np.pi/100)
        xrim = np.cos(th)
        yrim = np.sin(th)
        
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Calculate to check boundary condition"
        print( header," Calculate BC matching..." )
        sys.stdout.flush()
        
        # if self.bc_stress != None:
        #    del self.bc_stress
        
        if self.param['sap_type'] == 'bow-tie':
            self.bc_stress = bowtie_fiber.Stress(xrim,yrim, u=th/np.pi, param=self.param, blank_outside=False)
        else:
            self.bc_stress = panda_fiber.Stress(xrim,yrim, u=th/np.pi, param=self.param, blank_outside=False)
        
        # plotting BC
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Plotting boundary"
        print( header," plot boundary" )
        sys.stdout.flush()
        self.bc_stress.plot_bc( ax=self.ax_bc )
        self.display.tabs.select( self.display.tab_bc )
        
        print( header," BC error s_r ",  np.abs(self.bc_stress.sr).max() ) 
        print( header," BC error tt_r ", np.abs(self.bc_stress.tt).max() )
        sys.stdout.flush()
        
        return
    
    # calculate disc stress ---------------------------------------------------------------
    def calc_disc(self):
        header = "\tcalc_disc(): "
        
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Calculate stress on disc"
        print( header," Calculate disc...")
        sys.stdout.flush()
        
        # if self.disc_stress != None:
        #    del self.disc_stress
        
        # clean up previous first!
        if self.disc_stress != None:
            self.ax_rstress.cla()
            self.ax_tstress.cla()
            self.ax_xstress.cla()
            self.ax_ystress.cla()
            self.disc_stress.clear_xy_colorbars()
            self.disc_stress.clear_rt_colorbars()
        
        rext = 1.1
        Ns = self.param['gridsize']
        x, y = np.mgrid[ -rext:rext:1j*Ns, -rext:rext:1j*Ns ]
        
        # save snapshot
        # self.disc_stress = panda_fiber.Stress(x,y, param=self.param, save=True)
        if self.param['sap_type'] == 'bow-tie':
            self.disc_stress = bowtie_fiber.Stress(x,y, param=self.param, save=self.savesnapshot)
        else:
            self.disc_stress = panda_fiber.Stress(x,y, param=self.param, save=self.savesnapshot)
        
        print( header," birefringence range ",self.disc_stress.vmin," to ",self.disc_stress.vmax )
        
        # plotting r-theta
        if self.plot2drt_var.get() == True:
            
            self.note['background'] = '#FFAAAA'
            self.note['text'] = "Plotting in polar-coordinates (please wait....)"
            print( header," plot r-theta stress")
            sys.stdout.flush()
            
            self.disc_stress.plot_rt( ax1=self.ax_rstress,
                                     ax2=self.ax_tstress )
            self.display.tabs.select( self.display.tab_rtheta )
        
        # plotting X-Y
        if self.plot2dxy_var.get() == True:
            self.note['background'] = '#FFAAAA'
            self.note['text'] = "Plotting in X-Y coordinates (please wait...)"
            print( header," plot x-y stress" )
            sys.stdout.flush()
            
            self.disc_stress.plot_xy( ax1=self.ax_xstress, ax2=self.ax_ystress, ax3=self.ax_xystress,
                                     plotdn=self.plot_xydn_var.get() )
            self.display.tabs.select( self.display.tab_xy )
            
        # plotting stress tensor
        if self.plot_el_var.get() == True:
            self.note['background'] = '#FFAAAA'
            self.note['text'] = "Plotting stress tensor visualization (please wait...)"
            print( header," plot stress tensor" )
            sys.stdout.flush()
            
            sx  = self.disc_stress.sx
            sy  = self.disc_stress.sy
            txy = self.disc_stress.txy
            sub = self.subsample_var.get()
            color = self.zmap_var.get()
            
            ellipses.plot( sx,sy,txy, ax=self.ax_el, subsample=sub, color=color )
            self.display.tabs.select( self.display.tab_el )
        
        return
    
    # calculate X-Y cross-section stress ---------------------------------------------------------------
    def calc_cut(self):
        header = "\tcalc_cut(): "
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Calculate stress along X-Y cross-sections"
        print( header," Calculate stress along X-Y cross-sections..." )
        sys.stdout.flush()
        
        # parameters for calculating dn;
        # we divide by C because all stresses are already normalized by C
        param = self.param
        AA = (param['C1'] + param['nu'] * param['C2']) / param['C']
        BB = (param['C2'] * ( 1.0 + param['nu'] )) / param['C']
        
        # x-cross: with SAP
        u  = np.arange(-100,101, 1)*0.01
        zr = np.zeros_like(u)
        if self.param['sap_type'] == 'bow-tie':
            self.xcut_stress = bowtie_fiber.Stress(u,zr, param=self.param)
        else:
            self.xcut_stress = panda_fiber.Stress(u,zr, param=self.param)
        
        self.cut_xlist = u
        
        # y-cross: with SAP
        if self.param['sap_type'] == 'bow-tie':
            self.ycut_stress = bowtie_fiber.Stress(zr,u, param=self.param)
        else:
            self.ycut_stress = panda_fiber.Stress(zr,u, param=self.param)
        
        # without SAP --------------------------------------------------------------
        param2 = self.param.copy()
        param2['beta3'] = 0.0
        
        # x-cross: without SAP
        if self.param['sap_type'] == 'bow-tie':
            self.xcut_nosap = bowtie_fiber.Stress(u,zr, param=param2)
        else:
            self.xcut_nosap = panda_fiber.Stress(u,zr, param=param2)
            
        # y-cross: without SAP (by symmetry this should be identical to xcut_nosap, with x and y axes exchanged)
        # self.ycut_nosap = self.xcut_nosap
        
        if self.plot_xydn_var.get(): # plot delta n instead of stress
            # .... x-axis .....
            lsx = self.xcut_stress.sx
            lsy = self.xcut_stress.sy
            XLINE = dnx = -AA * lsx - BB * lsy
            
            # .... y-axis .....
            lsx = self.ycut_stress.sx
            lsy = self.ycut_stress.sy
            YLINE = dny = -BB * lsx - AA * lsy
            
            # .... x-axis (no sap).....
            lsx = self.xcut_nosap.sx
            lsy = self.xcut_nosap.sy
            BASE = dnx_nosap = -AA * lsx - BB * lsy
            
            title1 = r'$\Delta n \times 10^{-4}$ along X-Y'
            label1 = r'$\Delta n_x$ along X-axis'
            label2 = r'$\Delta n_y$ along Y-axis'
        
        else:
            XLINE = self.xcut_stress.sx
            YLINE = self.ycut_stress.sy
            BASE  = self.xcut_nosap.sx
            title1 = r'Stress scaled to $\Delta n \times 10^{-4}$ along X-Y'
            label1 = r'$\sigma_x$ along X-axis (scaled)'
            label2 = r'$\sigma_y$ along Y-axis (scaled)'
        
        XLINE_NOSAP = XLINE - BASE
        YLINE_NOSAP = YLINE - BASE
        
        # plot ----------------------------------------------------------------------
        self.note['background'] = '#FFAAAA'
        self.note['text'] = "Plot stress along X-Y cross-sections"
        print( header," Plot stress along X-Y cross-sections...")
        sys.stdout.flush()
        
        self.ax_cut1.cla()
        self.ax_cut1.set_title(title1)
        self.ax_cut1.set_xlabel('radius (normalized to cladding radius)')
        self.ax_cut1.set_ylabel(r'$\Delta n \times 10^{-4}$')
        
        self.ax_cut1.plot( u, XLINE, label='X-axis cross-section' )
        self.ax_cut1.plot( u, YLINE, label='Y-axis cross-section' )
        
        self.ax_cut1.plot( u, BASE, label='X/Y-axis cross-section (without SAP)' )
        # self.ax_cut1.plot( u, ycut_nosap.sy, label='Y-axis cross-section (without SAP)' )
        
        self.ax_cut1.legend( loc='upper right' )
        self.ax_cut1.figure.canvas.draw()
        
        # plot difference (birefringence), but only in neighborhood of
        # core (it doesn't make sense extending to the stress rods)
        ww = np.where( np.abs(u) < 2.0*self.param['a'] )
        
        # Following Eq. (21) of Chu and Sammut
        biref = self.xcut_stress.sx[ww] - self.ycut_stress.sy[ww]
        
        self.ax_cut2.cla()
        
        self.ax_cut2.set_title(r'Change due to SAP')
        self.ax_cut2.set_xlabel('radius (normalized to cladding radius)')
        self.ax_cut2.set_ylabel(r'$\Delta n \times 10^{-4}$')
        
        self.ax_cut2.plot( u, XLINE_NOSAP, label=label1 )
        self.ax_cut2.plot( u, YLINE_NOSAP, label=label2 )
        self.ax_cut2.plot( u[ww], biref, linewidth=2, 
                          label=r'Birefringence ($C(\sigma_x-\sigma_y)$)' )
        
        # add core boundary lines
        self.ax_cut2.axvline( self.param['a'], linestyle='--', color='#888888')
        self.ax_cut2.axvline( -self.param['a'], linestyle='--', color='#888888')
        
        self.ax_cut2.legend( loc='lower right' )
        self.ax_cut2.figure.canvas.draw()
        
        # select tab
        self.display.tabs.select( self.display.tab_cut )
        
        self.calc_cut_run = True
        
        return
    
    # ------------------------------------------------------------------
    def best_config(self, *args):
        header = "best_config(): "
        
        self.get_param()
        self.calc_design()
        
        # r is distance from the fiber center to the edge of the SAP
        r = self.r_var.get()/self.param['b']                                        
        tb = 2.0 * self.bir.best( r ) # tb is now diameter!
        dia = tb * self.param['b']    # convert to mm
        # print( header," DEBUG ",tb, dia )
        print( header," For this core/cladding ratio of %.3f" % self.param['a'] )
        print( header," and r/a of %.3f (fiber center to the edge of the SAP)" % (r/self.param['a']) )
        print( header," the best SAP diameter is %.2f mm (preform scale)" % dia )
        
        self.note['background'] = '#AAFFAA'
        self.note['text'] = "Ready"
        sys.stdout.flush()
        
        # pop-up
        msg  = "For this core/cladding ratio of %.3f\n" % self.param['a']
        msg += "and r/a of %.3f (fiber center to the\n" % (r/self.param['a'])
        msg += "edge of the SAP) the best SAP diameter is %.2f mm (preform scale)\n" % dia
        
        tkMessageBox.showinfo("Panda Fiber Design", msg)
        
        return
    
    # ------------------------------------------------------------------
    def save_config(self, *args):
        header = "save_config(): "
        
        print( header," saving...")
        sys.stdout.flush()
        
        # read parameters from GUI first!
        self.get_param()
        
        config.save(self.param)
        
        print( header," done!"  )
        sys.stdout.flush()
        
        return

    # ------------------------------------------------------------------
    def load_config(self, *args):
        header = "load_config(): "
        newparam = config.get()
        if newparam == None:
            return
        
        print( header," Load configuration from file" )
        sys.stdout.flush()
        
        param = self.param = newparam
        
        # populate
        self.DT_var.set( param['DT'] )
        
        tmp = "%.2f" % (param['E']*1e-10)
        self.E_var.set( float(tmp) )
        self.nu_var.set( param['nu'] )
        self.C1_var.set( param['C1']*1e12 )
        self.C2_var.set( param['C2']*1e12 )
        
        self.alpha12_var.set( param['alpha12']*1e6 )
        self.alpha32_var.set( param['alpha32']*1e6 )
        
        b = param['b'] # radius
        self.clad_var.set( param['b']*2.0 )   # to diameter
        self.core_var.set( param['a']*b*2.0 ) # to diameter
        self.t1_var.set( param['d1']*b*2.0 )  # to diameter
        self.d2_var.set( param['d2']*b ) # fiber center to SAP center
        
        self.N_var.set( param['N'] )
        self.G_var.set( param['gridsize'] )
        self.subsample_var.set( param['subsample'] )
        
        self.preform_to_fiber() # rescale
        
        # clear all graphs
        self.display.clear_all()
        
        print( header," done!" )
        sys.stdout.flush()
        
        return
    
    # ------------------------------------------------------------------
    def export1(self, *args):
        pre = "export1(): "
        
        if self.calc_cut_run == False:
            self.get_param()
            self.calc_cut()
        
        filename = "fiber_birefringence.csv"
        filename = tkf.asksaveasfilename(filetypes=svFormats, title='Save to file', initialfile=filename)
        if len(filename) <=0:
            print( pre,"! Canceled" )
            return
        
        # gather parameters
        param = self.param
        header = "Parameters:"
        header += " DT = %.1f deg C\n" % param['DT']
        header += " E  = %.2g N/m^2\n" % param['E']
        header += " nu = %.3f\n" % param['nu']
        header += " C1 = %.3g m^2/N\n" % param['C1']
        header += " C2 = %.3g m^2/N\n" % param['C2']
        header += " core-to-cladding thermal elastic coeff. diff. = %.3g 1/deg-C\n" % param['alpha12']
        header += " SAP-to-cladding thermal elastic coeff. diff.  = %.3g 1/deg-C\n" % param['alpha32']
        header += " cladding diameter = %.2f mm\n" % (param['b']*2.0)
        header += " core diameter     = %.3f mm\n" % (param['a']*2.0*param['b'])
        header += " SAP diameter      = %.3f mm\n" % (param['d1']*2.0*param['b'])
        header += " SAP center to fiber center = %.3f mm\n" % (param['d2']*param['b'])
        
        # column names
        header += "radius,x cross-section,y cross-section,x/y-axis cross-section without SAP,"
        header += "dn along x, dn along y,birefringence (x-y)"
        # no \n for last line
        
        AA = (param['C1'] + param['nu'] * param['C2']) / param['C']
        BB = (param['C2'] * ( 1.0 + param['nu'] )) / param['C']
        
        # .... x-axis .....
        lsx = self.xcut_stress.sx - self.xcut_nosap.sx
        lsy = self.xcut_stress.sy - self.xcut_nosap.sy
        dnx = -AA * lsx - BB * lsy
        
        # .... y-axis .....
        # not a typo: we did not calculate self.ycut_nosap, but it
        # will be the same as self.xcut_nosap because of symmetry
        lsx = self.ycut_stress.sx - self.xcut_nosap.sx
        lsy = self.ycut_stress.sy - self.xcut_nosap.sy
        dny = -BB * lsx - AA * lsy
        
        biref = -self.xcut_stress.sx + self.ycut_stress.sy

        data = np.column_stack( (self.cut_xlist, -self.xcut_stress.sx,
                                 -self.ycut_stress.sy, -self.xcut_nosap.sx, 
                                 dnx,dny,biref) )
        fmt = '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f'
        np.savetxt(filename, data, delimiter=',',
                   fmt=fmt, header=header)
        print( pre," data saved to ",filename )
        
        return
    
    # ----------------------------------------------------------
    def preform_to_fiber(self, *args):
        b    = self.clad_var.get() # preform
        b_fb = self.clad_fb_var.get() # fiber
        scale = b_fb/b
        
        core_fb = "%.2f" % (self.core_var.get() * scale)
        self.core_fb_var.set( float(core_fb) )
        
        t1_fb = "%.2f" % (self.t1_var.get() * scale)
        self.t1_fb_var.set( float(t1_fb) )
        
        d2_fb = "%.2f" % (self.d2_var.get() * scale)
        self.d2_fb_var.set( float(d2_fb) )
        
        r1_fb = "%.2f" % (self.r1_var.get() * scale )
        self.r1_fb_var.set( float(r1_fb) )
        
        r2_fb = "%.2f" % (self.r2_var.get() * scale )
        self.r2_fb_var.set( float(r2_fb) )
        
        return
    
    # ----------------------------------------------------------
    # set up GUI for estimating TEC (thermal expansion coefficient)
    def tec_setup(self):
        te_frame = self.te_frame
        titlefont = self.titlefont
        boldfont  = self.boldfont
        
        rownum = 0
        
        h2area = ttk.Frame(te_frame, padding=(200,0,200,0)) # title
        h2area.grid( column=0, row=rownum, sticky=(N,W,E,S) )
        
        rownum += 1
        warea = ttk.Frame(te_frame, relief='solid')
        warea.grid(column=0,row=rownum,sticky=(N,W,E,S), padx=15,pady=5, ipadx=15,ipady=10)
        
        xx = ttk.Label(h2area,text='Thermal Expansion Coefficient of SAP',font=titlefont,foreground='blue')
        xx.grid(column=0, row=0,sticky=(N,S,E,W),ipadx=15, ipady=5, padx=10, pady=10)
        
        wrow = -1
        
        wrow += 1 # ......................................................
        ttk.Label(warea,text='Material').\
        grid(column=0,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Label(warea,text='TEC\n[x10^(-6)/deg-C]').\
        grid(column=1,row=wrow, sticky=(E), padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Label(warea,text='Molar Conc.\n[%]').\
        grid(column=2,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Label(warea,text='Doped TEC\n[x10^(-6)/deg-C]').\
        grid(column=3,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        wrow += 1 # ......................................................
        ttk.Label(warea,text='Silica', font=boldfont).\
        grid(column=0,row=wrow, padx=5,pady=2,ipadx=5,ipady=2, sticky=(W))
        
        self.tec_si_var = DoubleVar( value=0.5 )
        xx = ttk.Entry(warea,textvariable = self.tec_si_var, width=6)
        xx.grid(column=1,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        wrow += 1 # ......................................................
        # ttk.Label(warea,text='B2O3', font=boldfont).\
        # grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
        xx = Text(warea,width=5,height=1, borderwidth=0, font=boldfont, background='#DDDDDD')
        xx.tag_configure("subscript", offset=-4)
        xx.insert("insert", "B","","2","subscript","O","","3","subscript")
        xx.configure(state='disabled')
        xx.grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
        
        self.tec_b2o3_var = DoubleVar( value=10 )
        xx = ttk.Entry(warea,textvariable = self.tec_b2o3_var, width=6)
        xx.grid(column=1,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        # default 20% boron dopant
        self.dop_b2o3_var = DoubleVar( value=20 )
        xx = ttk.Entry(warea,textvariable = self.dop_b2o3_var, width=6)
        xx.grid(column=2,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        self.mix_b2o3_var = DoubleVar()
        ttk.Label(warea,textvariable = self.mix_b2o3_var, width=6, relief='sunken', background="#CCCCCC").\
        grid(column=3,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as SAP',command=lambda: self.use_as_sap('b2o3'), style='green.TButton').\
        grid(column=4,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as core',command=lambda: self.use_as_core('b2o3'), style='green.TButton').\
        grid(column=5,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        wrow += 1 # ......................................................
        # ttk.Label(warea,text='GeO2', font=boldfont).\
        # grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
        xx = Text(warea,width=5,height=1, borderwidth=0, font=boldfont, background='#DDDDDD')
        xx.tag_configure("subscript", offset=-4)
        xx.insert("insert", "GeO","","2","subscript")
        xx.configure(state='disabled')
        xx.grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
        
        self.tec_geo2_var = DoubleVar( value=7 )
        xx = ttk.Entry(warea,textvariable = self.tec_geo2_var, width=6)
        xx.grid(column=1,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        self.dop_geo2_var = DoubleVar( value=3 )
        xx = ttk.Entry(warea,textvariable = self.dop_geo2_var, width=6)
        xx.grid(column=2,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        self.mix_geo2_var = DoubleVar()
        ttk.Label(warea,textvariable = self.mix_geo2_var, width=6, relief='sunken', background="#CCCCCC").\
        grid(column=3,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as SAP',command=lambda: self.use_as_sap('geo2'), style='green.TButton').\
        grid(column=4,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as core',command=lambda: self.use_as_core('geo2'), style='green.TButton').\
        grid(column=5,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        wrow += 1 # ......................................................
        # ttk.Label(warea,text='P2O5', font=boldfont).\
        # grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
        xx = Text(warea,width=5,height=1, borderwidth=0, font=boldfont, background='#DDDDDD')
        xx.tag_configure("subscript", offset=-4)
        xx.insert("insert", "P","","2","subscript","O","","5","subscript")
        xx.configure(state='disabled')
        xx.grid(column=0,row=wrow, sticky=(W), padx=5,pady=2,ipadx=5,ipady=2)
                
        self.tec_p2o5_var = DoubleVar( value=14 )
        xx = ttk.Entry(warea,textvariable = self.tec_p2o5_var, width=6)
        xx.grid(column=1,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        self.dop_p2o5_var = DoubleVar( value=3 )
        xx = ttk.Entry(warea,textvariable = self.dop_p2o5_var, width=6)
        xx.grid(column=2,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        xx.bind('<Return>', self.calc_tec)
        
        self.mix_p2o5_var = DoubleVar()
        ttk.Label(warea,textvariable = self.mix_p2o5_var, width=6, relief='sunken', background="#CCCCCC").\
        grid(column=3,row=wrow, sticky=(E,W), padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as SAP',command=lambda: self.use_as_sap('p2o5'), style='green.TButton').\
        grid(column=4,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        ttk.Button(warea,text='Use as core',command=lambda: self.use_as_core('p2o5'), style='green.TButton').\
        grid(column=5,row=wrow, padx=5,pady=2,ipadx=5,ipady=2)
        
        return
    
    # ----------------------------------------------------------
    def calc_tec(self, *args):
        # calculate TEC
        # header = "calc_tec(): "
        # print( header," Estimate thermal expansion coefficient" )
        
        s0 = self.tec_si_var.get()
        
        # GeO2
        sd = self.tec_geo2_var.get()
        d  = self.dop_geo2_var.get()/100.0 # percent -> fractional
        x  = d*sd + (1-d)*s0
        self.mix_geo2_var.set( float("%.4f" % x ) )
        
        # P2O5
        sd = self.tec_p2o5_var.get()
        d  = self.dop_p2o5_var.get()/100.0 # percent -> fractional
        x  = d*sd + (1-d)*s0
        self.mix_p2o5_var.set( float("%.4f" % x ) )
        
        # b2o3
        sd = self.tec_b2o3_var.get()
        d  = self.dop_b2o3_var.get()/100.0 # percent -> fractional
        x  = d*sd + (1-d)*s0
        self.mix_b2o3_var.set( float("%.4f" % x ) )
        
        return
    
    # ----------------------------------------
    def get_tec(self, material):
        header = "get_tec(): "
        # print( header," ",material )
        self.calc_tec()
        
        if material == 'b2o3':
            x = self.mix_b2o3_var.get()
        elif material == 'geo2':
            x = self.mix_geo2_var.get()
        elif material == 'p2o5':
            x = self.mix_p2o5_var.get()
        else:
            print( header,"Un-recognized material ",material )
            return None
        
        si = self.tec_si_var.get()
        
        # return difference
        return x - si
    
    # ----------------------------------------------------------
    def use_as_sap(self, material):
        x = self.get_tec(material)
        if x == None:
            return
        
        self.alpha32_var.set( float( "%.5f" % x) )
        self.nt.select( 0 )
        
        return
    # ----------------------------------------------------------
    def use_as_core(self, material):
        x = self.get_tec(material)
        if x == None:
            return
        
        self.alpha12_var.set( float( "%.5f" % x) )
        self.nt.select( 0 )
        
        return
    
    # ----------------------------------------------
    def change_N(self, *args):
        if self.sap_type_var.get() == 'panda':
            self.N_var.set(20)
        else:
            self.N_var.set(150)
        
        msg = "Note: changing number of terms to truncate in series"
        tkMessageBox.showinfo("N value change", msg)
        
        return
        
# =================================================================================
if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1]=='-h':
        print( "USAGE: ",sys.argv[0]," [save]" )
        print( "\tsave snapshot if 'save' is set to 1" )
        sys.exit()
    
    if len(sys.argv)>=2:
        save = (sys.argv[1] == '1')
    else:
        save = False
    
    # create log file
    print( "* Create/open log file" )
    try:
        logfile = open("logfile.txt","w")
    except IOError:
        print( "Error creating log file logfile.txt" )
        logfile = None
    
    print( "* Initializing; please wait...." )
    control = Controller(save=save)
    
    print( "* Main loop closed. Bye!" )
    
    sys.exit()
    
    
    
