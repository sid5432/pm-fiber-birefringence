#!/usr/bin/python
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import pylab as plt

# python 2
import Tkinter as Tk

from Tkinter import *
import ttk

class Display:
    def __init__(self, controller, win_main=""):
        self.controller = controller
        
        self.root = Tk()
        self.root.title("Graphs - "+win_main)
        self.tabs = n = ttk.Notebook( self.root, width=1400, height=650 )
        n.grid(column=0, row=0, padx=5, pady=5)
        
        frm = []
        
        idx = -1
        
        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="Cross Section")
        self.tab_cut = idx
        
        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="Pand Design")
        self.tab_design = idx
        
        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="Boundary Condition")
        self.tab_bc = idx
        
        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="X-Y Delta n or Stress")
        self.tab_xy = idx
        
        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="R-theta Stress")
        self.tab_rtheta = idx

        idx += 1
        frm.append( ttk.Frame(n) )
        n.add(frm[idx], text="Stress Tensor Visualization")
        self.tab_el = idx

        # -----------------------------------------------------------------
        idx = -1
        fig = []
        canvas = []
        
        # frame #1: X-Y stress
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_xy = fig[idx]
        controller.ax_xstress = fig[idx].add_subplot(131, aspect='equal')
        controller.ax_ystress = fig[idx].add_subplot(132, aspect='equal')
        controller.ax_xystress = fig[idx].add_subplot(133, aspect='equal')
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_xy]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_xy] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------------------------------------------------------
        # frame #2: R-theta stress
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_rtheta = fig[idx]
        controller.ax_rstress = fig[idx].add_subplot(121, aspect='equal')
        controller.ax_tstress = fig[idx].add_subplot(122, aspect='equal')
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_rtheta]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_rtheta] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------------------------------------------------------
        # frame #3: boundary condition
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_bc = fig[idx]
        controller.ax_bc = fig[idx].add_subplot(111)
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_bc]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_bc] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------------------------------------------------------
        # frame #4: cross-section 
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_cut = fig[idx]
        controller.ax_cut1 = fig[idx].add_subplot(121)
        controller.ax_cut2 = fig[idx].add_subplot(122)
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_cut]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_cut] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------------------------------------------------------
        # frame #5: design
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_design = fig[idx]
        controller.ax_design = fig[idx].add_subplot(111)
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_design]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_design] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------------------------------------------------------
        # frame #6: stress tensor visualization
        idx += 1
        
        fig.append( Figure() )
        # naming....
        controller.fig_el = fig[idx]
        controller.ax_el = fig[idx].add_subplot(111, aspect='equal')
        
        canvas.append( FigureCanvasTkAgg(fig[idx], master=frm[self.tab_el]) )
        canvas[idx].show()
        
        canvas[idx].get_tk_widget().pack( side=TOP )
        
        NavigationToolbar2TkAgg( canvas[idx], frm[self.tab_el] ).update()
        canvas[idx]._tkcanvas.pack(side=TOP, fill=BOTH, expand=1 )
        
        # -----------------
        n.pack(fill=BOTH, expand=1)
        
        # self.root.mainloop()
        
        # self.init_fill() # doesn't work
        
        return
    
    # -------------------------------------------
    def quit(self):
        self.root.destroy()
        
        return
    
    # -------------------------------------------
    def clear_all(self):
        """clear all graphs"""
        header = "clear_all(): "
        
        print header," x-y plot"
        self.controller.ax_xstress.cla()
        self.controller.ax_ystress.cla()
        self.controller.ax_xystress.cla()
        self.controller.ax_xstress.figure.canvas.draw()
        
        print header," r-theta plot"
        self.controller.ax_rstress.cla()
        self.controller.ax_tstress.cla()
        self.controller.ax_rstress.figure.canvas.draw()
        
        print header," BC plot"
        self.controller.ax_bc.cla()
        self.controller.ax_bc.figure.canvas.draw()
        
        print header," cross-section plot"
        self.controller.ax_cut1.cla()
        self.controller.ax_cut2.cla()
        self.controller.ax_cut1.figure.canvas.draw()
        
        print header," design plot"
        self.controller.ax_design.cla()
        self.controller.ax_design.figure.canvas.draw()
        
        print header," tensor visual plot"
        self.controller.ax_el.cla()
        self.controller.ax_el.figure.canvas.draw()
                
        return
    
    # ----------------------------------------------------
    # plotting the initial flat z profile fixes the glitch in removing
    # the colorbar; otherwise the first time the profiles are plotted
    # it is a little smaller (thereafter, subsequent plots are all the
    # same size, but slightly larger).
    def init_fill(self):
        rext = 1.0
        Ns = 51
        x, y = np.mgrid[ -rext:rext:1j*Ns, -rext:rext:1j*Ns ]
        z = np.zeros_like(x)
        
        self.controller.ax_xstress.pcolor(x,y,z, cmap=plt.cm.coolwarm)
        self.controller.ax_ystress.pcolor(x,y,z, cmap=plt.cm.coolwarm)
        self.controller.ax_xystress.pcolor(x,y,z, cmap=plt.cm.coolwarm)
        
        self.controller.ax_rstress.pcolor(x,y,z, cmap=plt.cm.coolwarm)
        self.controller.ax_tstress.pcolor(x,y,z, cmap=plt.cm.coolwarm)
        
        return
    
# =======================
if __name__ == '__main__':
    import sys
    
    class Controller: pass
    controller = Controller()
    
    # main
    root = Tk()
    # create display
    display = Display(controller, win_main='TEST WINDOW')
    
    def quit():
        try:
            display.quit()
        except:
            pass
        root.destroy()
        return
    
    def focus():
        display.tabs.select( display.tab_bc )
        return

    mf = ttk.Frame(root, padding="3 3 12 12")
    mf.grid( column=0, row =0, sticky=( N, W,E, S))
    
    carea = ttk.Frame(mf)
    carea.grid( column=0, row=0, sticky=(N,S) )
    
    ttk.Button(carea, text='Quit', command=quit).grid( column=0, row=0, sticky=(E,W))
    ttk.Button(carea, text='focus to BC', command=focus).grid( column=0, row=1, sticky=(E,W))
    ttk.Button(carea, text='clear all', command=display.clear_all).grid( column=0, row=2, sticky=(E,W))
    print "MAIN: got here!"
    
    # --------------------------- display ---------------
    for child in mf.winfo_children():
        # print "child ", child
        child.grid_configure(padx=5, pady=5)
        
    root.mainloop()
    print "MAIN: done!"
    
    sys.exit()
