#!/usr/bin/python
from __future__ import print_function
import sys
sys.path.insert(0, './')
sys.path.insert(0, '../')

import numpy as np
import panda_fiber
import config

# =======================================
if __name__ == '__main__':
    import os
    import pylab as plt
    
    filename = "lib/default.ini"
    if len(sys.argv) >= 2:
        filename = sys.argv[1]
        print( "TEST: using config file ",filename )
    
    if os.path.isfile(filename):
        pass
    else:
        filename2 = "../lib/default.ini"
        # try again
        if os.path.isfile(filename2):
            filename = filename2
            pass
        else:
            print( "TEST: !!! can not find file ",filename,filename2 )
            sys.exit()
    
    # set up grid
    rext = 1.1
    Ns = 201
    
    x, y = np.mgrid[ -rext:rext:1j*Ns, -rext:rext:1j*Ns ]
    # print( "TEST: shape ",x.shape )
    
    print( "TEST: calculate stress on disc" )
    sys.stdout.flush()
    param = config.get( filename=filename )
    if param == None:
        print( "TEST: can not get configuration file ",filename )
        sys.exit()
    
    st1 = panda_fiber.Stress(x,y, param=param, debug=True)
    sx = st1.sx
    sy = st1.sy
    # vmax = st1.vmax
    # vmin = st1.vmin
    print( "TEST: disc value range ",st1.vmin," to ",st1.vmax )
    
    # ------------------------------------
    print( "TEST: calculate stress on rim" )
    sys.stdout.flush()
    th = np.arange(-np.pi, np.pi, np.pi/100)
    xrim = np.cos(th)
    yrim = np.sin(th)
    
    st2 = panda_fiber.Stress(xrim,yrim, u=th/np.pi, param=param, blank_outside=False, debug=True)
    
    sr_r = st2.sr
    tt_r = st2.tt
    sr_a = st2.sr_a
    tt_a = st2.tt_a
    
    print( "TEST: rim error s_r ",  np.abs(sr_r).max() )
    print( "TEST: rim error tt_r ", np.abs(tt_r).max() )
    
    # ------------------ plot ------------
    mycmap = plt.cm.coolwarm
    # mycmap = plt.cm.afmhot
    # mycmap = plt.cm.jet
    print( "TEST: now plotting (this will take a while)" )
    sys.stdout.flush()
    
    fig1 = plt.figure("Stress", figsize=(14,6))
    ax1 = fig1.add_subplot(121,aspect='equal')
    ax2 = fig1.add_subplot(122,aspect='equal')

    fig2 = plt.figure("Stress cross-sections", figsize=(14,6))
    ax3 = fig2.add_subplot(111)
    # ax4 = fig2.add_subplot(122)

    # .....................
    st1.plot_xy( ax1=ax1, ax2=ax2)
    st2.plot_bc( ax=ax3 )
    
    print( "TEST: remove colorbar!" )
    ax1.figure.delaxes( st1.cbar_xy1.ax )
    # ax1.figure.delcolorbar()
    
    # print "TEST: try deleting st1"
    # del st1
    
    plt.draw()
    plt.show()
    
    
