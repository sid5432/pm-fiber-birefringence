#!/usr/bin/python
from __future__ import print_function
import sys
sys.path.insert(0, './')
sys.path.insert(0, '../')

import numpy as np
import bowtie_fiber
import config
import os
import pylab as plt

# =======================================
def test_bowtie():
    
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
    
    st1 = bowtie_fiber.Stress(x,y, param=param, debug=True)
    
    print( "TEST: disc value range ",st1.vmin," to ",st1.vmax )
    
    # histogram
    # hist = np.histogram( st1.sx, bins=20 )
    # print( "TEST: histogram ", hist )
    # np.savetxt("hist.dat", hist, delimiter='\t')
    
    # revise
    print( "TEST: force revision of range to +/-5" )
    st1.vmin = -5
    st1.vmax = 5
    
    # ------------------------------------
    print( "TEST: calculate stress on rim")
    sys.stdout.flush()
    th = np.arange(-np.pi, np.pi, np.pi/100)
    xrim = np.cos(th)
    yrim = np.sin(th)
    
    st2 = bowtie_fiber.Stress(xrim,yrim, u=th/np.pi, param=param, blank_outside=False, debug=True)
    
    sr_r = st2.sr
    tt_r = st2.tt
    sr_a = st2.sr_a
    tt_a = st2.tt_a
    
    print( "TEST: rim error s_r ",  np.abs(sr_r).max() )
    print( "TEST: rim error tt_r ", np.abs(tt_r).max() )
    assert( np.abs(sr_r).max() < 1e-15 )
    assert( np.abs(tt_r).max() < 1e-15 )
    
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
    # ax3 = fig2.add_subplot(121)
    # ax4 = fig2.add_subplot(122)
    
    # .....................
    st2.plot_bc( ax=ax3 )
    # print( "TEST: force range to +/-5" )
    # st1.plot_rt( ax1=ax1, ax2=ax2, vmin=-5, vmax=5)
    mu = np.mean( st1.sx[st1.w_clad] )
    sigma = np.std( st1.sx[st1.w_clad] )
    print( "TEST: stats: %f +/- %f " % (mu, sigma) )
    lb = mu - 3*sigma
    ub = mu + 3*sigma
    print( "TEST: stats: %f +/- %f " % (mu, sigma) )
    sys.stdout.flush()
    
    st1.plot_xy( ax1=ax1, ax2=ax2)
    
    # ax4.hist( st1.sx[st1.w_clad] )
    
    # print( "TEST: remove colorbar!")
    # ax1.figure.delaxes( st1.cbar_xy1.ax )
    ## ax1.figure.delcolorbar()
    
    # print( "TEST: try deleting st1" )
    # del st1
    
    plt.draw()
    plt.show()
    
    
# =======================================
if __name__ == '__main__':
    test_bowtie()
