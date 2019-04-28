#!/usr/bin/python
# test tensor display
from __future__ import print_function
import sys
import numpy as np
import pylab as plt
from matplotlib.patches import Ellipse

# ==============================================================
if __name__ == '__main__':
    sx,sy,txy = np.load("snapshot.npy")
    print( "recover shape ", sx.shape, sy.shape, txy.shape )
    
    # generate grid
    rext = 1.1
    Nx,Ny = sx.shape
    shrink = 0.5
    x, y = np.mgrid[ -rext:rext:1j*Nx, -rext:rext:1j*Ny ]
    
    xnew = x[::2,::2]
    ynew = y[::2,::2]
    print( "shape ", xnew.shape, ynew.shape)
    
    fig = plt.figure("TEST", figsize=(15,6))
    ax1 = fig.add_subplot(131, aspect='equal')
    ax2 = fig.add_subplot(132, aspect='equal')
    ax3 = fig.add_subplot(133, aspect='equal')
    
    s_trace = sx + sy
    
    vmin = s_trace.min()
    vmax = s_trace.max()
    print( "TEST: range ",vmin," to ",vmax )
    print( "TEST: now plotting ....." )
    sys.stdout.flush()
    
    # ------------------------------------------------------
    cax1 = ax1.pcolor(x[::4,::4],y[::4,::4], sx[::4,::4], cmap=plt.cm.coolwarm)
    cb1  = ax1.figure.colorbar( cax1, ax=ax1, shrink=shrink )
    ax1.set_xlim(-rext,rext)
    ax1.set_ylim(-rext,rext)
    ax1.set_title('sx (subsample-4)')
    
    cax2 = ax2.pcolor(xnew,ynew, sy[::2,::2], cmap=plt.cm.coolwarm)
    cb2  = ax2.figure.colorbar( cax2, ax=ax2, shrink=shrink )
    ax2.set_xlim(-rext,rext)
    ax2.set_ylim(-rext,rext)
    ax2.set_title('sy (subsample-2)')
    
    cax3 = ax3.pcolor(xnew,ynew, s_trace[::2,::2], cmap=plt.cm.coolwarm)
    cb3  = ax3.figure.colorbar( cax3, ax=ax3, shrink=shrink )
    ax3.set_xlim(-rext,rext)
    ax3.set_ylim(-rext,rext)
    ax3.set_title('Total')

    # cax = ax.pcolor(x,y, s_trace, cmap=plt.cm.coolwarm, vmin=1, vmax=2.0)
    # ------------------------------------------------------
    print( "TEST: add three ellipses on top of 'Total'" )
    sys.stdout.flush()
    el = Ellipse(xy=[0.0,0.0], width=1.4, height=0.2, angle=30.0)
    ax2.add_artist(el)
    el.set_clip_box( ax2.bbox )
    el.set_alpha( 0.5 )
    el.set_facecolor( np.array([0.1, 0.0, 0.9]) )
    
    el = Ellipse(xy=[0.1,-0.2], width=0.15, height=0.08, angle=-30.0)
    ax2.add_artist(el)
    el.set_clip_box( ax2.bbox )
    # el.set_alpha( 0.5 )
    # el.set_facecolor( np.array([0.1, 0.0, 0.9]) )
    
    # el = Ellipse(xy=[0.1,0.0], width=0.15, height=0.08, angle=90.0, fill=False)
    el = Ellipse(xy=[0.1,0.0], width=0.15, height=0.08, angle=90.0)
    ax2.add_artist(el)
    el.set_clip_box( ax2.bbox )
    # el.set_alpha( 0.5 )
    el.set_facecolor( plt.cm.coolwarm(0.1) )
    
    plt.draw()
    plt.show()
    
    
