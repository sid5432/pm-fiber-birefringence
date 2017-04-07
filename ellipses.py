#!/usr/bin/python
# plot stress/strain tensor as ellipses
import sys
import numpy as np
import pylab as plt
from matplotlib.patches import Ellipse

def plot( sxr, syr, txyr, rext=1.1, ax=None, color='volume',
         subsample=4, debug=False ):
    """
    sxr,syr,txyr: are the sigma_x, sigma_y, tau_xy 2-D stress tensors;
    rext: marks the 2D region that they cover (x=-rext to +rext, y=-rext to rext)
    ax: is a matplotlib/pylab axes; 
    subsample: the number of points to skip along x and y
    color: is coloring scheme to use for the plotted ellipses; 'x' uses the sx 
    stress values, 'y' uses the sy stress values, and 'volume' uses the trace
    sx+sy
    """
    header = "ellipses.plot(): "
    shrink = 0.5
    
    # assume sy and txy are the same shape; no sanity check
    sx  = sxr[::subsample, ::subsample]
    sy  = syr[::subsample, ::subsample]
    txy = txyr[::subsample, ::subsample]
    
    Nx, Ny = sx.shape
    x, y = np.mgrid[ -rext:rext:1j*Nx, -rext:rext:1j*Ny ]
    
    if ax != None:
        ax.cla()
        ax.set_axis_bgcolor('#222222')
        # ax.figure.patch.set_facecolor('blue') # set surrounding
        
        # gradient
        # plotlim = plt.xlim() + plt.ylim()
        # ax.imshow([[-rext,-rext],[rext,rext]], cmap=plt.cm.Greens, interpolation='bicubic', extent=plotlim)
        
        
    
    # NOTE: ex + ey is the thermoelastic strain = volumetric strain (assuming ez=0)
    # in the absence of all other stress, ex = sx/E, where E is Young's modulus.
    # s_tr (trace) is an invariant (of coordinate rotation)
    s_tr = sx + sy
    zmap = s_tr
    label = r'$\sigma_{trace}$'
    if color=='volume':
        pass
    elif color=='x':
        zmap = sx
        label = r'$\sigma_x$'
    elif color=='y':
        zmap = sy
        label = r'$\sigma_y$'
    else:
        print header," unrecognized color mapping ",color,"; using sx+sy"
    
    if ax != None:
        ax.set_title("Stress distribution; color represents "+label)
    
    Nx, Ny = sx.shape
    x, y = np.mgrid[ -rext:rext:1j*Nx, -rext:rext:1j*Ny ]
    
    # get eigenvalues
    dxy = sx - sy
    axy = sx + sy
    c2  = 2.0*txy
    D = np.sqrt( dxy**2 + c2**2 )
    ev1 = (axy + D)/2.0
    ev2 = (axy - D)/2.0 # ev1 > ev2
    
    th = np.arctan2( c2 , dxy )/2.0
    # NOTE: the angles thus calculated might be off by PI
    # (depending on how we choose the direction of the eigenvector)
    # but this is not important as far as depicting the ellipses
    # is concerned.
    #
    # print "DEBUG: th range ",th.min()," to ",th.max()
    
    # for Ellipse we need to use degrees!
    ang = (180.0/np.pi) * th
    
    # print "DEBUG: th range ",th.min(), th.max()
    # print "DEBUG: ev1 range ",ev1.min(), ev1.max()
    # print "DEBUG: ev2 range ",ev2.min(), ev2.max()
    emax = np.abs(ev1).max()
    tmp  = np.abs(ev2).max()
    if emax < tmp:
        emax = tmp
    
    # baseline = np.abs(ev2).min()*2.0
    baseline = emax*1.5
    
    # scaling factor for size of ellipses
    efactor = 1.0*rext/Nx/emax
    
    # this is added to ev1 and ev2 so that the ellipses have positive major/minor axes
     
    # -------------------------------------------------------
    # iterate through all grid points inside unit circle
    r2 = x*x + y*y
    w = np.where(r2<1)
    
    # for fill-color mapping
    vmin = zmap[w].min()
    vmax = zmap[w].max()
    norm = vmax - vmin
    
    # ax.pcolor(x,y, s_tr, cmap=plt.cm.coolwarm)
    
    if ax == None:
        # nothing else to do
        return
    
    ax.set_xlim( -rext, rext )
    ax.set_ylim( -rext, rext )
    
    for i in range(len(w[0])):
        # print "point-0 ",w[0]
        # print "point-1 ",w[1]
        ix = w[0][i]
        iy = w[1][i]
        xp = x[ix,iy]
        yp = y[ix,iy]
        zp = zmap[ix,iy]
        # print "coord x=", xp," y=",yp," z=",zp
        angle = ang[ix,iy]
        
        # NOTE: concerning the sign of stress, etc., see notes
        # in help.py
        A = (baseline + ev1[ix,iy] ) * efactor
        B = (baseline + ev2[ix,iy] ) * efactor
        
        el = Ellipse( xy=[xp,yp], width=A, height=B, angle=angle, linewidth=0)
        ax.add_artist(el)
        # el.set_clip_box( ax.bbox )
        el.set_alpha( 0.8 )
        el.set_facecolor( plt.cm.coolwarm( (zp-vmin)/norm ) )
    
    ax.figure.canvas.draw()
    
    return

# ==============================================================
if __name__ == '__main__':
    sx,sy,txy = np.load("tests/snapshot.npy")
    print "TEST: recover shape ", sx.shape, sy.shape, txy.shape
    
    fig = plt.figure("TEST", figsize=(15,8))
    # ax1 = fig.add_subplot(111, aspect='equal', axisbg='black')
    ax1 = fig.add_subplot(111, aspect='equal')
    
    plot( sx,sy,txy, ax=ax1, debug=True, subsample=4, color='x' )
    
    plt.show()
    
