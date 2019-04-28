#!/usr/bin/python
# Reference: see info in help.py
# this module implements the stress calculation of panda fibers
from __future__ import print_function
import numpy as np
import sys
import pylab as plt

class Stress:
    def __init__(self, x, y, u=None, param=None, blank_outside=True, 
                 save=False, debug=False):
        self.debug = debug
        self.blank_outside = blank_outside
        self.u = u # implicit parameter that generates x and y (used for boundary check)
        
        # we assume x,y are numpy arrays
        self.x = x
        self.y = y
        self.param = param
        
        a = param['a'] # core radius
        self.a2 = a*a
        
        self.R2 = x**2 + y**2
        self.R  = np.sqrt(self.R2)
        self.th = np.arctan2(y, x)
        
        self.cs = np.cos(self.th)
        self.sn = np.sin(self.th)
        self.cs2 = self.cs**2
        self.sn2 = self.sn**2
        self.sncs = self.cs*self.sn
        
        # initialize
        self.ax_xy1 = None
        self.ax_xy2 = None
        self.ax_xy3 = None
        self.ax_rt1 = None
        self.ax_rt2 = None
        
        self.cbar_xy1 = None
        self.cbar_xy2 = None
        self.cbar_xy3 = None
        self.cbar_rt1 = None
        self.cbar_rt2 = None
        
        self.w_outside = np.where( self.R2 > 1.0 )
        self.w_clad = np.where( self.R2 <= 1.0 )
        self.w_core = np.where( self.R2 > self.a2 )
        
        # calculate stress distribution        
        self.total(save=save)
        
        return
    
    # ----------------------------------------------------
    # we use this to calculate the stress potential contributions from
    # the core and the two circular SAP regions
    def sap(self, xc, yc, a, beta ):
        # (xc,yc) is the center of the circular SAP region
        # a is the normalized radius (wrt cladding radius)
        x = self.x
        y = self.y
        Dx = (x - xc)
        Dy = (y - yc)
        rr2 = Dx**2 + Dy**2
        th2 = 2.0*np.arctan2(Dy, Dx) # angle wrt (xc,yc)!
        a2 = a**2
        
        w = np.where( rr2 > a2 ) # outside core
        
        # inside core values:
        sx  = -np.ones_like(x) # sigma_r; but r and theta are with respect to (xc,yc), not (0,0)!
        sy  = -np.ones_like(x)
        txy = np.zeros_like(x)
        
        # ...outside
        sx[w]  = -a2/rr2[w]*np.cos(th2[w])
        sy[w]  = -sx[w] # +a2/rr2[w]*np.cos(th2[w])
        txy[w] = -a2/rr2[w]*np.sin(th2[w])
        
        param = self.param
        # multiply by E/(1+nu)*beta*DT/2.0;
        # beta is beta1 or beta3, depending on region
        E  = param['E']    # Young's modulus
        nu = param['nu']   # Poisson's ratio
        DT = -param['DT']   # difference between ambient and glass melting/drawing temperature
        C  = param['C']    # photoelastic coefficient
        
        # NOTE: we are calculating stress, but scaled to delta n; however
        # delta n is of opposite sign than stress!
        Q  = C*E/(1.0 + nu)*beta*DT/2.0 * 1e4
        sx  = Q*sx
        sy  = Q*sy
        txy = Q*txy
        
        # now calculate sr (sigma_r) and tt (tau_r,theta) wrt to (0,0)
        th = self.th # np.arctan2(y, x)
        cs = self.cs # np.cos(th)
        sn = self.sn # np.sin(th)
        cs2 = self.cs2 # cs**2
        sn2 = self.sn2 # sn**2
        crs = self.sncs # cs*sn
        
        sr = sx*cs2 + sy*sn2 + 2.0*txy*crs  # sigma_r
        tt = (sy-sx)*np.sin( 2.0*th )/2.0 + txy*np.cos( 2.0*th )
        
        return sx, sy, txy, sr, tt
    
    # ----------------------------------------------------
    # we use this to calculate the Airy stress function contribution
    def airy(self):
        # x,y are numpy arrays;
        # N is the order at which we truncate the sum
        # all dimensions normalized wrt cladding radius
        
        x = self.x
        y = self.y
        param = self.param
        N = param['N']
        
        beta1 = param['beta1'] # beta1 - beta2
        beta3 = param['beta3'] # beta3 - beta2
        a     = param['a']     # core radius
        d1    = param['d1']    # SAP radius
        d2    = param['d2']    # distance from SAP center to fiber center
        E     = param['E']     # Young's modulus
        nu    = param['nu']    # Poisson's ratio
        C     = param['C']     # Elasto-optic coefficient
        DT    = -param['DT']   # difference between ambient and glass melting/drawing temperature
        
        R2 = self.R2 # x**2 + y**2
        R  = self.R  # np.sqrt(R2)
        th = self.th # np.arctan2(y, x)
        
        # A = 1.0 + 0.5*beta1/beta3*(a/d1)**2
        
        # sr = A*np.ones_like(x)  # sigma_r
        # st = A*np.ones_like(x)  # sigma_theta
        sr = np.ones_like(x)  # sigma_r
        st = np.ones_like(x)  # sigma_theta
        tt = np.zeros_like(x)  # tau_r,theta
        
        d22 = d2*d2
        th2 = 2.0*th
        theta_t = th2
        mult = d22 * R2
        
        for n in range(1,N+1):
            m1 = (4*n*n - 1)*d22
            m2 = 2*(2*n+1)*(n-1)
            m3 = 2*(2*n+1)*(n+1)
            m4 = 2*(2*n+1)*n
            
            theta_t = 2.0*n*th
            cs = np.cos(theta_t)
            sn = np.sin(theta_t)
            
            xp1 = mult**(n-1)
            xp2 = mult**n
            
            sr = sr + (m1*xp1 - m2*xp2)*cs
            st = st - (m1*xp1 - m3*xp2)*cs
            tt = tt - (m1*xp1 - m4*xp2)*sn
            
            # next round
            # theta_t = theta_t + th2
            # theta_t = 2.0*n*th
            # xp1 = xp2
            # xp2 = xp2*mult

        # NOTE: we are calculating stress, but scaled to delta n; however
        # delta n is of opposite sign than stress!
        Q = C*E*DT/(1+nu) * 1e4
        ff = beta3*d1*d1
        gg = 0.5*beta1*a*a
        
        sr = Q*( ff*sr + gg )
        st = Q*( ff*st + gg )
        tt = Q*( ff*tt )
        
        # convert to x y components
        cs = self.cs # np.cos(th)
        sn = self.sn # np.sin(th)
        cs2 = self.cs2 # cs**2
        sn2 = self.sn2 # sn**2
        # cr2 = 2.0*self.sncs # 2.0*cs*sn
        
        # components in cylindrical coordinates
        self.sr_a = sr # $\sigma_r$
        self.tt_a = tt # $\tau_{r,\theta}$
        self.st_a = st # $\sigma_\theta$
        
        # convert to components in rectangular coordinates
        tmp = tt*2.0*self.sncs
        
        self.sx_a = sr*cs2 + st*sn2 - tmp
        self.sy_a = sr*sn2 + st*cs2 + tmp
        self.txy_a = (sr-st)*self.sncs + tt*(cs2 - sn2)
        
        return
    
    # ----------------------------------------------------
    # all all stress contributions
    def total(self, save=False):
        # total stress: contributions from the core, the SAP's, and the Airy function
        header = "stress.total(): "
        
        x = self.x
        y = self.y
        param = self.param
        N = param['N']
        blank_outside=self.blank_outside
        
        print( header," calculate core contribution" ) # ------------------------------------
        sys.stdout.flush()
        xc = 0.0
        yc = 0.0
        # param['beta'] = param['beta1']
        sx_core, sy_core, txy_core, sr_core, tt_core = self.sap(xc,yc,
                                                                param['a'], param['beta1'])
        
        print( header," calculate SAP#1 contribution" ) # ------------------------------------
        sys.stdout.flush()
        
        # print( header," SAP center ",param['d2']," radius ",param['d1'] )
        xc = param['d2']
        yc = 0.0
        sx_sap1, sy_sap1, txy_sap1, sr_sap1, tt_sap1 = self.sap(xc,yc, param['d1'], param['beta3'])
        
        print( header," calculate SAP#2 contribution" ) # ------------------------------------
        sys.stdout.flush()
        
        xc = -param['d2']
        yc = 0.0
        sx_sap2, sy_sap2, txy_sap2, sr_sap2, tt_sap2 = self.sap(xc,yc, param['d1'], param['beta3'])
        
        # ------------------------------------
        print( header," calculate Airy stress function (background), N=", N )
        sys.stdout.flush()
        self.airy()
        # results: self.sx_a, self.sy_a, self.sr_a, self.tt_a, self.st_a
        
        print( header," calculate total" ) # ------------------------------------
        sys.stdout.flush()
        sx  = sx_core + sx_sap1 + sx_sap2 + self.sx_a
        sy  = sy_core + sy_sap1 + sy_sap2 + self.sy_a
        txy = txy_core + txy_sap1 + txy_sap2 + self.txy_a
        
        self.sr = sr_core + sr_sap1 + sr_sap2 + self.sr_a
        self.tt = tt_core + tt_sap1 + tt_sap2 + self.tt_a
        
        # blank outside
        if blank_outside:
            w = np.where( self.R2 > 1.0 )
            sx[w] = 0.0
            sy[w] = 0.0
            self.sr[w] = 0.0
            self.tt[w] = 0.0
        
        # get min and max for x-y
        vmin = sx.min()
        tmp  = sy.min()
        if vmin > tmp:
            vmin = tmp
        
        vmax = sx.max()
        tmp  = sy.max()
        if vmax < tmp:
            vmax = tmp
        
        vmin2 = self.sr.min()
        vmin3 = self.tt.min()
        
        vmax2 = self.sr.max()
        vmax3 = self.tt.max()
        
        # ... again...
        if blank_outside:
            sx[w] = vmin
            sy[w] = vmin
            self.sr[w] = vmin2
            self.tt[w] = vmin3
        
        self.vmin = vmin
        self.vmax = vmax
        
        self.vmin_r = vmin2
        self.vmax_r = vmax2
        
        self.vmin_t = vmin3
        self.vmax_t = vmax3

        print( header," x-y val range ",vmin," to ",vmax )
        print( header," r val range ",vmin2," to ",vmax2 )
        print( header," t val range ",vmin3," to ",vmax3 )

        self.sx = sx
        self.sy = sy
        self.txy = txy
        
        if save:
            print( header," saving snapshot...." )
            np.save("snapshot", (self.sx, self.sy, self.txy))
        
        #  sx, sy, sr, tt, sr_a, tt_a, vmin, vmax
        return
    
    #---------------------------------------------------------
    def plot_bc(self, ax=None):
        if ax == None:
            return
        
        ax.cla()
        ax.set_xlabel(r'$\theta/\pi$')
        ax.set_ylabel(r'Birefringence ($\times 10^{-4}$)')
        
        ax.set_title("Birefringence along cladding boundary (with N=%d Airy stress function terms)" % self.param['N'])
        ax.plot(self.u, self.sr, label=r'$\sigma_r$')
        ax.plot(self.u, self.tt, label=r'$\sigma_{r\theta}$')
        ax.legend( loc='best' )
        
        ax.plot(self.u, self.sr_a, label=r'$\sigma_r$ (Airy)')
        ax.plot(self.u, self.tt_a, label=r'$\sigma_{r\theta}$ (Airy)')
        ax.legend( loc='best' )
        ax.figure.canvas.draw()
        
        return
    
    #---------------------------------------------------------
    def clear_xy_colorbars(self):
        if self.ax_xy1 != None:
            self.ax_xy1.figure.delaxes( self.cbar_xy1.ax )
            self.ax_xy1.figure.subplots_adjust(right=0.90)
            
        if self.ax_xy2 != None:
            self.ax_xy2.figure.delaxes( self.cbar_xy2.ax )
            self.ax_xy2.figure.subplots_adjust(right=0.90)
            
        if self.ax_xy3 != None:
            self.ax_xy3.figure.delaxes( self.cbar_xy3.ax )
            self.ax_xy3.figure.subplots_adjust(right=0.90)
        
        return
            
    def clear_rt_colorbars(self):
        if self.ax_rt1 != None:            
            self.ax_rt1.figure.delaxes( self.cbar_rt1.ax )
            self.ax_rt1.figure.subplots_adjust(right=0.90)
        
        if self.ax_rt2 != None:            
            self.ax_rt2.figure.delaxes( self.cbar_rt2.ax )
            self.ax_rt2.figure.subplots_adjust(right=0.90)
        
        return
    
    #---------------------------------------------------------
    def plot_xy(self, ax1=None, ax2=None, ax3=None, plotdn=True):
        rext = 1.1 # fixed
        mycmap = plt.cm.coolwarm
        # mycmap = plt.cm.afmhot
        # mycmap = plt.cm.jet
        
        param = self.param
        
        if plotdn:
            # parameters for calculating dn;
            # we divide by C because all stresses are already normalized by C
            AA = (param['C1'] + param['nu'] * param['C2']) / param['C']
            BB = (param['C2'] * ( 1.0 + param['nu'] )) / param['C']
            XSURF = np.zeros_like(self.x)
            YSURF = np.zeros_like(self.x)
            
            XSURF[self.w_clad] = -AA*self.sx[self.w_clad] - BB*self.sy[self.w_clad]
            YSURF[self.w_clad] = -BB*self.sx[self.w_clad] - AA*self.sy[self.w_clad]
            vmin = XSURF.min()
            tmp = YSURF.min()
            if tmp < vmin:
                vmin = tmp
            vmax = XSURF.max()
            tmp = YSURF.max()
            if tmp > vmax:
                vmax = tmp
            
            XSURF[self.w_outside] = vmin
            YSURF[self.w_outside] = vmin
            
            title_x = r'$\Delta n_x \times 10^{-4}$'
            title_y = r'$\Delta n_y \times 10^{-4}$'
        else:
            XSURF = self.sx
            YSURF = self.sy
            vmin = self.vmin
            vmax = self.vmax
            title_x = r'$\sigma_x$ stress (scaled to $\Delta n \times 10^{-4}$)'
            title_y = r'$\sigma_y$ stress (scaled to $\Delta n \times 10^{-4}$)'
        
        TSURF = np.zeros_like(self.x)
        TSURF[self.w_clad] = self.txy[self.w_clad]
        tmin = TSURF.min()
        tmax = TSURF.max()
        TSURF[self.w_outside] = tmin
        
        # ............ now plot ............................
        if ax1 != None:
            self.ax_xy1 = ax1
            ax1.cla()
            
            ax1.set_title(title_x)
            ax1.set_xlim(-rext,rext)
            ax1.set_ylim(-rext,rext)
            cax1 = ax1.pcolor(self.x, self.y, XSURF, cmap=mycmap, vmin=vmin, vmax=vmax )
            self.cbar_xy1 = ax1.figure.colorbar( cax1, ax=ax1, shrink=0.6, use_gridspec=False )
            ax1.figure.canvas.draw()
        
        if ax2 != None:
            self.ax_xy2 = ax2
            ax2.cla()
                
            ax2.set_title(title_y)
            ax2.set_xlim(-rext,rext)
            ax2.set_ylim(-rext,rext)
            cax2 = ax2.pcolor(self.x, self.y, YSURF, cmap=mycmap, vmin=vmin, vmax=vmax )
            self.cbar_xy2 = ax2.figure.colorbar( cax2, ax=ax2, shrink=0.6, use_gridspec=False )
            ax2.figure.canvas.draw()
            
        if ax3 != None:
            self.ax_xy3 = ax3
            ax3.cla()
                
            ax3.set_title(r'$\sigma_{xy}$ stress (scaled to $\Delta n \times 10^{-4}$)')
            ax3.set_xlim(-rext,rext)
            ax3.set_ylim(-rext,rext)
            cax3 = ax3.pcolor(self.x, self.y, TSURF, cmap=mycmap, vmin=tmin, vmax=tmax )
            self.cbar_xy3 = ax3.figure.colorbar( cax3, ax=ax3, shrink=0.6, use_gridspec=False )
            ax3.figure.canvas.draw()
        
        
        return
    
    #---------------------------------------------------------
    def plot_rt(self, ax1=None, ax2=None):
        rext = 1.1 # fixed
        mycmap = plt.cm.coolwarm
        # mycmap = plt.cm.afmhot
        # mycmap = plt.cm.jet

        if ax1 != None:
            self.ax_rt1 = ax1
            ax1.cla()
            
            ax1.set_title(r'$\sigma_r$ stress (scaled to $\Delta n \times 10^{-4}$)')
            ax1.set_xlim(-rext,rext)
            ax1.set_ylim(-rext,rext)
            cax1 = ax1.pcolor(self.x, self.y, self.sr, cmap=mycmap, vmin=self.vmin_r, vmax=self.vmax_r )
            self.cbar_rt1 = ax1.figure.colorbar( cax1, ax=ax1, shrink=0.8, use_gridspec=False )
            ax1.figure.canvas.draw()
        
        if ax2 != None:
            self.ax_rt2 = ax2
            ax2.cla()
                
            ax2.set_title(r'$\sigma_{r\theta}$ stress (scaled to $\Delta n \times 10^{-4}$)')
            ax2.set_xlim(-rext,rext)
            ax2.set_ylim(-rext,rext)
            cax2 = ax2.pcolor(self.x, self.y, self.tt, cmap=mycmap, vmin=self.vmin_t, vmax=self.vmax_t )
            self.cbar_rt2 = ax2.figure.colorbar( cax2, ax=ax2, shrink=0.8, use_gridspec=False )
            ax2.figure.canvas.draw()
        
        return
    

# =======================================
if __name__ == '__main__':
    print( "Test with test_panda.py" )
