#!/usr/bin/python
# Reference: see info in help.py;
# this module implements the stress calculation of bow-tie fibers
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
        self.ax_xy1 = None # plot for sigma_x
        self.ax_xy2 = None # plot for sigma_y
        self.ax_xy3 = None # plot for sigma_xy
        self.ax_rt1 = None # plot for sigma_r
        self.ax_rt2 = None # plot for sigma_theta
        
        self.cbar_xy1 = None # colorbars for each plot above
        self.cbar_xy2 = None
        self.cbar_xy3 = None
        self.cbar_rt1 = None
        self.cbar_rt2 = None
        
        # regions:
        self.r1 = param['r1']     # bow-tie radius r1
        self.r2 = param['r2']     # bow-tie radius r2
        btr1 = self.r1**2
        btr2 = self.r2**2
        
        self.w_outside = np.where( self.R2 > 1.0 )
        self.w_clad = np.where( self.R2 <= 1.0 )
        self.w_core = np.where( self.R2 > self.a2 )
        self.w_lt_r1 = np.where( self.R2 < btr1 )
        # give it a little margin...
        self.w_gt_r2 = np.where( np.logical_and(self.R2 > btr2, self.R2<1.05) )
        self.w_between = np.where( np.logical_and(self.R2 >=btr1, self.R2<=btr2) )
        
        # calculate stress distribution
        self.total(save=save)
        
        return
    
    # ----------------------------------------------------
    # we use this to calculate the stress potential contributions from
    # the core
    def core(self):
        # (xc,yc) is the center of the circular SAP region
        # a is the normalized radius (wrt cladding radius)
        xc, yc = 0, 0
        a = self.param['a']
        beta = self.param['beta1']
        
        x = self.x
        y = self.y
        Dx = (x - xc)
        Dy = (y - yc)
        rr2 = Dx**2 + Dy**2
        th2 = 2.0*np.arctan2(Dy, Dx)
        a2 = a**2
        
        cs2 = np.cos(th2)
        w = np.where( rr2 > a2 ) # outside core
        sx = -np.ones_like(x)
        sy = -np.ones_like(x)
        
        sr = -np.ones_like(x) # sigma_r; but r and theta are with respect to (xc,yc), not (0,0)!
        sr[w] = -a2/rr2[w]
        
        # sx[w] = sr[w] * cs2[w]  # sx[w] = -a2 * cs2[w]/rr2[w]
        # sy[w] = -sx[w]          # sy[w] = a2 * cs2[w]/rr2[w]
        
        sx[w] = -a2 * cs2[w]/rr2[w]
        sy[w] = a2 * cs2[w]/rr2[w]
        
        # we also need this to verify the boundary condition is good
        txy = np.zeros_like(x)
        txy[w] = sr[w]*np.sin(th2[w])
        
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
        
        sx = Q*sx
        sy = Q*sy
        txy = Q*txy # extra (see above)
        
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
    
    # ----------------------------------------------------------
    # this calculates the stress contribution from the bow-tie SAP
    # regions
    def sap(self):
        header = "bowtie.sap(): "
        
        sr = np.zeros_like(self.R)
        st = np.zeros_like(self.R)
        tt = np.zeros_like(self.R)
        
        x = self.x
        y = self.y
        param = self.param
        N = param['N']
        
        beta3 = param['beta3'] # beta3 - beta2
        
        r1     = param['r1']     # bow-tie radius r1
        r2     = param['r2']     # bow-tie radius r2
        theta1 = param['theta1'] # bow-tie half-angle
        
        E     = param['E']     # Young's modulus
        nu    = param['nu']    # Poisson's ratio
        C     = param['C']     # Elasto-optic coefficient
        DT    = -param['DT']   # difference between ambient and glass melting/drawing temperature
        
        R2 = self.R2 # x**2 + y**2
        R  = self.R  # np.sqrt(R2)
        th = self.th # np.arctan2(y, x)
        
        r1sq = r1**2
        r2sq = r2**2
        
        th2 = 2.0*th
        theta_t = th2.copy()
        # NOTE: we keep th2 fixed, and only increment theta_t as we go through summation
        pw1 = 1.0
        pw2 = 1.0
        xth = theta1t = theta1*2
        
        fact_cos = np.zeros_like(x)
        fact_sin = np.zeros_like(x)
        xp1      = np.zeros_like(x)
        xp2      = np.zeros_like(x)
        
        # round one: r < r1 ...................................................
        w = self.w_lt_r1
        # for n == 1:
        tmp = -(np.log(r2) - np.log(r1)) * np.sin(xth)
        sr[w] = sr[w] + tmp * np.cos(th2[w])
        tt[w] = tt[w] - tmp * np.sin(th2[w])
        
        for n in range(2,N+1):
            # prep
            theta_t = theta_t + th2 # $2 n \theta$
            xth = xth + theta1t     # $2 n \theta_1$
            pw1 /= r1sq # $1/r_1^{2(n-1)}$
            pw2 /= r2sq # $1/r_2^{2(n-1)}$
            
            # calc
            mm = 0.5*(2*n-1)/float( n*(n-1) )
            common = mm*(pw2-pw1)*np.sin( xth )     # $(r_2^{2(n+1)} - r_1^{2(n+1)}) \sin( 2n\theta_1 )$
            fact_cos[w] = common* np.cos( theta_t[w] ) # $\cos(2 n \theta)$
            fact_sin[w] = common* np.sin( theta_t[w] ) # $\sin(2 n \theta)$
            
            xp1[w] = R2[w]**(n-1)
            
            sr[w] = sr[w] + fact_cos[w] *xp1[w]
            tt[w] = tt[w] - fact_sin[w] *xp1[w]
        
        st[w] = -sr[w] # note: r<r1 is outside the "source"
        
        # round two: r > r2 ...................................................
        w = self.w_gt_r2
        tmp = np.ones_like(x)
        tmp[w] = theta1*(r2sq-r1sq)/R2[w]
        sr[w] = sr[w] - tmp[w]
        st[w] = st[w] + tmp[w]
        
        theta_t = th2.copy()
        pw1 = r1sq
        pw2 = r2sq
        xth = theta1t = theta1*2
        
        for n in range(1,N+1):
            pw1 *= r1sq # $r_1^{2(n+1)}$
            pw2 *= r2sq # $r_2^{2(n+1)}$
            
            mm = 0.5*(2*n+1)/float( n*(n+1) )
            common = mm*(pw2-pw1)*np.sin( xth )     # $(r_2^{2(n+1)} - r_1^{2(n+1)}) \sin( 2n\theta_1 )$
            fact_cos[w] = common* np.cos( theta_t[w] ) # $\cos(2 n \theta)$
            fact_sin[w] = common* np.sin( theta_t[w] ) # $\sin(2 n \theta)$
            
            xp1[w] = 1.0/R2[w]**(n+1)
            
            sr[w] = sr[w] - fact_cos[w] *xp1[w]
            tt[w] = tt[w] - fact_sin[w] *xp1[w]
            
            # next round
            theta_t = theta_t + th2 # $2 n \theta$
            xth = xth + theta1t     # $2 n \theta_1$
        
        st[w] = -sr[w] # note r>r2 is outside the "source"
        
        # round thre: r1 < r < r2 ...................................................
        w = self.w_between
        sr[w] = sr[w] + theta1*( r1sq/R2[w] -1.0 )
        st[w] = st[w] - theta1*( r1sq/R2[w] +1.0 )
        
        theta_t = th2.copy()
        pw1 = r1sq
        pw2 = 1.0/r2sq
        xth = theta1t = theta1*2

        # for n == 1:
        sr[w] = sr[w] - ( np.log(r2)-np.log(R[w]) + 0.75*(1-r1**4/R2[w]**2) ) * np.sin(xth) * np.cos( th2[w] )
        st[w] = st[w] + ( np.log(r2)-np.log(R[w]) - 1.25 - 0.75 * r1**4/R2[w]**2 ) * np.sin(xth) * np.cos( th2[w] )
        tt[w] = tt[w] + ( np.log(r2)-np.log(R[w]) - 0.75*(1-r1**4/R2[w]**2) ) * np.sin(xth) * np.sin( th2[w] )
        
        for n in range(2,N+1):
            # prep
            theta_t = theta_t + th2 # $2 n \theta$
            xth = xth + theta1t     # $2 n \theta_1$
            
            pw1 *= r1sq # $r_1^{2(n+1)}$
            pw2 *= r2sq # $1/r_2^{2(n+1)}$
            
            xp1[w] = 1.0/R2[w]**(n+1)
            xp2[w] = R2[w]**(n-1)
            
            common = np.sin( xth ) # $\sin( 2n\theta_1 )$
            fact_cos[w] = common* np.cos( theta_t[w] ) # $\sin( 2n\theta_1 ) \cos(2 n \theta)$
            fact_sin[w] = common* np.sin( theta_t[w] ) # $\sin( 2n\theta_1 ) \sin(2 n \theta)$
            
            m11 = float( 2*n**2-1 )/float( n*(n**2 - 1) )
            m12 = -float(2*n+1)/float( 2*n*(n+1) )
            m13 = -float(2*n-1)/float( 2*n*(n-1) )
            
            m21 = 1.0/float(n*n**2-1)
            m22 = m12 # -float(2*n+1)/float( 2*n*(n+1) )
            m23 = m13 # -float(2*n-1)/float( 2*n*(n-1) )
            
            m31 = m21  # 1.0/float(n*n**2-1)
            m32 = -m12 # float(2*n+1)/float( 2*n*(n+1) )
            m33 = m13  # -float(2*n-1)/float( 2*n*(n-1) )
            
            sr[w] = sr[w] - (m11 + m12*pw1*xp1[w] + m13*pw2*xp2[w])* fact_cos[w]
            st[w] = st[w] + (m21 + m22*pw1*xp1[w] + m23*pw2*xp2[w])* fact_cos[w]
            tt[w] = tt[w] + (m31 + m32*pw1*xp1[w] + m33*pw2*xp2[w])* fact_sin[w]
            
        
        # NOTE: we are calculating stress, but scaled to delta n; however
        # delta n is of opposite sign than stress!
        Q  = C*E*DT/(1.0 + nu) * beta3/np.pi * 1e4
        sr = sr * Q
        st = st * Q
        tt = tt * Q
        
        # convert to x y components
        cs = self.cs # np.cos(th)
        sn = self.sn # np.sin(th)
        cs2 = self.cs2 # cs**2
        sn2 = self.sn2 # sn**2
        
        tmp = tt*2.0*self.sncs
        
        sx  = sr*cs2 + st*sn2 - tmp
        sy  = sr*sn2 + st*cs2 + tmp
        txy = (sr-st)*self.sncs + tt*(cs2 - sn2)
        
        return sx, sy, txy, sr, tt
    
    # ----------------------------------------------------
    # we use this to calculate the Airy stress function contribution
    def airy(self):
        # x,y are numpy arrays;
        # N is the order at which we truncate the sum
        # all dimensions normalized wrt cladding radius
        header = "bowtie.airy(): "
        
        x = self.x
        y = self.y
        param = self.param
        N = param['N']
        
        beta1 = param['beta1'] # beta1 - beta2
        beta3 = param['beta3'] # beta3 - beta2
        a     = param['a']     # core radius
        
        r1     = param['r1']     # bow-tie radius r1
        r2     = param['r2']     # bow-tie radius r2
        theta1 = param['theta1'] # bow-tie half-angle
        
        E     = param['E']     # Young's modulus
        nu    = param['nu']    # Poisson's ratio
        C     = param['C']     # Elasto-optic coefficient
        DT    = -param['DT']   # difference between ambient and glass melting/drawing temperature
        
        R2 = self.R2 # x**2 + y**2
        R  = self.R  # np.sqrt(R2)
        th = self.th # np.arctan2(y, x)
        
        sr = np.ones_like(x)  # sigma_r
        st = np.ones_like(x)  # sigma_theta
        tt = np.zeros_like(x) # tau_r,theta = sigma_r,theta
        
        r1sq = r1**2
        r2sq = r2**2
        th2 = 2.0*th
        theta_t = th2.copy()
        pw1 = r1sq
        pw2 = r2sq
        xth = theta1t = theta1*2
        
        tmp = theta1*(r2sq - r1sq)
        sr *= tmp
        st *= tmp
        
        for n in range(1,N+1):
            pw1 *= r1sq # $r_1^{2(n+1)}$
            pw2 *= r2sq # $r_2^{2(n+1)}$
            
            common = (pw2-pw1)*np.sin( xth )     # $(r_2^{2(n+1)} - r_1^{2(n+1)}) \sin( 2n\theta_1 )$
            fact_cos = common* np.cos( theta_t ) # $\cos(2 n \theta)$
            fact_sin = common* np.sin( theta_t ) # $\sin(2 n \theta)$
            
            m11 = float(4*n*n - 1)/float(2*n*(n+1))
            m12 = -float((n-1)*(2*n+1))/float( n*(n+1) )
            
            m21 = m11
            m22 = -float(2*n+1)/n
            
            m31 = m11
            m32 = -float(2*n+1)/(n+1)
            
            xp1 = R2**(n-1)
            xp2 = R2**n
            
            sr = sr + (m11*xp1 + m12*xp2)*fact_cos
            st = st - (m21*xp1 + m22*xp2)*fact_cos
            tt = tt - (m31*xp1 + m32*xp2)*fact_sin
            
            # next round
            theta_t = theta_t + th2 # $2 n \theta$
            xth = xth + theta1t     # $2 n \theta_1$
            
        # NOTE: we are calculating stress, but scaled to delta n; however
        # delta n is of opposite sign than stress!
        Q = C*E*DT/(1+nu) * 1e4
        ff = beta3/np.pi
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
        self.tt_a = tt # $\tau_{r,\theta} = \sigma_{r,\theta}$
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
        
        print header," calculate core contribution" # ------------------------------------
        sys.stdout.flush()
        xc = 0.0
        yc = 0.0
        # param['beta'] = param['beta1']
        sx_core, sy_core, txy_core, sr_core, tt_core = self.core()
        
        print header," calculate bow-tie SAP contribution" # ------------------------------------
        sys.stdout.flush()
        
        # print header," SAP center ",param['d2']," radius ",param['d1']
        xc = param['d2']
        yc = 0.0
        sx_sap, sy_sap, txy_sap, sr_sap, tt_sap = self.sap()
        
        # ------------------------------------
        print header," calculate Airy stress function (background), N=", N
        sys.stdout.flush()
        self.airy()
        # results: self.sx_a, self.sy_a, self.sr_a, self.tt_a, self.st_a
        
        print header," calculate total" # ------------------------------------
        sys.stdout.flush()
        sx  = sx_core + sx_sap + self.sx_a
        sy  = sy_core + sy_sap + self.sy_a
        txy = txy_core + txy_sap + self.txy_a
        
        self.sr = sr_core + sr_sap + self.sr_a
        self.tt = tt_core + tt_sap + self.tt_a
        
        # blank outside
        if blank_outside:
            # w = np.where( self.R2 > 1.0 )
            w = self.w_outside
            sx[w] = 0.0
            sy[w] = 0.0
            self.sr[w] = 0.0
            self.tt[w] = 0.0
        
        # clip sx,sy,txy values to +/- 3 sigma
        if blank_outside:
            mu    = np.mean( sx[self.w_clad] )
            sigma = np.std(  sx[self.w_clad] )
            vmin = mu - 3*sigma
            vmax = mu + 3*sigma
            w = np.where( sx > vmax )
            sx[w] = vmax
            w = np.where( sx < vmin )
            sx[w] = vmin
            
            mu    = np.mean( sy[self.w_clad] )
            sigma = np.std(  sy[self.w_clad] )
            vmin = mu - 3*sigma
            vmax = mu + 3*sigma
            w = np.where( sy > vmax )
            sy[w] = vmax
            w = np.where( sy < vmin )
            sy[w] = vmin
            
            mu    = np.mean( txy[self.w_clad] )
            sigma = np.std(  txy[self.w_clad] )
            vmin = mu - 3*sigma
            vmax = mu + 3*sigma
            w = np.where( txy > vmax )
            txy[w] = vmax
            w = np.where( txy < vmin )
            txy[w] = vmin
        
        vmin = sx.min()
        tmp  = sy.min()
        if vmin > tmp:
            vmin = tmp
        
        vmax = sx.max()
        tmp  = sy.max()
        if vmax < tmp:
            vmax = tmp
        
        # clip sr, tt to +/- 3 sigma (ignoring st)
        if blank_outside:
            mu    = np.mean( self.sr[self.w_clad] )
            sigma = np.std(  self.sr[self.w_clad] )
            vmin = mu - 3*sigma
            vmax = mu + 3*sigma
            w = np.where( self.sr > vmax )
            self.sr[w] = vmax
            w = np.where( self.sr < vmin )
            self.sr[w] = vmin
            
            mu    = np.mean( self.tt[self.w_clad] )
            sigma = np.std(  self.tt[self.w_clad] )
            vmin = mu - 3*sigma
            vmax = mu + 3*sigma
            w = np.where( self.tt > vmax )
            self.tt[w] = vmax
            w = np.where( self.tt < vmin )
            self.tt[w] = vmin
        
        vmin2 = self.sr.min()
        vmin3 = self.tt.min()
        
        vmax2 = self.sr.max()
        vmax3 = self.tt.max()
        
        # ... again...
        if blank_outside:
            w = self.w_outside
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
        
        print header," clipped x-y val range ",vmin," to ",vmax
        print header," clipped r val range ",vmin2," to ",vmax2
        print header," clipped t val range ",vmin3," to ",vmax3
        sys.stdout.flush()
        
        self.sx = sx
        self.sy = sy
        self.txy = txy
        
        if save:
            print header," saving snapshot...."
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
    def plot_rt(self, ax1=None, ax2=None,
               vmin=None, vmax=None):
        rext = 1.1 # fixed
        mycmap = plt.cm.coolwarm
        # mycmap = plt.cm.afmhot
        # mycmap = plt.cm.jet
        
        if ax1 != None:
            self.ax_rt1 = ax1
            ax1.cla()
            
            if vmin == None:
                vmin = self.vmin_r
            if vmax == None:
                vmax = self.vmax_r
            
            ax1.set_title(r'$\sigma_r$ stress (scaled to $\Delta n \times 10^{-4}$)')
            ax1.set_xlim(-rext,rext)
            ax1.set_ylim(-rext,rext)
            cax1 = ax1.pcolor(self.x, self.y, self.sr, cmap=mycmap, vmin=vmin, vmax=vmax )
            self.cbar_rt1 = ax1.figure.colorbar( cax1, ax=ax1, shrink=0.8, use_gridspec=False )
            ax1.figure.canvas.draw()
        
        if ax2 != None:
            self.ax_rt2 = ax2
            ax2.cla()
                
            if vmin == None:
                vmin = self.vmin_t
            if vmax == None:
                vmax = self.vmax_t
            
            ax2.set_title(r'$\sigma_{r\theta}$ stress (scaled to $\Delta n \times 10^{-4}$)')
            ax2.set_xlim(-rext,rext)
            ax2.set_ylim(-rext,rext)
            cax2 = ax2.pcolor(self.x, self.y, self.tt, cmap=mycmap, vmin=vmin, vmax=vmax )
            self.cbar_rt2 = ax2.figure.colorbar( cax2, ax=ax2, shrink=0.8, use_gridspec=False )
            ax2.figure.canvas.draw()
        
        return
    

# =======================================
if __name__ == '__main__':
    print "Test with test_stress2.py"
