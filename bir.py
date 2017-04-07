#!/usr/bin/python
import numpy as np;
import pylab as plt
import config

# function
def func1(x, A, B_c):
    y = (x**2)*(1 - 3.0/16.0*(x + A)**4)/( (x + A)**2 )
    return B_c*y

# derivative
def deriv1(x, A):
    y = -16*A[1] + 3*A[5] + x*( 18*A[4] + x*( 42*A[3] + x*(48*A[2] + x*(27*A[1] + x*6))))
    return y

def poly1(A):
    # pp = np.array([ 6, 27, 48, 42, 18, 3-16/A**4 ])
    pp = np.array([ 1, 4.5, 8, 7, 3, 0.5-8.0/3.0/A**4 ])
    return pp

# ------------------------------------------------------
class Design:
    def __init__(self, param, ax=None):
        header = "Design(): "
        self.param = param # save for later
        self.ax = ax
        
        Dalpha = param['alpha32'] # diff in TEC between SAP and cladding [1/deg C]
        DT     = param['DT']      # [deg C]
        E      = param['E']       # Young's modulus [N/m^2]
        C      = param['C']       # photoelastic coefficient [m^2/N]
        nu     = param['nu']      # Poisson ratio
        
        B_c = 2.0*C*E*Dalpha*DT/ ( 1.0 - nu ) * 1e4
        # print "# B_c = %.4e" % B_c;
        
        # geometry
        ab = param['a']   # core radius a/b
        
        x1 = []
        y1 = []
        
        self.ralist = ralist = np.arange(3.0, 5.5, 0.5)
        idx = 0
        for ra in ralist:
            # r = d2 - d1 = sap_center - sap_radius
            # ra = r/a
            A = 2.0*ra*ab;
            # Alist = [ 1.0, A, A*A, A*A*A, A*A*A*A, A*A*A*A*A ]
            
            print header," calculate birefringe for r/a = ",ra;
            x1.append( np.arange(0, 0.9, 0.01) )
            y1.append( np.array([func1( x, A, B_c) for x in x1[idx]]) )
            idx += 1
        
        # self.brcurves = np.column_stack( (x1, y1) )
        self.brcurves_x = x1
        self.brcurves_y = y1
        
        # ------- max curve ---------
        x2 = []
        y2 = []
        print header," calculate maxima curve"
        for ra in np.arange(2.8, 5.5, 0.1 ):
            A = 2.0*ra*ab;
            pp = poly1(A)
            roots = np.roots( pp )
            x = roots[4].real * A
            y = func1(x, A, B_c);
            x2.append( x )
            y2.append( y )
        
        self.maxvalues_x = np.array(x2)
        self.maxvalues_y = np.array(y2)
        
        
        # --- 95% lines --------------
        idx = 0
        x3 = []
        y3 = []
        for ra in ralist:
            A = 2.0*ra*ab;
            pp = poly1(A)
            roots = np.roots( pp )
            x = roots[4].real * A
            y = func1(x, A, B_c) # peak value
            y = y*0.95
            x3.append( np.arange( x-0.2, x+0.2, 0.01 ) )
            y3.append( np.zeros_like( x3[idx] ) + y )
            idx += 1
        
        self.lines95_x = x3
        self.lines95_y = y3
        
        return
    
    # ---- optimum pick ----------------
    def best(self, r):
        """
        Given the distance from fiber center to the edge of SAP (r), calculate the optimum SAP radius.
        The result is independent of core radius a; all distances (including the answer) are
        normalized with respect to the cladding radius (i.e., b=1)
        """
        # a = self.param['a'] # core_radius
        # r = sap_center - sap_radius
        # ra = r/a
        # ab = a
        # A = 2.0*ra*ab;
        A = 2.0*r
        
        pp = poly1(A)
        roots = np.roots( pp )
        x = roots[4].real * A
        
        # t = sap_radius * 2.0 (diameter)
        # best radius = t/2
        return x/2.0
    
    # ----- export -------------------------
    def export(self, filename):
        # TBD
        return
    
    # ----- plotting -------------------------
    def plot(self):
        header = "bir.Design.plot(): "
        if self.ax == None:
            print header," plotting axis not defined; aborting"
            return
        
        print header," generating birefingence vs. SAP diameter plot"
        self.ax.cla()
        self.ax.set_xlabel("Normalized SAP diameter t/b")
        self.ax.set_ylabel(r'Modal birefringence B ($\times 10^{-4}$)')
        self.ax.set_title("Birefringence at core center")
        
        idx = 0
        for idx in range(len(self.ralist)):
            x = self.brcurves_x[idx]
            y = self.brcurves_y[idx]
            self.ax.plot(x,y, label=("r/a = %.1f" % self.ralist[idx]))
            
            x = self.lines95_x[idx]
            y = self.lines95_y[idx]
            self.ax.plot(x,y,'--',color='black')
        
        
        self.ax.plot(self.maxvalues_x, self.maxvalues_y, '-.', color='red')
        self.ax.legend(loc='best')
        
        self.ax.figure.canvas.draw()
        
        return

# ===============================================================
if __name__ == '__main__':
    import sys
    
    cfgfile = "lib/test1.ini"
    print "TEST: Using configuration file ",cfgfile
    sys.stdout.flush()
    param = config.get(cfgfile)
    print "TEST: E = %.1e" % param['E']
    print "TEST: nu = %.3f" % param['nu']
    print "TEST: Delta alpha = %.3g" % param['alpha32']
    print "TEST: C = %.3g" % param['C']
    print "TEST: DT = %.1f" % param['DT']
    print "TEST: a = %.1f" % param['a']
    
    print "TEST: setting up plotting area..."
    sys.stdout.flush()
    
    fig1 = plt.figure("Birefringence", figsize=(10,6))
    ax1 = fig1.add_subplot(111)
    
    print "TEST: calculating..."
    sys.stdout.flush()
    
    bir = Design(param, ax=ax1)
    
    r = 3.0*param['a']
    tb = 2.0*bir.best( r ) # NOTE: return is radius!
    print "TEST: est best pick for r = 3a is t/b = %.3f" % tb
    
    print "TEST: calculations done! Now plotting..."
    sys.stdout.flush()
    
    bir.plot()
    
    plt.show()

