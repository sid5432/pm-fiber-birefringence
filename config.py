#!/usr/bin/env python
from __future__ import print_function
# may need these for pyinstaller
import six
import numpy as np
# import tkMessageBox
# import FileDialog

import sys
if sys.version_info > ( 3, 0 ):
    import tkinter.filedialog as tkf
else:
    import tkFileDialog as tkf

import iniparse

iniFormats = [ 
        ('All Files','*.*'),
        ('INI Files','*.ini'),
        ]

def get(filename=None):
    header = "config.get(): "
    # open INI file
    if filename == None:
        filename = tkf.askopenfilename(filetypes=iniFormats, title='Open parameters file')
    
    if len(filename) <= 0:
        print( header," * Cancel" )
        return None
    
    print( header," loading configuration from ", filename )
    cfg = iniparse.INIConfig( open(filename) )
    
    # NOTE: core is region 1, SAP's are region 3, everywhere else (cladding) is region 2
    param = {}
    param['DT'] = float(cfg.material.DT) # temperature diff between ambient and melt/draw
    param['E']  = float(cfg.material.E)  # Young's modulus of glass [N/m^2]
    param['nu'] = float(cfg.material.nu) # Poisson's ratio
    # param['C']  = float(cfg.material.C)  # photoelastic coefficient [m^2/N]
    param['C1']  = float(cfg.material.C1)  # C1=B2; photoelastic coefficient [m^2/N]
    param['C2']  = float(cfg.material.C2)  # C2=B1; photoelastic coefficient [m^2/N]
    
    param['C'] = param['C2'] - param['C1']

    param['alpha12'] = float(cfg.material.alpha12) # difference in the thermal elasticity coefficient between core and cladding
    param['alpha32'] = float(cfg.material.alpha32) # difference in the thermal elasticity coefficient between SAP and cladding
    
    # derived
    param['beta1'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha12'] # core
    param['beta3'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha32'] # SAP
    
    # geometry parameters
    param['b']  = float(cfg.geometry.cladding_radius)  # cladding radius (the only one not normalized)
    param['a']  = float(cfg.geometry.core_radius)  # core radius (normalized wrt cladding radius)
    param['sap_type'] = cfg.geometry.sap_type
    
    # panda: we get and save all parameters, but don't use all of them
    param['d1'] = float(cfg.geometry.sap_radius)   # SAP radius (normalized wrt cladding radius)
    param['d2'] = float(cfg.geometry.sap_center)   # center of SAP to center of fiber/core (normalized wrt cladding radius)
    
    # bow-tie
    param['r1'] = float(cfg.geometry.r1) # small radius of bow-tie (normalized wrt cladding radius)
    param['r2'] = float(cfg.geometry.r2) # large radius of bow-tie (normalized wrt cladding radius)
    param['theta1'] = float(cfg.geometry.theta1) # half-angle extent of bow-tie (degrees)
    param['theta1'] *= np.pi/180.0 # convert to radians
    
    if param['sap_type'] != 'bow-tie' and param['sap_type'] != 'panda':
        print( header,"!!! Unrecognized SAP type ",param['sap_type'])

    # others
    param['N'] = int(cfg.calc.N) # Number of terms to take in Airy function summation
    param['gridsize'] = int(cfg.calc.gridsize) # size of grid for computation; should be an odd integer
    param['subsample'] = int(cfg.calc.subsample) # subsample (skips) of grid for tensor visualization
    
    return param

# ----------------------------------------------------------------------
def save(param, filename=None):
    """save param to file"""
    header = "config.save(): "
    
    cfg = iniparse.INIConfig()
    
    cfg.material.DT = param['DT'] # temperature diff between ambient and melt/draw
    cfg.material.E  = param['E']  # Young's modulus of glass [m^2/M]
    cfg.material.nu = param['nu'] # Poisson's ratio
    cfg.material.C1 = param['C1']  # photoelastic coefficient [m^2/N = Pa]
    cfg.material.C2 = param['C2']  # photoelastic coefficient [m^2/N = Pa]
    
    cfg.material.alpha12 = param['alpha12'] # difference in the thermal elasticity coefficient between core and cladding
    cfg.material.alpha32 = param['alpha32'] # difference in the thermal elasticity coefficient between SAP and cladding
    
    # derived
    # param['beta1'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha12'] # core
    # param['beta3'] = (1 + param['nu'])/(1 - param['nu']) * param['alpha32'] # SAP
    
    # geometry parameters
    cfg.geometry.cladding_radius = param['b']  # cladding radius in um (the only one not normalized!)
    cfg.geometry.core_radius = param['a']   # core radius (normalized wrt cladding radius)
    
    cfg.geometry.sap_type = param['sap_type']
    
    cfg.geometry.sap_radius  = param['d1']  # SAP radius (normalized wrt cladding radius)
    cfg.geometry.sap_center  = param['d2']  # center of SAP to center of fiber/core (normalized wrt cladding radius)
    cfg.geometry.r1 = param['r1']
    cfg.geometry.r2 = param['r2']
    # convert to degrees
    cfg.geometry.theta1 = param['theta1'] * 180.0/np.pi
    
    if param['sap_type'] != 'panda' and param['sap_type'] != 'bow-tie':
        print( header,"!!! Unrecognized SAP type ",param['sap_type'] )
    
    # order
    cfg.calc.N = param['N'] # Number of terms to take in Airy function summation
    cfg.calc.gridsize = param['gridsize'] # size of grid for computation; should be an odd integer
    cfg.calc.subsample = param['subsample'] # subsample/skip of grid in tensor visualization
    
    # save to file
    if filename == None:
        filename = tkf.asksaveasfilename(filetypes=iniFormats, title='Save parameters to file')
    
    if len(filename) <= 0:
        print( header,"* Cancel" )
        return False
    
    try:
        f = open(filename,'w')
        print( cfg, file=f )
        f.close()
    except IOError:
        print( header,"! can not save parameters to file" )
    
    return True

# =============================================================================
if __name__ == '__main__':
    filename = "lib/default.ini"
    
    param = get(filename)
    
    print( "TEST: material" )
    print( "TEST:\tE  = %.3g" % param['E'] )
    print( "TEST:\tC1 = %.3g" % param['C1'] )
    print( "TEST:\tC2 = %.3g" % param['C2'] )
    
    print( "TEST:\tSAP type = %s" % param['sap_type'] )
    
    status = save(param)
    
