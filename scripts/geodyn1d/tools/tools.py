#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 12:15:25 2018

@author: thomastheunissen

----------------------------------


"""
from scipy.interpolate import interp1d
#from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import os,errno
from scipy.interpolate import interp2d
from scipy.interpolate import InterpolatedUnivariateSpline

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def pos2value_1d(x,v,xs):
    """ position to value 1D uses
     a spline interpolation 2nd order. """
    f = interp1d(x, v, kind='quadratic')
    return f(xs)

def pos2value_2d(x,y,v,xs,ys):
    """ position to value 2D uses
     a linear interpolation around xs. """
    sampling = 20
    dx = abs(x[0,1]-x[0,0])
    dy = abs(y[1,0]-y[0,0])
    x_tmp = x[0,:]
    y_tmp = y[:,0]
    nx = int(len(x_tmp[(x_tmp>=xs-sampling*dx) & (x_tmp<=xs+sampling*dx)]))
    ny = int(len(y_tmp[(y_tmp>=ys-sampling*dy) & (y_tmp<=ys+sampling*dy)]))
    
    if ((nx>0)&(ny>0)):
        x_small = x[ (x>=xs-sampling*dx) &
                     (x<=xs+sampling*dx) &
                     (y>=ys-sampling*dy) &
                     (y<=ys+sampling*dy) ].reshape(ny,nx)
        y_small = y[ (x>=xs-sampling*dx) &
                     (x<=xs+sampling*dx) &
                     (y>=ys-sampling*dy) &
                     (y<=ys+sampling*dy) ].reshape(ny,nx)
        v_small = v[ (x>=xs-sampling*dx) &
                     (x<=xs+sampling*dx) &
                     (y>=ys-sampling*dy) &
                     (y<=ys+sampling*dy) ].reshape(ny,nx)
        
        vmin = np.max(v_small)
        vmax = np.max(v_small)
        if (isclose(vmin,vmax)):
            return vmin
        
        # we normalise
        xmax = np.max(x_small)
        ymax = np.max(y_small)
        x_small = x_small / xmax
        y_small = y_small / ymax
        #f = interp2d(x_small, y_small, v_small, kind='linear')
        f = Rbf(x_small, y_small, v_small, function='linear')
        return f(xs/xmax,ys/ymax)
    else:
        return -9.9e12
    
def pos2value_regular2d(x,y,v,xs,ys):
    """ position to value 2D avoids expensive triangulation of the input data 
    by taking advantage of the regular grid structure.. """
    
    xmax = np.max(x)
    ymax = np.max(y) 
        
    f = RegularGridInterpolator((np.array(x)/xmax, np.array(y)/ymax), v)
    
    return f([xs/xmax,ys/ymax])[0]

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

