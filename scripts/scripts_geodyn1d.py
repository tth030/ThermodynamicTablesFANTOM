#!/usr/bin/env python

from . import geodyn1d
from .geodyn1d.tools import *
from .book_deltarho import *
from copy import deepcopy
import sys

def compute_rhoref_along_geotherm(perplex,
                                  lithos=None,
                                  temperature_shift=0,
                                  use_readrho=False):
    ''' 
        Compute reference density, i.e. at surface conditions (use_readrho=False), with its standard deviation.
        along a given geotherm based
    '''
        
    Tref  = 273
    Pref  = 0
    rho   = deepcopy(perplex.rho)
    T     = deepcopy(perplex.T)
    P     = deepcopy(perplex.P)
    alpha = deepcopy(perplex.alpha)
    beta  = deepcopy(perplex.beta)
    melt  = deepcopy(perplex.melt)
    
    factor1 = 1 + np.multiply(beta,P-Pref)
    factor2 = 1 - np.multiply(alpha,T-Tref)
    factor  = np.multiply(factor1,factor2)
    rhoref = np.divide(rho,factor)
    if (use_readrho):
        rhoref = rho
    
    if ( lithos is None ):
        print('Lithosphere lith125 used as default to compute rhoref')
        lithos = geodyn1d.lithosphere('lith125')

    if ( isinstance(lithos, str) ):
        lithos = geodyn1d.lithosphere(lithos)
    pmax = np.max(perplex.P)
    tmax = np.max(perplex.T)
    T = np.linspace(np.min(perplex.T),
                    tmax,
                    perplex.nt)
    P = np.linspace(np.min(perplex.P),
                    pmax,
                    perplex.np)
    f      = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), rhoref.reshape(perplex.np,perplex.nt))
    f_melt = RegularGridInterpolator((np.array(P)/pmax, np.array(T)/tmax), perplex.melt.reshape(perplex.np,perplex.nt))
    
    # we check for geotherm and pressure profile
    if ( not lithos.geotherm):
        lithos.get_steady_state_geotherm(printOut=False)
    if ( not lithos.P_rho_DepthProfile):
        lithos.get_pressure_density_profile(printOut=False)
    
    sort_rhoref = []
    sort_melt     = []
    for i in range(0,len(lithos.geotherm.Z)):
        pressure = lithos.P_rho_DepthProfile.P[i]
        temp     = lithos.geotherm.T[i] + temperature_shift
        melt_read = f_melt([pressure/pmax,temp/tmax])
        if (melt_read>0):
            sort_rhoref.append(f([pressure/pmax,temp/tmax]))
            sort_melt.append(melt_read)
    sort_melt = [i * 100 for i in sort_melt]
        
    return sort_melt, sort_rhoref


def update_mantle_density_to_reference(name,book,perplex_name='slm_sp2008',rhoc0=2860):
    '''
       Update mantle parameters to use perplex grids to compute density
    '''
    newbook    = book.copy()
    layers = { 'matLM1', 'matLM2', 'matLM3', 'matSLMd', 'matSLM' }

    if perplex_name:
        newbook = update_book(newbook,'use_tidv',7,['LM1','LM2','LM3','SLMd','SLM'])
        newbook = update_perplex_name(layers,newbook,perplex_name)
        newbook = update_deltarho(layers,newbook,name,perplex_name,use_perplex=True,rhoc0=rhoc0) # should be run after update_perplex_name if use_perplex is True
    else:
        newbook = update_book(newbook,'use_tidv',1,['LM1','LM2','LM3','SLMd','SLM'])
        newbook = update_deltarho(layers,newbook,name,perplex_name,use_perplex=False,rhoc0=rhoc0)
        newbook = update_book(newbook,'rho',3311,['LM1','LM2','LM3','SLMd','SLM'])
    return newbook


def update_perplex_name(layers,book,perplex_name):
    '''
      name: perplex grid name that refers to the composition
    '''
    newbook = book.copy()
    for layer in layers:
        if layer in newbook:
            if layer == 'matSLM':
                newbook[layer]['perplex_name'] = 'slm_sp2008' # default reference fertile sub-lithospheric mantle composition
            else:
                newbook[layer]['perplex_name'] = perplex_name
        #else:
        #    print('No "{}" in this book'.format(layer))
    return newbook

def update_deltarho(layers,book,name,perplex_name='slm_sp2008',use_perplex=False,rhoc0=2860):
    '''
      name: lithosphere name
    '''
    newbook = book.copy()
    for layer in layers:
        if layer == 'matSLM':
            deltarho = 0
        else:
            deltarho = find_drho(perplex_name,name,use_perplex,rhoc0)
        if layer in newbook:
            newbook[layer]['deltarho'] = deltarho
        #else:
        #    print('No "{}" in this book'.format(layer))
    return newbook

def find_drho(perplex_name='slm_sp2008',name='lith125',use_perplex=False,rhoc0=2860):
    '''
       type = 'perplex' or 'onlytempdep'
    '''
    book=book_deltarho # database valuable for rhoc0 == 2860 kg/m3 
    if rhoc0==2860:
        if use_perplex:
            if name not in book['perplex']:
                print("ERROR: Lithosphere name is unknown ("+name+")")
                sys.exit()
            if perplex_name in book['perplex'][name]:
                deltarho = book['perplex'][name][perplex_name]
            else:
                print("WARNING: reference deltarho for this composition ("+perplex_name+") and this lithosphere ("+name+") is unknown. We use deltarho = 0")
                deltarho = 0
        else:
            if name not in book['onlytempdep']:
                print("ERROR: Lithosphere name is unknown ("+name+")")
                sys.exit()
            deltarho = book['onlytempdep'][name]
    else:
        print("ERROR: Not implemented yet with different crustal density than 2860 kg/m3")
        sys.exit()
    return deltarho

def update_book(book,name,value,layers):
    newbook = book.copy()
    for layer in layers:
        if 'mat'+layer in newbook:
            newbook['mat'+layer][name] = value
        #else:
        #    print('No "{}" in this book'.format('mat'+layer))
    return newbook

