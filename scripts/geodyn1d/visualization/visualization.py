# some of visualization routines are from "pyFantom" by
# Romain Beucher, Sebastian Wolf, Leo Zijerveld
# https://bitbucket.org/SebWolf/pyfantom/src/master
#
#

# ---------------------------------------------------------------
# -------Definition of global variables--------------------------
# ---------------------------------------------------------------

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from functools import wraps
from .. import core as geodyn1d

from ..geotherm import *
from ..isostasy import *
from ..thermodyn import *
from ..tools import *


# ---------------------------------------------------------------
# -------decorators----------------------------------------------
# ---------------------------------------------------------------

def add_subplt():
    def wrap(f):
        @wraps(f)
        def new_function(self,*args, **kwargs):
            if 'inner' in kwargs :
                pass
            else:
                if 'subplot' in kwargs and 'fig' in kwargs :
                    subplot = kwargs["subplot"]
                    fig = kwargs["fig"]
                    if 'sharexaxis' in kwargs and 'axis2share' in kwargs:
                        axis2share = kwargs['axis2share']
                        if len(str(subplot)) == 4:
                            ax = fig.add_subplot(int(str(subplot)[0]),int(str(subplot)[1]),int(str(subplot)[2::]), sharex = axis2share)#,adjustable='box'
                        else:
                            ax = fig.add_subplot(subplot, sharex = axis2share,adjustable='box-forced')#,adjustable='box'
                    else:
                        if len(str(subplot)) == 4:
                            ax = fig.add_subplot(int(str(subplot)[0]),int(str(subplot)[1]),int(str(subplot)[2::]))#,adjustable='box'
                        else:
                            ax = fig.add_subplot(subplot)#,adjustable='box'
                    #print('Subplot created in fig \n')
                # this is from pyFantom, no need of those conditions here now
                #elif not 'material_colors' in kwargs and not 'PlotInThisAx' in kwargs:
                else:
                    #fig = plt.figure(frameon=False) # ,figsize=[8,6]
                    #fig = plt.figure()
                    #ax  = fig.add_subplot(111)
                    fig, ax = plt.subplots()
                    #print('Plotting without subplots \n')
            a = f(self,*args, **kwargs)
            return a
        return new_function
    return wrap


@ add_subplt()
def plot_points(perplex,
               data = 'rho',
               xlim=[0,1650], ylim=[0,6e9],
               contours = None,
               vmin = None, vmax = None, colorbar="right",cmap = "jet",
               display_title=True,printOut=False,display_contour = True,
               geotherms=False,lithospheres=None,
               **kwargs):
 
    #T  = np.linspace(perplex.T[0],max(perplex.T),perplex.nt)
    #T  = T - 273
    #P  = np.linspace(perplex.P[0],max(perplex.P),perplex.np)
    T = perplex.T -273.15
    P = perplex.P

    if (data=='rhoref'):
        Tref  = 273
        Pref  = 0
        factor1        = 1 + np.multiply(perplex.betaresidual,perplex.P-Pref)
        factor2        = 1 - np.multiply(perplex.alpharesidual,perplex.T-Tref)
        factor         = np.multiply(factor1,factor2)
        perplex.rhoref = np.divide(perplex.rhoresidual,factor)
    
    book = {'rho':  perplex.rho,
            'rhoref':  perplex.rhoref,
            'alpha':perplex.alpha,
            'beta': perplex.beta,
            'cp': perplex.cp,
            'melt': perplex.melt,
            'deltarho': perplex.deltarho,
            'rhomelt':  perplex.rhomelt,
            'rhoresidual':  perplex.rhoresidual,
            'alphamelt':perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt': perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt': perplex.cpmelt,
            'cpresidual': perplex.cpresidual }
    
    if ( not data in book ):
        print('data "{}" does not belong to perplex object'.format(data))
        return
    
    value = deepcopy(book[data])
    print(value)
    
    book_scaling = {'rho':  1.0e0,
                    'rhoref':  1.0e0,
                    'alpha':1e5,
                    'beta': 1e11,
                    'cp': 1.0e0,
                    'melt': 1.0e2,
                    'deltarho': 1.0e0,
                    'rhomelt':  1.0e0,
                    'rhoresidual':  1.0e0,
                    'alphamelt':1e5,
                    'alpharesidual':1e5,
                    'betamelt': 1e11,
                    'betaresidual': 1e11,
                    'cpmelt': 1.0e0,
                    'cpresidual': 1.0e0}
    value = value * book_scaling[data]
    
    #if ( not perplex.rho_ref ):
    #    lith1 = geodyn1d.lithosphere('lith125')
    #    perplex.rho_ref, rho_ref_std, rho_ref = compute_rho_ref_perplex(perplex,domainBased=False,lithosphere=lith1)
    book_colorbar_label = {'rho':  "Density ($kg.m^{-3}$)",
                           'rhoref':  "Density ($kg.m^{-3}$)",
                           'alpha':"Thermal expansion x$10^{5}$ ($1/K$)",
                           'beta': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'cp': "Specific Heat Capacity ($J/K/Kg$)",
                           'melt': "Melt percentage (%)",
                           'deltarho': "Density variations ($kg.m^{-3}$) \n $rho_{ref}=$ "+str(perplex.rho_ref),
                           'rhomelt':  "Density ($kg.m^{-3}$)",
                           'rhoresidual':  "Density ($kg.m^{-3}$)",
                           'alphamelt': "Thermal expansion x$10^{5}$ ($1/K$)",
                           'alpharesidual': "Thermal expansion x$10^{5}$ ($1/K$)",
                           'betamelt': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'betaresidual': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'cpmelt': "Specific Heat Capacity ($J/K/Kg$)",
                           'cpresidual': "Specific Heat Capacity ($J/K/Kg$)"}
    colorbar_label = book_colorbar_label[data]
    
    book_fmt = {'rho':  '%d',
                'rhoref':  '%d',
                'alpha':'%3.1f',
                'beta': '%3.1f',
                'cp': '%6.1f',
                'melt': '%d',
                'deltarho': '%d',
                'rhomelt':  '%d',
                'rhoresidual':  '%d',
                'alphamelt': '%3.1f',
                'alpharesidual': '%3.1f',
                'betamelt': '%3.1f',
                'betaresidual': '%3.1f',
                'cpmelt': '%6.1f',
                'cpresidual': '%6.1f'}
    fmt = book_fmt[data]
   
    #raster = plt.pcolormesh(T,P,value,rasterized=True, vmin=vmin, vmax=vmax,cmap = cmap,shading='nearest')
    raster = plt.scatter(T,P,s=5,c=value,rasterized=True, vmin=vmin, vmax=vmax,cmap = cmap)
    
    if (geotherms):
        if ( not lithospheres):
            print('Cannot display any geotherms. "Lithospheres" argument is required.')
        else:
            for lith in lithospheres:
                lith1 = geodyn1d.lithosphere(lith=lith)
                lith1.get_steady_state_geotherm(printOut=False)
                lith1.get_pressure_density_profile(printOut=False)
                plt.plot(lith1.geotherm.T-273.15,lith1.P_rho_DepthProfile.P)
    
    ax = plt.gca()
    
    if(colorbar):
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        orientation = {"left": "vertical", "right": "vertical",
                       "bottom": "horizontal", "top": "horizontal"}
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes(colorbar, size="5%", pad=0.05)
        #cb  = plt.colorbar(orientation=orientation[colorbar], cax=cax,
        #                  label=colorbar_label)
        #cb.solids.set_rasterized(True)
        plt.colorbar(raster,orientation=orientation[colorbar],label=colorbar_label)
    
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    plt.gca().invert_yaxis()
    
    return ax



@ add_subplt()
def plot_table(perplex,
               data = 'rho',
               xlim=[0,1650], ylim=[0,6e9],
               contours = None,
               vmin = None, vmax = None, colorbar="right",cmap = "jet",
               display_title=True,printOut=False,display_contour = True,
               geotherms=False,lithospheres=None,book_lith=None,
               **kwargs):
 
    #dT = (max(perplex.T)-perplex.T[0])/(perplex.nt-1)
    T  = np.linspace(perplex.T[0],max(perplex.T),perplex.nt)
    T  = T - 273
    #dP = (max(perplex.P)-perplex.P[0])/(perplex.np-1)
    P  = np.linspace(perplex.P[0],max(perplex.P),perplex.np)
   
    #if (data=='cp'):
    #    print('bah alaors, ou est cp ??.............')
    #    print( perplex.cp )
 
    if (data=='rhoref'):
        Tref  = 273
        Pref  = 0
        factor1        = 1 + np.multiply(perplex.beta,perplex.P-Pref)
        factor2        = 1 - np.multiply(perplex.alpha,perplex.T-Tref)
        factor         = np.multiply(factor1,factor2)
        perplex.rhoref = np.divide(perplex.rhoresidual,factor)
    
    book = {'rho':  perplex.rho,
            'rhoref':  perplex.rhoref,
            'alpha':perplex.alpha,
            'beta': perplex.beta,
            'cp': perplex.cp,
            'melt': perplex.melt,
            'deltarho': perplex.deltarho,
            'rhomelt':  perplex.rhomelt,
            'rhoresidual':  perplex.rhoresidual,
            'alphamelt':perplex.alphamelt,
            'alpharesidual':perplex.alpharesidual,
            'betamelt': perplex.betamelt,
            'betaresidual': perplex.betaresidual,
            'cpmelt': perplex.cpmelt,
            'cpresidual': perplex.cpresidual }
    
    if ( not data in book ):
        print('data "{}" does not belong to perplex object'.format(data))
        return
    
    value = deepcopy(book[data]).reshape(perplex.np, perplex.nt)
    
    book_scaling = {'rho':  1.0e0,
                    'rhoref':  1.0e0,
                    'alpha':1e5,
                    'beta': 1e11,
                    'cp': 1.0e0,
                    'melt': 1.0e2,
                    'deltarho': 1.0e0,
                    'rhomelt':  1.0e0,
                    'rhoresidual':  1.0e0,
                    'alphamelt':1e5,
                    'alpharesidual':1e5,
                    'betamelt': 1e11,
                    'betaresidual': 1e11,
                    'cpmelt': 1.0e0,
                    'cpresidual': 1.0e0}
    value = value * book_scaling[data]
    
    book_contours = {'rho':  [3100,3150,3200,3250,3300,3350,3400],
                     'rhoref':  [3100,3150,3200,3225,3250,3275,3300,3325,3350,3375,3400,3425,3450],
                     'alpha':[2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6],
                     'beta': [0.7,0.8,0.9,1.0,1.1,1.2],
                     'cp': [600,650,700,750,800,850,900],
                     'melt': [10,20,30,40,50,60],
                     'deltarho':     [-125,-100,-75,-50,-25,0,25,50,75,100,125],
                     'rhomelt':      [2700,2750,2800,2850,2900,2950,3000],
                     'rhoresidual':  [3100,3150,3200,3250,3300,3350,3400],
                     'alphamelt':    [4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5],
                     'alpharesidual':[2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6],
                     'betamelt':     [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],
                     'betaresidual': [0.7,0.8,0.9,1.0,1.1,1.2],
                     'cpmelt': [600,650,700,750,800,850,900],
                     'cpresidual': [600,650,700,750,800,850,900]}
    contours = book_contours[data]
    
    book_manual_contour = {'rho':  False,
                           'rhoref':  False,
                           'alpha':False, #[ (t,2e9) for t in tcontour ] ,
                           'beta': False,
                           'cp': False,
                           'melt': False,
                           'deltarho': False,
                           'rhomelt':  False,
                           'rhoresidual':  False,
                           'alphamelt': False,
                           'alpharesidual': False,
                           'betamelt': False,
                           'betaresidual': False,
                           'cpmelt': False,
                           'cpresidual': False}
    manual_contour = book_manual_contour[data]
    
    if (isclose(np.min(value),np.max(value))):
        display_contour = False
    
    book_colorbar_label = {'rho':  "Density ($kg.m^{-3}$)",
                           'rhoref':  "Density ($kg.m^{-3}$)",
                           'alpha':"Thermal expansion x$10^{5}$ ($1/K$)",
                           'beta': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'cp': "Specific Heat Capacity ($J/K/Kg$)",
                           'melt': "Melt percentage (%)",
                           'deltarho': "Density variations ($kg.m^{-3}$) \n $rho_{ref}=$ "+str(perplex.rho_ref),
                           'rhomelt':  "Density ($kg.m^{-3}$)",
                           'rhoresidual':  "Density ($kg.m^{-3}$)",
                           'alphamelt': "Thermal expansion x$10^{5}$ ($1/K$)",
                           'alpharesidual': "Thermal expansion x$10^{5}$ ($1/K$)",
                           'betamelt': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'betaresidual': "Compressibility x$10^{11}$ ($1/Pa$)",
                           'cpmelt': "Specific Heat Capacity ($J/K/Kg$)",
                           'cpresidual': "Specific Heat Capacity ($J/K/Kg$)"}
    colorbar_label = book_colorbar_label[data]
    
    book_fmt = {'rho':  '%d',
                'rhoref':  '%d',
                'alpha':'%3.1f',
                'beta': '%3.1f',
                'cp': '%6.1f',
                'melt': '%d',
                'deltarho': '%d',
                'rhomelt':  '%d',
                'rhoresidual':  '%d',
                'alphamelt': '%3.1f',
                'alpharesidual': '%3.1f',
                'betamelt': '%3.1f',
                'betaresidual': '%3.1f',
                'cpmelt': '%6.1f',
                'cpresidual': '%6.1f'}
    fmt = book_fmt[data]
    
    raster = plt.pcolormesh(T,P,value,rasterized=True, vmin=vmin, vmax=vmax,cmap = cmap,shading='nearest')
    
    if (geotherms):
        if ( not lithospheres):
            print('Cannot display any geotherms. "Lithospheres" argument is required.')
        else:
            i = 0
            for lith in lithospheres:
                lith1 = geodyn1d.lithosphere(lith=lith)
                if (book_lith):
                    lith1.update_materials(book_lith[i])
                lith1.get_steady_state_geotherm(printOut=False)
                lith1.get_pressure_density_profile(printOut=False)
                plt.plot(lith1.geotherm.T-273.15,lith1.P_rho_DepthProfile.P)
                i += 1
    
    ax = plt.gca()
    
    # contours
    if (display_contour):
        cont = ax.contour(T, P, value,contours , colors="black",linewidths=1.)
        if (np.nanstd(value)>0):
            plt.clabel(cont, inline=1,fmt=fmt,manual=manual_contour)

    book_title = {'rho':  "Density ($kg.m^{-3}$)",
                  'rhoref':  "Density $\rho(P,T)_{0}$ ($kg.m^{-3}$)",
                  'alpha':"Thermal expansion x$10^{5}$ ($1/K$)",
                  'beta': "Compressibility x$10^{11}$ ($1/Pa$)",
                  'cp': "Specific Heat Capacity ($J/K/Kg$)",
                  'melt': "Melt percentage (%)",
                  'deltarho': "Density variations with respect to an average ($kg.m^{-3}$) \n $rho_{ref}=$ "+str(perplex.rho_ref),
                  'rhomelt':  "Density of the melt ($kg.m^{-3}$)",
                  'rhoresidual':  "Density of the residual ($kg.m^{-3}$)",
                  'alphamelt': "Thermal expansion of the melt x$10^{5}$ ($1/K$)",
                  'alpharesidual': "Thermal expansion of the residual x$10^{5}$ ($1/K$)",
                  'betamelt': "Compressibility of the melt x$10^{11}$ ($1/Pa$)",
                  'betaresidual': "Compressibility of the residual x$10^{11}$ ($1/Pa$)",
                  'cpmelt': "Specific Heat Capacity of the melt ($J/K/Kg$)",
                  'cpresidual': "Specific Heat Capacity of the residual ($J/K/Kg$)"}
 
    if display_title:
        plt.title(book_title[data])

    if(colorbar):
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        orientation = {"left": "vertical", "right": "vertical",
                       "bottom": "horizontal", "top": "horizontal"}
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes(colorbar, size="5%", pad=0.05)
        #cb  = plt.colorbar(orientation=orientation[colorbar], cax=cax,
        #                  label=colorbar_label)
        #cb.solids.set_rasterized(True)
        plt.colorbar(raster,orientation=orientation[colorbar],label=colorbar_label)
    
    plt.xlim(xlim)
    plt.ylim(ylim)
    
    plt.gca().invert_yaxis()
    
    return ax

@ add_subplt()
def plot_delta_rho_vertical_gradient_perplex(perplex,
                  xlim=[0,1650], ylim=[0,6e9],
                  vmin = 0, vmax = 10, colorbar="right",cmap = "jet",
                  display_title=True,printOut=False,
                  ax=None,fig=None,subplot=None,
                  **kwargs):
    
    compute_delta_rho_perplex(perplex)
    deltarho_vert_gradient = compute_delta_rho_vertical_gradient_perplex(perplex)
    
    ax = plt.gca()
    plt.pcolormesh(perplex.T.reshape(perplex.np,perplex.nt)-273,
                   perplex.P.reshape(perplex.np,perplex.nt),
                   deltarho_vert_gradient,
                   rasterized=True, cmap = cmap,vmin=vmin,vmax=vmax,shading='nearest')
    plt.ylim(ylim)
    plt.xlim(xlim)
    
    colorbar_label = "Vertical "+r"$\Delta\rho$"+" gradient"
    cax = None
    
    if(colorbar):
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        orientation = {"left": "vertical", "right": "vertical",
                       "bottom": "horizontal", "top": "horizontal"}
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(colorbar, size="5%", pad=0.05)
        cb = plt.colorbar(orientation=orientation[colorbar], cax=cax,
                          label=colorbar_label)
        cb.solids.set_rasterized(True)
    
    return ax

@ add_subplt()
def plot_density_depth_profile(rho,Z,xlim=None, ylim=None,
                               display_xlabel=True,
                               display_title=True,title_in_plot=False,description='unknown lithos.',
                               printOut=False,*args, **kwargs):
    """ Plot density (Kg/m3) - depth (Km) profile """
    
    plt.plot(rho,np.multiply(Z,1e-3))
    ax = plt.gca()
    
    if ( display_title ):
        title = 'Density profile ({})'
        plt.title(title.format(description))
        if title_in_plot:
            title = '{}'
            ax.text(0.97, 0.1, title.format(description),
                verticalalignment='center', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=16, bbox =dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.5) )#{'facecolor':'white', 'alpha':0.6, 'pad':5}), pad=5
            plt.title("")
    if ( display_xlabel ):
        plt.xlabel("Density ($kg.m^{-3}$)")
    plt.ylabel("Depth (km)")
    if (xlim):
        plt.xlim(xlim)
    if (ylim):
        plt.ylim(ylim)
        
    ax.invert_yaxis()
            
    return ax

@ add_subplt()
def plot_temperature_depth_profile(T,Z,xlim=None, ylim=None,
                        display_title=True,title_in_plot=False,description='unknown lithos.',
                        printOut=False,*args, **kwargs):
    """ Plot temperature (ºC) - depth (Km) profile """
    
    plt.plot(T-273.15,Z*1e-3)
    ax = plt.gca()
    
    if ( display_title ):
        title = 'Temperature profile at steady-state ({})'
        plt.title(title.format(description))
        if title_in_plot:
            title = '{}'
            ax.text(0.97, 0.1, title.format(description),
                verticalalignment='center', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=16, bbox =dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.5) )#{'facecolor':'white', 'alpha':0.6, 'pad':5}), pad=5
            plt.title("")
    plt.xlabel("Temperature (º Celsius)")
    plt.ylabel("Depth (km)")
    if (xlim):
        plt.xlim(xlim)
    if (ylim):
        plt.ylim(ylim)
    
    ax.invert_yaxis()
            
    return ax

@ add_subplt()
def plot_temperature_pressure_profile(T,P,xlim=None, ylim=None,
                        display_title=True,title_in_plot=False,
                        printOut=False,description='unknown lithos.',*args, **kwargs):
    """ PLot temperature (ºC) - pressure (Pa) profile """
    
    plt.plot(T-273.15,P)
    ax = plt.gca()
    
    if ( display_title ):
        title = 'Temperature profile at steady-state ({})'
        plt.title(title.format(description))
        if title_in_plot:
            title = '{}'
            ax.text(0.97, 0.1, title.format(description),
                verticalalignment='center', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=16, bbox =dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.5) )#{'facecolor':'white', 'alpha':0.6, 'pad':5}), pad=5
            plt.title("")
    plt.xlabel("Temperature (º Celsius)")
    plt.ylabel("Depth (km)")
    if (xlim):
        plt.xlim(xlim)
    if (ylim):
        plt.ylim(ylim)
    
    ax.invert_yaxis()
    
    return ax
