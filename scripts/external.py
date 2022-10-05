#!/usr/bin/env python

import numpy as np
import colorsys
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator,ScalarFormatter)
import netCDF4
import os
import pyproj
import pickle
from statsmodels.stats.weightstats import DescrStatsW
from scipy.interpolate import RegularGridInterpolator

from matplotlib import rcParams, rcParamsDefault

import xarray as xr

geod = pyproj.Geod(ellps='WGS84')

def define_rcParams():
    rcParams.update({
    "text.usetex": False,
    "font.size": 6})
    rcParams['axes.titlesize'] = 8
    rcParams['axes.labelsize'] = 6
    rcParams['lines.linewidth'] = 1
    rcParams['lines.markersize'] = 4
    rcParams['xtick.labelsize'] = 6
    rcParams['ytick.labelsize'] = 6
    rcParams['figure.figsize'] = [10/2.54, 8/2.54]
    return rcParams

def gmtColormap(fileName,GMTPath = None):
    ''' https://scipy-cookbook.readthedocs.io/items/Matplotlib_Loading_a_colormap_dynamically.html
        modified from James Boyle and Andrew Straw - Thomas Theunissen 2021'''
    if type(GMTPath) == type(None):
        filePath = "./"+ fileName+".cpt"
    else:
        filePath = GMTPath+"/"+ fileName +".cpt"
    try:
        f = open(filePath)
    except:
        print("file "+filePath+"not found")
        return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if (len(l)>1): 
            if l[0] == "#":
               if ls[-1] == "HSV":
                   colorModel = "HSV"
                   continue
               else:
                   continue
            if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
               pass
            else:
                x.append(float(ls[0]))
                r.append(float(ls[1]))
                g.append(float(ls[2]))
                b.append(float(ls[3]))
                xtemp = float(ls[4])
                rtemp = float(ls[5])
                gtemp = float(ls[6])
                btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = np.array( x , np.float32)
    r = np.array( r , np.float32)
    g = np.array( g , np.float32)
    b = np.array( b , np.float32)
    if colorModel == "HSV":
       for i in range(r.shape[0]):
           rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
           r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
       for i in range(r.shape[0]):
           rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
           r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r/255.
        g = g/255.
        b = b/255.
    xNorm = (x - x[0])/(x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    return (colorDict)

def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    https://www.geophysique.be/2011/02/20/matplotlib-basemap-tutorial-09-drawing-circles/
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):

        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi	

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    ''' Plotting circles on a map '''
    glon1 = centerlon
    glat1 = centerlat
    X = [] ;
    Y = [] ;
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    X=np.asarray(X)
    Y=np.asarray(Y)
    diff = 999999999999999
    if (np.min(X)>0):
        diff = np.max(X)-np.min(X)
    elif (np.min(X)<0 and np.max(X)>0):
        diff = np.max(X)-np.min(X)
    elif (np.min(X)<0 and np.max(X)<0):
        diff = np.abs(np.min(X))-np.abs(np.max(X))
    # Wrapping the circle correctly across map bounds
    # simple fix enough here
    if (diff>300):
        X2 = X[X>0]
        Y2 = Y[X>0]
        Y  = Y[X<=0]
        X  = X[X<=0]
        X2,Y2 = m(X2,Y2)
        p = plt.plot(X2,Y2,**kwargs)
        X,Y = m(X,Y)
        if 'color' in kwargs:
            plt.plot(X,Y,**kwargs)
        else:
            plt.plot(X,Y,color=p[0].get_color(),**kwargs)
    else:
        #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
        X,Y = m(X,Y)
        plt.plot(X,Y,**kwargs)

def plot_histo(data,title,
               filename='histo.pdf',
               xlabel='Elevation (m)',unit='m',
               GaussianModel=False,sigmamodel=270,meanmodel=-2950,
               legends="upper right",text="left",
               approximation_display="int",
               xlim=None,
               savefig=False,
               fig_x=9.5,
               fig_y=9,
               weights=None,nbins=40):
    '''
       Plot histogram with some statistics
    '''
    define_rcParams()
    plt.figure(figsize=(fig_x/2.54,fig_y/2.54))

    if xlim:
        data = np.ma.MaskedArray(data, mask=( (data<xlim[0]) | (data>xlim[1]) ))

    n, bins, patches = plt.hist(x=data,bins=nbins, color='#0504aa',alpha=0.7, rwidth=0.85,weights=weights)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.title(title)

    mean   = np.nanmean(data)
    median = np.ma.median(data)
    sigma  = np.nanstd(data)

    if weights is not None:
        ma      = np.ma.MaskedArray(data, mask=np.isnan(data))
        meanW   = np.ma.average(ma, weights=weights)
        dsw     = DescrStatsW(ma, weights=weights)
        stdW    = dsw.std  # weighted std
        mean    = meanW
        sigma   = stdW

    xval  = np.linspace(np.nanmin(data),np.nanmax(data),1000)
    yval  = np.exp(-(xval-mean)**2/(2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    yval  = yval*n.max()/np.nanmax(yval)

    plt.plot(xval,yval,label='Data Gaussian model')

    if GaussianModel:
        yval2 = np.exp(-(xval-meanmodel)**2/(2*sigmamodel**2)) / np.sqrt(2*np.pi*sigmamodel**2)
        yval2 = yval2*n.max()/np.nanmax(yval2)
        p = plt.plot(xval,yval2,label='Estimated Gaussian model')

    if approximation_display=="int":
        mean  = int(np.ceil(mean))
        med   = int(np.ceil(median))
        sigma = int(np.ceil(sigma))
    else:
        mean  = np.around(mean,1)
        med   = np.around(median,1)
        sigma = np.around(sigma,1)

    if legends=="upper right":
        plt.legend(loc=legends, bbox_to_anchor=(1,1),framealpha=0.4)
    elif legends=="upper left":
        plt.legend(loc=legends, bbox_to_anchor=(0,1),framealpha=0.4)
    ax = plt.gca()
    if text=="right":
        xtext = 0.7
    else:
        xtext = 0.05
    ytext = 0.7 ; voff = 0.035 ; hoff = 0.1
    if weights is not None:
        plt.text(xtext,ytext       ,r'$\overline{elev}_W$',transform=ax.transAxes) ; plt.text(xtext+hoff,ytext       ,r'$=$'+str(mean)  +' '+unit,transform=ax.transAxes)
    else:
        plt.text(xtext,ytext       ,r'$\overline{elev}$',transform=ax.transAxes)   ; plt.text(xtext+hoff,ytext       ,r'$=$'+str(mean)  +' '+unit,transform=ax.transAxes)
    plt.text(xtext,ytext-voff  ,r'$median$',transform=ax.transAxes)          ; plt.text(xtext+hoff,ytext-voff  ,r'$=$'+str(med)   +' '+unit,transform=ax.transAxes)
    if weights is not None:
        plt.text(xtext,ytext-2*voff,r'$\sigma_{W}$',transform=ax.transAxes)          ; plt.text(xtext+hoff,ytext-2*voff,r'$=$'+str(sigma) +' '+unit,transform=ax.transAxes)
    else:
        plt.text(xtext,ytext-2*voff,r'$\sigma$',transform=ax.transAxes)          ; plt.text(xtext+hoff,ytext-2*voff,r'$=$'+str(sigma) +' '+unit,transform=ax.transAxes)
    
    if GaussianModel:
        # model
        if weights is None:
            vshift = 0.15
        else:
            vshift = 0.18
        plt.text(xtext,ytext-vshift     ,r'$\overline{elev}$',transform=ax.transAxes,color=p[0].get_color())
        plt.text(xtext+hoff,ytext-vshift,r'$=$'+str(meanmodel)  +' '+unit,transform=ax.transAxes,color=p[0].get_color())
        plt.text(xtext,ytext-vshift-voff     ,r'$\sigma$',transform=ax.transAxes,color=p[0].get_color())
        plt.text(xtext+hoff,ytext-vshift-voff,r'$=$'+str(sigmamodel) +' '+unit,transform=ax.transAxes,color=p[0].get_color())

    maxfreq = n.max()
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    if xlim:
        plt.xlim(xlim)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    if (savefig):
        plt.savefig(filename,dpi=300)

    return ax


def plot_correlation(x,y,title,filename='correlation.pdf',
                     nbins = 20,
                     xlabel='x',ylabel='y',unit='m',text="right",plottext=True,
                     xlim=None,ylim=None,ticks=None,
                     savefig=False,
                     fig_x=9.5,
                     fig_y=9,
                     weights=None):
    ''' ticks = [tick_dx1,tick_dx2,tick_dy1,tick_dy2]
        unit = from "y" variable
    '''
    define_rcParams()
   
    plt.figure(figsize=(fig_x/2.54,fig_y/2.54))
    ax = plt.gca()
    
    plt.plot(x,y,'ko',markersize=0.5,zorder=0,rasterized=True)

    if (xlim):
        xstat = np.linspace(xlim[0],xlim[1],nbins)
        nb12 = (xlim[1]-xlim[0])/(2*nbins)
    else:
        xstat = np.linspace(np.nanmin(x),np.nanmax(x),nbins)
        nb12 = (np.nanmax(x)-np.nanmin(x))/(2*nbins)

    mean = [] ; median = [] ; sigma = [] ; xused = []
    for rate in xstat:
        data  = y[( (x>=rate-nb12) & (x<=rate+nb12))]
        if weights is not None:
            selweights = weights[( (x>=rate-nb12) & (x<=rate+nb12))]
        if data!=[]:
            xused.append(rate)
            med    = np.ma.median(data)
            if weights is None:
                avg    = np.nanmean(data)
                std    = np.nanstd(data)
            else:
                ma      = np.ma.MaskedArray(data, mask=np.isnan(data))
                avgW    = np.ma.average(ma, weights=selweights)
                dsw     = DescrStatsW(ma, weights=selweights)
                stdW    = dsw.std  # weighted std
                avg     = avgW
                std     = stdW
            mean.append(avg)
            median.append(med)
            sigma.append(std)
    mean = np.asarray(mean) ; sigma = np.asarray(sigma) ; median = np.asarray(median) ; xused = np.asarray(xused)

    plt.plot(xused,median,color='r',zorder=3,linewidth=2,label='median')
    plt.plot(xused,np.add(median,sigma),color='g',linewidth=2,zorder=4)
    plt.plot(xused,np.subtract(median,sigma),color='g',linewidth=2,zorder=5)

    plt.plot(xused,mean,color='b',zorder=3,linewidth=2,label='mean')

    if plottext:
        if (xlim):
            xstat = np.linspace(xlim[0],xlim[1],1000)
        else:
            xstat = np.linspace(np.nanmin(x),np.nanmax(x),1000)
        themedian = np.zeros_like(xstat)
        themedian = themedian + np.ma.median(median)
        mean  = np.around(np.nanmean(median),1)
        med   = np.around(np.ma.median(median),1)
        sigma = np.around(np.nanstd(median),1)
        plt.plot(xstat,themedian)

        if text=="left":
            xtext = 0.25
        elif text=="right":
            xtext = 0.7
        else:
            xtext = 0.35
        ytext = 0.85 ; voff = 0.03 ; hoff = 0.1
        plt.text(xtext,ytext       ,r'$\overline{elev}$',transform=ax.transAxes) ; plt.text(xtext+hoff,ytext       ,r'$=$'+str(mean)  +' '+unit,transform=ax.transAxes)
        plt.text(xtext,ytext-voff  ,r'$median$',transform=ax.transAxes)          ; plt.text(xtext+hoff,ytext-voff  ,r'$=$'+str(med)   +' '+unit,transform=ax.transAxes)
        plt.text(xtext,ytext-2*voff,r'$\sigma$',transform=ax.transAxes)          ; plt.text(xtext+hoff,ytext-2*voff,r'$=$'+str(sigma) +' '+unit,transform=ax.transAxes)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    if not ticks:
        if (xlim):
            tick_dx1 = int(np.ceil((xlim[1]-xlim[0])/6))
            tick_dx2 = int(np.ceil((xlim[1]-xlim[0])/24))
        else:
            tick_dx1 = int(np.ceil((np.nanmax(x)-np.nanmin(x))/6))
            tick_dx2 = int(np.ceil((np.nanmax(x)-np.nanmin(x))/24))
        if (ylim):
            tick_dy1 = int(np.ceil((ylim[1]-ylim[0])/5))
            tick_dy2 = int(np.ceil((ylim[1]-ylim[0])/20))
        else:
            tick_dy1 = int(np.ceil((np.nanmax(y)-np.nanmin(y))/5))
            tick_dy2 = int(np.ceil((np.nanmax(y)-np.nanmin(y))/20))
    else:
        tick_dx1 = ticks[0]
        tick_dx2 = ticks[1]
        tick_dy1 = ticks[2]
        tick_dy2 = ticks[3]

    ax.yaxis.set_major_locator(MultipleLocator(tick_dy1))
    ax.yaxis.set_minor_locator(MultipleLocator(tick_dy2))
    ax.xaxis.set_major_locator(MultipleLocator(tick_dx1))
    ax.xaxis.set_minor_locator(MultipleLocator(tick_dx2))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if (xlim):
        plt.xlim(xlim)
    if (ylim):
        plt.ylim(ylim)

    plt.legend()

    if (savefig):
        plt.savefig(filename,dpi=300)
    return ax

#--------------------------- READING DATA ---------------------

def load_spreading_rate(path='./data/topography/'):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Half Spreading rate
    nc_spreading_rate = netCDF4.Dataset(path+'rate.3.2.nc')
    #print(nc_spreading_rate.variables.keys())
    x=nc_spreading_rate.variables['lon'][:]
    y=nc_spreading_rate.variables['lat'][:]
    spreading_rate=nc_spreading_rate.variables['z'][:]
    spreading_rate=2*spreading_rate/1000
    print("Full seafloor spreading rate (Muller et al, 2008) cm/yr min/max {} {}".format(np.nanmin(spreading_rate),np.nanmax(spreading_rate)))
    return x,y,spreading_rate

def load_seafloor_ages(path='./data/topography/',sampling=None):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Seafloor age
    nc_age = netCDF4.Dataset(path+'age.3.2.nc')
    #print(nc_age.variables.keys())
    x=nc_age.variables['lon'][:]
    y=nc_age.variables['lat'][:]
    age=nc_age.variables['z'][:]
    age=age/100
    print("Seafloor age (Muller et al, 2008) Myrs min/max {} {}".format(np.nanmin(age),np.nanmax(age)))
    return x,y,age

def load_strain_rate(path='./data/topography/',sampling=None):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Strain rate
    nc_strain_rate = netCDF4.Dataset(path+'GSRM_strain_2m.grd')
    #print(nc_strain_rate.variables.keys())
    x=nc_strain_rate.variables['lon'][:]
    y=nc_strain_rate.variables['lat'][:]
    strain_rate=nc_strain_rate.variables['z'][:]
    strain_rate=strain_rate*1e-9/(365*24*3600)
    strain_rate=np.where(np.isnan(strain_rate),1e-21,strain_rate)
    print("Strain rate (Kreemer et al, 2014) s-1 min/max {} {}".format(np.nanmin(strain_rate),np.nanmax(strain_rate)))
    if sampling is not None:
        def values(x,y,x_old,y_old,values_old):
            X,Y = np.meshgrid(x,y)
            f_interp = RegularGridInterpolator((y_old, x_old), values_old)
            values = np.zeros_like(X)
            nlon = values.shape[0]
            nlat = values.shape[1]
            minx=np.min(x_old)
            maxx=np.max(x_old)
            miny=np.min(y_old)
            maxy=np.max(y_old)
            for i in np.arange(0,nlon):
                for j in np.arange(0,nlat):
                    if ((X[i,j]<=minx) | (X[i,j]>=maxx) | (Y[i,j]<=miny) | (Y[i,j]>=maxy)):
                        ilon   = np.argmin(np.abs(x_old - x[i]))
                        ilat   = np.argmin(np.abs(y_old - y[j]))
                        values[i,j] = strain_rate[ilat,ilon]
                    else:
                        values[i,j] = f_interp([Y[i,j],X[i,j]])
            return values
        nlon = int(np.floor(360/sampling))
        nlat = int(np.floor(180/sampling))
        xnew = np.linspace(-180+sampling/2,180-sampling/2,nlon)
        ynew = np.linspace(-90+sampling/2,90-sampling/2,nlat)
        strain_rate = values(xnew,ynew,x,y,strain_rate)
    return x,y,strain_rate

def load_etopo(path='./data/topography/',filtered=True,resolution=2,corrected_from_ice=False):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Etopo 1
    collection = False
    if (filtered==True):
        if corrected_from_ice:
            filename  = 'ETOPO1_BedIceCorrected_g_gmt4_filtered.grd'
        else:
            filename  = 'ETOPO1_Bed_g_gmt4_filtered.grd'
        if (os.path.isfile(path+filename)):
            nc_etopo1 = netCDF4.Dataset(path+filename)
        else:
            name      = filename.split('.grd')[0]
            nc_etopo1 = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                          concat_dim=['lon'], combine='nested',engine='netcdf4')
            collection = True
    else:
        if resolution==1:
            filename = 'ETOPO1_Bed_g_gmt4.grd'
        else:
            filename = 'ETOPO1_Bed_g_gmt4_2m.grd'
        if ( not os.path.isfile(path+filename)):
            if filename == 'ETOPO1_Bed_g_gmt4.grd':
                print("Unfiltered raw 1 arc-min ETOPO1 dataset not available (check if you have downloaded FigShare Dataset) - Not available on Binder ({})".format(path+filename))
                x=[] ; y=[] ; elev=[]
                return x,y,elev
            name      = filename.split('.grd')[0]
            nc_etopo1 = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                          concat_dim=['lon'], combine='nested',engine='netcdf4')
            collection = True
        else:
            nc_etopo1 = netCDF4.Dataset(path+filename)
    #print(nc_etopo1.variables.keys())
    x=nc_etopo1.variables['lon'][:]
    y=nc_etopo1.variables['lat'][:]
    elev=nc_etopo1.variables['z'][:]
    if collection:
        x=x.to_numpy()
        y=y.to_numpy()
        elev=elev.to_numpy()
    print("ETOPO 1 m ({}) min/max {} {}".format(filename,np.nanmin(elev),np.nanmax(elev)))
    return x,y,elev

def get_shape_etopo(path='./data/topography/'):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Etopo 1
    filename  = 'ETOPO1_Bed_g_gmt4_filtered.grd'
    if (os.path.isfile(path+filename)):
        nc_etopo1 = netCDF4.Dataset(path+filename)
    else:
        name      = filename.split('.grd')[0]
        nc_etopo1 = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['lon'], combine='nested',engine='netcdf4')
    return nc_etopo1

def load_hotspots(path='./data/topography/',write_id=False,write_grid=False,sampling=None):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Hot spots
    filename = 'dist_closest_hs.grd'
    if write_grid:
        if ( not os.path.isfile(path+filename)):
            # Very slow to compute this grid
            filename2='Morgan-Morgan_2007_hs.txt'
            data = np.loadtxt(path+filename2,dtype='float',delimiter=' ')
            lat_read = data[:,1]
            lon_read = data[:,0]
            nc_etopo1 = get_shape_etopo()
            x=nc_etopo1.variables['lon'][:]
            y=nc_etopo1.variables['lat'][:]
            dist_closest_hs = np.zeros([nc_etopo1.variables['lat'].shape[0], nc_etopo1.variables['lon'].shape[0]])
            id_closest_hs   = np.zeros([nc_etopo1.variables['lat'].shape[0], nc_etopo1.variables['lon'].shape[0]])
            ipt = 0
            npt = len(x)*len(y)
            for ilon in np.arange(0,len(x)):
                lon1 = x[ilon]
                for ilat in np.arange(0,len(y)):
                    min_dist = 9e9
                    min_id   = 0
                    lat1     = y[ilat]
                    ipt = ipt + 1
                    for ihs in np.arange(0,len(lon_read)):
                        lon2 = lon_read[ihs]
                        lat2 = lat_read[ihs]
                        #if (np.abs(lon1-lon2)<25 and np.abs(lat1-lat2)<25):
                        azimuth1, azimuth2, dist = geod.inv(lon1, lat1, lon2, lat2)
                        if (dist<min_dist):
                            min_dist = dist
                            min_id   = data[ihs,3]
                        #print('Hot spot {}'.format(data[ihs,3]))
                    dist_closest_hs[ilat,ilon] = min_dist
                    id_closest_hs[ilat,ilon]   = min_id
                    if ipt%10000 == 0:
                        print("{}/{}".format(ipt,npt))
            # Create a new netcdf file with distances to closest hotspots
            if write_id:
                with netCDF4.Dataset(path+'id_closest_hs.grd', "w", format="NETCDF4_CLASSIC") as f:
                    f.description = 'ID of the closest hotspot'
                    # dimensions
                    f.createDimension('x', nc_etopo1.variables['lon'].shape[0])
                    f.createDimension('y', nc_etopo1.variables['lat'].shape[0])
                    # variables
                    xnew = f.createVariable('x', 'f4', ('x',))
                    ynew = f.createVariable('y', 'f4', ('y',))
                    znew = f.createVariable('z', 'f4', ('y', 'x',))
                    xnew.units = "degrees east"
                    ynew.units = "degrees north"
                    znew.units = "-"
                    # data
                    xnew[:] = x
                    ynew[:] = y
                    znew[:,:] = id_closest_hs
    
            with netCDF4.Dataset(path+filename, "w", format="NETCDF4_CLASSIC") as f:
                f.description = 'Distance to the closest hotspot'
                # dimensions
                f.createDimension('x', nc_etopo1.variables['lon'].shape[0])
                f.createDimension('y', nc_etopo1.variables['lat'].shape[0])
                # variables
                xnew = f.createVariable('x', 'f4', ('x',))
                ynew = f.createVariable('y', 'f4', ('y',))
                znew = f.createVariable('z', 'f4', ('y', 'x',))
                xnew.units = "degrees east"
                ynew.units = "degrees north"
                znew.units = "m"
                # data
                xnew[:] = x
                ynew[:] = y
                znew[:,:] = dist_closest_hs
    
    if (os.path.isfile(path+filename)):
        nc_dist_closest_hs = netCDF4.Dataset(path+filename)
    else:
        name               = filename.split('.grd')[0]
        nc_dist_closest_hs = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['x'], combine='nested',engine='netcdf4')
    #print(nc_dist_closest_hs.variables.keys())
    x=nc_dist_closest_hs.variables['x'][:]
    y=nc_dist_closest_hs.variables['y'][:]
    dist_closest_hs=nc_dist_closest_hs.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        dist_closest_hs=dist_closest_hs.to_numpy()
    print("Distances closest hot spot (Morgan and Morgan, 2007) m min/max {} {}".format(np.nanmin(dist_closest_hs),np.nanmax(dist_closest_hs)))
    
    return x,y,dist_closest_hs

def load_crustal_thickness(path='./data/topography/'):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # Crustal thickness
    filename = 'crust_thickness_GEMMA_2m.grd'
    if (os.path.isfile(path+filename)):
        nc_crust = netCDF4.Dataset(path+filename)
    else:
        name      = filename.split('.grd')[0]
        nc_crust = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['lon'], combine='nested',engine='netcdf4')
    #print(nc_crust.variables.keys())
    x=nc_crust.variables['lon'][:]
    y=nc_crust.variables['lat'][:]
    crustal_thickness=nc_crust.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        crustal_thickness=crustal_thickness.to_numpy()
    print("Crustal thickness (GEMMA 2 arc-min) km min/max {} {}".format(np.nanmin(crustal_thickness),np.nanmax(crustal_thickness)))
    return x,y,crustal_thickness

def load_lithospheric_thickness(path='./data/topography/'):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # lithospheric thickness
    filename = 'lith_ave_no_slabs_-180_180_2m.grd'
    if (os.path.isfile(path+filename)):
        nc_liththick = netCDF4.Dataset(path+filename)
    else:
        name         = filename.split('.grd')[0]
        nc_liththick = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['lon'], combine='nested',engine='netcdf4')
    #print(nc_liththick.variables.keys())
    x=nc_liththick.variables['lon'][:]
    y=nc_liththick.variables['lat'][:]
    lithospheric_thickness=nc_liththick.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        lithospheric_thickness=lithospheric_thickness.to_numpy()
    print("Lithospheric thickness (SteinBerger 2016) km min/max {} {}".format(np.nanmin(lithospheric_thickness),np.nanmax(lithospheric_thickness)))
    lithospheric_thickness=np.where(np.isnan(lithospheric_thickness),-999,lithospheric_thickness)
    return x,y,lithospheric_thickness

def load_age_lithos(path='./data/topography/'):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    # thermal age continent
    filename = 'mant_age_map_-180_180_qgis_2m.grd'
    if (os.path.isfile(path+filename)):
        nc_agelith = netCDF4.Dataset(path+filename)
    else:
        name       = filename.split('.grd')[0]
        nc_agelith = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['lon'], combine='nested',engine='netcdf4')
    #print(nc_agelith.variables.keys())
    x=nc_agelith.variables['lon'][:]
    y=nc_agelith.variables['lat'][:]
    age_lithos=nc_agelith.variables['z'][:]*1e3
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        age_lithos=age_lithos.to_numpy()
    print("Age lithosphere (Poupinet_Shapiro_2008) Ma min/max {} {}".format(np.nanmin(age_lithos),np.nanmax(age_lithos)))
    age_lithos=np.where(np.isnan(age_lithos),-999,age_lithos)
    return x,y,age_lithos

def define_MOR_pts(path='./data/topography/',selection_name='MOR_pts_all',distance_between_pts_along_ridges = 25):
    '''
        Define points with elevation and spreading rates along MOR

        selection_name =  'MOR_pts_all' or 'MOR_pts_far_from_hs' or 'MOR_pts_close_to_hs'
    '''
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'

    filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'

    distance_between_pts_along_ridges = distance_between_pts_along_ridges*1e3

    if not (os.path.isfile(path+filename)):
        print("#################### Loading datasets #########################")
        x,y,elev                   = load_etopo(path=path)
        x,y,spreading_rate         = load_spreading_rate(path=path)
        x,y,age                    = load_seafloor_ages(path=path)
        x,y,strain_rate            = load_strain_rate(path=path)
        x,y,dist_closest_hs        = load_hotspots(path=path)

        print("#################### Applying mask on datasets #########################")
        # data selection
        min_dist_hs                     = 1000000 # m
        max_seafloor_age                = 10      # Myrs the width depends on spreading rate - 10 Myrs is a good compromise for computational reason 
                                                  #                                            It gives 50 km from ridge axis for ultra-slow 1 cm/yr full spreading rate MOR
        max_seafloor_age_for_ridge_axis = 0.5
        threshold_strain_rate           = 1e-16   # s-1
        
        xx, yy = np.meshgrid(x, y)
      
        if selection_name=='MOR_pts_all':
            mask      = ( (elev<0) & (age<=max_seafloor_age) )
            mask_axis = ( (elev<0) & (age<=max_seafloor_age_for_ridge_axis) )
        elif selection_name=='MOR_pts_far_from_hs':
            mask      = ( (elev<0) & (age<=max_seafloor_age) & (dist_closest_hs > min_dist_hs) )
            mask_axis = ( (elev<0) & (age<=max_seafloor_age_for_ridge_axis) & (dist_closest_hs > min_dist_hs) )
        elif selection_name=='MOR_pts_close_to_hs':
            mask      = ( (elev<0) & (age<=max_seafloor_age) & (dist_closest_hs <= min_dist_hs) )
            mask_axis = ( (elev<0) & (age<=max_seafloor_age_for_ridge_axis) & (dist_closest_hs <= min_dist_hs) )
        else:
            print("ERROR incorrect selection_name ")
            quit()
    
        # this array is used to define the localisation of the MOR
        active_MOR_elev           = elev[mask]
        active_MOR_x              = xx[mask]
        active_MOR_y              = yy[mask]
        active_MOR_x_axis         = xx[mask_axis]
        active_MOR_y_axis         = yy[mask_axis]
        active_MOR_spreading_rate = spreading_rate[mask]
        active_MOR_strain_rate    = strain_rate[mask]
   
        dd                                = 1.5   # Distance to look for points in the grid that belong to the same MOR segment 
                                                  # given in degrees for computational reason, could be function of the spreading rate
                                                  # Here, we define a constant that cover all cases (Fig. S4b) (~150-200 km)
                                                  # W    ~ 100-150 km ~ distance between rift flanks at ultra-slow spreading rates (Fig. 3)
                                                  #
    
        new_active_MOR_y = [] ; new_active_MOR_x = [] ; new_active_MOR_elev = [] ; new_active_MOR_spreading_rate = []
        ipt = 0

        print('Total #pts on the grid for age<={} Myrs = {} '.format(max_seafloor_age_for_ridge_axis,len(active_MOR_x_axis)))

        print("#################### Browsing all MOR points #########################")
        for xpt,ypt in zip(active_MOR_x_axis,active_MOR_y_axis):
            xsel = active_MOR_x_axis[ ( (np.abs(active_MOR_x_axis-xpt)<=dd/2) & (np.abs(active_MOR_y_axis-ypt)<=dd/2) ) ]
            ysel = active_MOR_y_axis[ ( (np.abs(active_MOR_x_axis-xpt)<=dd/2) & (np.abs(active_MOR_y_axis-ypt)<=dd/2) ) ]
            newx = np.median(xsel)
            newy = np.median(ysel)
            if (ipt==0):
                new_active_MOR_x.append(newx)
                new_active_MOR_y.append(newy)
                esel = active_MOR_elev[ ( (np.abs(active_MOR_x-newx)<=dd/2) & (np.abs(active_MOR_y-newy)<=dd/2) ) ]
                new_active_MOR_elev.append(np.max(esel))
                srsel = active_MOR_spreading_rate[ ( (np.abs(active_MOR_x-newx)<=dd/2) & (np.abs(active_MOR_y-newy)<=dd/2) ) ]
                new_active_MOR_spreading_rate.append(np.median(srsel))
            else:
                stsel = active_MOR_strain_rate[ ( (np.abs(active_MOR_x-newx)<=dd/2) & (np.abs(active_MOR_y-newy)<=dd/2) ) ]
                if (np.any(stsel>=threshold_strain_rate)):
                    azimuth1, azimuth2, dist = geod.inv(new_active_MOR_x[-1], new_active_MOR_y[-1], newx, newy)
                    if ( dist >= distance_between_pts_along_ridges ):
                        esel = active_MOR_elev[ ( (np.abs(active_MOR_x-newx)<=dd/2) & (np.abs(active_MOR_y-newy)<=dd/2) ) ]
                        srsel = active_MOR_spreading_rate[ ( (np.abs(active_MOR_x-newx)<=dd/2) & (np.abs(active_MOR_y-newy)<=dd/2) ) ]
                        new_active_MOR_x.append(newx)
                        new_active_MOR_y.append(newy)
                        new_active_MOR_elev.append(np.max(esel))
                        new_active_MOR_spreading_rate.append(np.median(srsel))
            ipt = ipt + 1
            if ipt%5000 == 0:
                print("{}/{}".format(ipt,len(active_MOR_x_axis)))

        new_active_MOR_x              = np.asarray(new_active_MOR_x)
        new_active_MOR_y              = np.asarray(new_active_MOR_y)
        new_active_MOR_elev           = np.asarray(new_active_MOR_elev)
        new_active_MOR_spreading_rate = np.asarray(new_active_MOR_spreading_rate)
        with open(path+filename, 'wb') as filehandle:
            pickle.dump([new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate],filehandle)

        print('Total defined pts along MOR = {} '.format(len(new_active_MOR_x)))

    else:
        print("This selection already exists ({})".format(path+filename))

        with open(path+filename, 'rb') as filehandle:
            # read the data as binary data stream
            [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)

        print('Total defined pts along MOR = {} '.format(len(new_active_MOR_x)))

