#!/usr/bin/env python

import numpy as np
import netCDF4
from matplotlib import rcParams, rcParamsDefault
import matplotlib.pylab as plt
import matplotlib as mpl
from scipy.interpolate import RegularGridInterpolator
import pyproj
from pygc import great_circle
#from mpl_toolkits.basemap import Basemap
import os

from .external import *
from .wintercg import *

rcParams.update({
    "text.usetex": False,
    "font.size": 6})
rcParams['axes.titlesize'] = 8
rcParams['axes.labelsize'] = 6
rcParams['lines.linewidth'] = 1
rcParams['lines.markersize'] = 4
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6

geod            = pyproj.Geod(ellps='WGS84')

def plot_MOR_profile(path='./data/topography/',name='SEIR_91E',use_fixed_width=True,width=300000,xlim=[0,300],ylim=[-4500,-1500],plot_wintercg_data=False):
    '''
       Plot MOR profiles (Some of them on Fig. 3).
       name= 'all' (to display all of them)
             'MAR_05N', 'MAR_12S', 'MAR_22N', 'MAR_27N', 'MAR_30S', 'MAR_46S', 'MAR_57N', 'EPR_10N', 'EPR_15S', 'EPR_40S', 'EPR_65S', 'SEIR_110E', 'SEIR_91E'
             'CIR_2N', 'CIR_11S', 'SWIR_24E', 'ARC_24E'
    '''

    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'
    
    allfilename = [] ; allsection_pts = [] ; alloutputname = []
    
    # INPUT DATA ----------------------------------------------------------------------
    if (name=='MAR_05N' or name=='all'):
        allfilename.append('GMRT3.6/MAR_05N_topo-mask.grd')
        allsection_pts.append([[-33.22,4.49],[-31.93,4.49]])
        alloutputname.append('MAR_05N.pdf')
    if (name=='MAR_12S' or name=='all'):
        allfilename.append('GMRT3.6/MAR_12S_topo-mask.grd')
        allsection_pts.append([[-15.254,-12.485],[-14.037,-12.496]])
        alloutputname.append('MAR_12S.pdf')
    if (name=='MAR_22N' or name=='all'):
        allfilename.append('GMRT3.6/MAR_22N_topo-mask.grd')
        allsection_pts.append([[-45.86,22.43],[-44.4,22.145]])
        alloutputname.append('MAR_22N.pdf')
    if (name=='MAR_27N' or name=='all'):
        allfilename.append('GMRT3.6/MAR_27N_topo-mask.grd')
        allsection_pts.append([[-45.32584584584585,26.476816816816818],[-44.183943943943945,26.010660660660662]])
        alloutputname.append('MAR_27N.pdf')
    if (name=='MAR_30S' or name=='all'):
        allfilename.append('GMRT3.6/MAR_30S_topo-mask.grd')
        allsection_pts.append([[-15.062,-30.232],[-12.375,-29.922]])
        alloutputname.append('MAR_30S.pdf')
    if (name=='MAR_46S' or name=='all'):
        allfilename.append('GMRT3.6/MAR_46S_topo-mask.grd')
        allsection_pts.append([[-15.538,-47.221],[-11.665675675675676,-46.00085585585585]])
        alloutputname.append('MAR_46S.pdf')
    if (name=='MAR_57N' or name=='all'):
        allfilename.append('GMRT3.6/MAR_57N_topo-mask.grd')
        allsection_pts.append([[-35.99,56.95],[-31.90,56.76]])
        alloutputname.append('MAR_57N.pdf')
    if (name=='EPR_10N' or name=='all'):
        allfilename.append('GMRT3.6/EPR_10N_topo-mask.grd')
        allsection_pts.append([[-104.608,9.31],[-103.75,9.34]])
        alloutputname.append('EPR_10N.pdf')
    if (name=='EPR_15S' or name=='all'):
        allfilename.append('GMRT3.6/EPR_15S_topo-mask.grd')
        allsection_pts.append([[-114.868,-16.020],[-110.953,-16.870]])
        alloutputname.append('EPR_15S.pdf')
    if (name=='EPR_40S' or name=='all'):
        allfilename.append('GMRT3.6/EPR_40S_topo-mask.grd')
        allsection_pts.append([[-91.93,-39.26],[-91.11,-39.28]])
        alloutputname.append('EPR_40S.pdf')
    if (name=='EPR_65S' or name=='all'):
        allfilename.append('GMRT3.6/EPR_65S_topo-mask.grd')
        allsection_pts.append([[-173.26,-63.36],[-167.14,-65.80]])
        alloutputname.append('EPR_65S.pdf')
    if (name=='SEIR_110E' or name=='all'):
        allfilename.append('GMRT3.6/SEIR_110E_topo-mask.grd')
        allsection_pts.append([[109.496,-50.036],[109.827,-49.169]])
        alloutputname.append('SEIR_110E.pdf')
    if (name=='SEIR_91E' or name=='all'):
        allfilename.append('GMRT3.6/SEIR_91E_topo-mask.grd')
        allsection_pts.append([[91.030,-43.560],[91.492,-42.891]])
        alloutputname.append('SEIR_91E.pdf')
    if (name=='CIR_2N' or name=='all'):
        allfilename.append('GMRT3.6/SEIR_91E_topo-mask.grd') #CIR_2N_topo.grd') no data
        allsection_pts.append([[65.863,1.8949],[66.752994,3.168003]])
        alloutputname.append('CIR_2N.pdf')
    if (name=='CIR_11S' or name=='all'):
        allfilename.append('GMRT3.6/SEIR_91E_topo-mask.grd') #CIR_2N_topo.grd') no data
        allsection_pts.append([[65.526129,-11.154338],[67.229010,-9.890581]])
        alloutputname.append('CIR_11S.pdf')
    if (name=='SWIR_24E' or name=='all'):
        allfilename.append('GMRT3.6/SWIR_24E_topo-mask.grd')
        allsection_pts.append([[23.678,-53.847],[24.305,-52.933]])
        alloutputname.append('SWIR_24E.pdf')
    if (name=='ARC_24E' or name=='all'):
        allfilename.append('none')
        allsection_pts.append([[23.678,-53.847],[24.305,-52.933]])
        alloutputname.append('ARC_24E.pdf')


    x,y,elev                = load_etopo(path=path)
    x_raw,y_raw,elev_raw    = load_etopo(path=path,filtered=False,resolution=1)
    x,y,spreading_rate      = load_spreading_rate(path=path)
    x,y,age                 = load_seafloor_ages(path=path)

    f_interp_age            = RegularGridInterpolator((y, x), age)
    f_interp_spreading_rate = RegularGridInterpolator((y, x), spreading_rate)
    f_interp_etopo1         = RegularGridInterpolator((y, x), elev)
    if (elev_raw!=[]):
        f_interp_etopo1_raw = RegularGridInterpolator((y_raw, x_raw), elev_raw)

    if plot_wintercg_data:
        maxdepth = -200 # max is -401
        alldepth = np.flip(np.arange(maxdepth,6,2))
        ylim     = [-1*ylim[1]/1e3,-1*maxdepth]
        #Load data
        temp    = []
        rhof    = []
        for depth in alldepth:
            x,y,T = load_wintercg_temperature(z=depth,verbose=False)
            f = RegularGridInterpolator((y, x), T.reshape(len(y),len(x)))
            temp.append(f)
            x,y,rho = load_wintercg_density(z=depth,verbose=False)
            f = RegularGridInterpolator((y, x), rho.reshape(len(y),len(x)))
            rhof.append(f)
        #xpt,ypt = map(*np.meshgrid(x,y))

    for filename,section_pts,outputname in zip(allfilename,allsection_pts,alloutputname):
    
        print("-----------------------")
        print(outputname.split('.')[0])
        print("-----------------------")
    
        # High-resolution regional dataset
        #print(filename)
        if ( not os.path.isfile(path+filename)):
            print("High-resolution regional dataset not available (check if you have downloaded FigShare Dataset) - Not available on Binder ({})".format(path+filename))
        else:
            nc = netCDF4.Dataset(path+filename)
            #print(nc.variables.keys())
    
            x_range=nc.variables['x_range'][:]
            y_range=nc.variables['y_range'][:]
            z_range=nc.variables['z_range'][:]
            spacing=nc.variables['spacing'][:]
            dimension=nc.variables['dimension'][:]
    
            z=nc.variables['z'][:]
            print("Elevation (GMRT3.6) m min/max {} {}".format(np.nanmin(z),np.nanmax(z)))
            
            z = np.flipud(z.reshape(dimension[1],dimension[0]))
            
            x        = np.linspace(np.nanmin(x_range),np.nanmax(x_range),dimension[0])
            y        = np.linspace(np.nanmin(y_range),np.nanmax(y_range),dimension[1])
            f_interp = RegularGridInterpolator((y, x),z)
            #print("Lon (GMRT3.6) degree min/max {} {}".format(np.nanmin(x),np.nanmax(x)))
            #print("Lat (GMRT3.6) degree min/max {} {}".format(np.nanmin(y),np.nanmax(y)))
            
        
        ################## PLOTTING ######################################
        
    #    plt.figure(figsize=(10,10))
    #    plt.imshow(z,origin='lower')
    #    plt.colorbar()
        
    #    plt.figure(figsize=(9/2.54,9/2.54))
    #    lat_ts = (np.ceil(np.nanmax(y_range))-np.floor(np.nanmin(y_range)))/2
    #    map = Basemap(projection='merc',llcrnrlat=np.nanmin(y_range),urcrnrlat=np.nanmax(y_range),
    #                                    llcrnrlon=np.nanmin(x_range),urcrnrlon=np.nanmax(x_range),lat_ts=lat_ts)
    #    map.drawcoastlines()
    #    map.drawparallels(np.arange(int(np.floor(np.nanmin(y_range))),int(np.ceil(np.nanmax(y_range))),1),labels=[1,0,0,0])
    #    map.drawmeridians(np.arange(int(np.floor(np.nanmin(x_range))),int(np.ceil(np.nanmax(x_range))),1),labels=[0,0,0,1])
    #    xx, yy = np.meshgrid(x, y)
    #    map.pcolormesh(xx, yy, z,latlon=True, cmap='RdBu_r')
        
        lon1 = section_pts[0][0]
        lat1 = section_pts[0][1]
        lon2 = section_pts[1][0]
        lat2 = section_pts[1][1]
        if use_fixed_width:
            azimuth1, azimuth2, distance = geod.inv(lon1, lat1, lon2, lat2)
            loncenter = (lon1+lon2)/2
            latcenter = (lat1+lat2)/2
            newpt2 = great_circle(distance=width/2, azimuth=azimuth1, latitude=latcenter, longitude=loncenter)
            newpt1 = great_circle(distance=width/2, azimuth=azimuth2, latitude=latcenter, longitude=loncenter)
            lon1   = newpt1['longitude'] ; lat1   = newpt1['latitude']
            lon2   = newpt2['longitude'] ; lat2   = newpt2['latitude']
        #print(lon1,lon2)
        #print(lat1,lat2)
        sectionX = np.linspace(lon1,lon2,1000)
        sectionY = np.linspace(lat1,lat2,1000)
        dist     = np.zeros_like(sectionX)
        elev     = np.zeros_like(sectionX)
        elev_etopo1 = np.zeros_like(sectionX)
        elev_etopo1_raw = np.zeros_like(sectionX)
       
        ipt      = 0
        for lon,lat in zip(sectionX,sectionY):
            if ( ( lat >= np.nanmin(y) and lat <= np.nanmax(y) and lon >= np.nanmin(x) and lon <= np.nanmax(x) ) and ( os.path.isfile(path+filename)) ):
                elev[ipt]                 = f_interp([lat,lon])[0]
            else:
                elev[ipt]                 = np.nan
            elev_etopo1[ipt]              = f_interp_etopo1([lat,lon])[0]
            if elev_raw==[]:
                elev_etopo1_raw[ipt]      = np.nan
            else:
                elev_etopo1_raw[ipt]      = f_interp_etopo1_raw([lat,lon])[0]
            azimuth1, azimuth2, dist[ipt] = geod.inv(lon1, lat1, lon, lat)
            ipt = ipt + 1
 
        if plot_wintercg_data:
            nz          = len(alldepth)
            nx          = len(dist/1e3)
            XX, ZZ      = np.meshgrid(dist/1e3,-1*alldepth)
            temperature = np.zeros_like(XX)
            density     = np.zeros_like(XX)
            ipt         = 0
            for lon,lat in zip(sectionX,sectionY):
                for k in np.arange(0,len(alldepth)):
                    temperature[k,ipt] = temp[k]([lat,lon])[0]
                    density[k,ipt]     = rhof[k]([lat,lon])[0]
                ipt = ipt + 1
        
        if not plot_wintercg_data:
            plt.figure(figsize=(10/2.54,4.5/2.54))
            plt.plot(dist/1e3,elev,color='black',linewidth=0.25)
            plt.plot(dist/1e3,elev_etopo1,color='r')
            plt.plot(dist/1e3,elev_etopo1_raw,color='b',linewidth=0.75)
        else:
            plt.figure(figsize=(20/2.54,15/2.54))
            ax = plt.gca()
            orientationbar = 'vertical'
            if orientationbar=='horizontal':
                shrink=0.5
                aspect=25
            else:
                shrink=1
                aspect=20
            #pcm = plt.pcolormesh(XX,ZZ,temperature,cmap=plt.cm.jet,vmin=0,vmax=1400,rasterized=True,shading='nearest')
            #plt.colorbar(pcm,extend='both',orientation=orientationbar,shrink=shrink,aspect=aspect)
            pcm = plt.pcolormesh(XX,ZZ,density,cmap=plt.cm.jet,vmin=3250,vmax=3400,rasterized=True,shading='nearest')
            cbar = plt.colorbar(pcm,extend='both',orientation=orientationbar,shrink=shrink,aspect=aspect)
            cbar.ax.set_ylabel(r'Density $kg/m^{3}$', rotation=270, labelpad=20)
            cont = ax.contour(XX, ZZ, temperature,np.linspace(200,2000,10) , colors="white")
            fmt='%4.0f'
            plt.clabel(cont, inline=1,fmt=fmt,manual=False)
            plt.xlabel('Distance (km)')
            plt.ylabel('Depth (km)')
            rcParams['font.size'] = 10
            rcParams['axes.labelsize'] = 8
            rcParams['axes.titlesize'] = 8
            rcParams['legend.fontsize'] = 8
            rcParams['xtick.labelsize'] = 8
            rcParams['ytick.labelsize'] = 8

            plt.plot(dist/1e3,-1*elev_etopo1_raw/1e3,color='b',linewidth=0.75)

        sectionX    = np.linspace(lon1,lon2,20)
        sectionY    = np.linspace(lat1,lat2,20)
        dist        = np.zeros_like(sectionX)
        elev_etopo1 = np.zeros_like(sectionX)
        ages        = np.zeros_like(sectionX)
        spreading_rates = np.zeros_like(sectionX)
        ipt         = 0
        for lon,lat in zip(sectionX,sectionY):
            elev_etopo1[ipt]              = f_interp_etopo1([lat,lon])[0]
            ages[ipt]                     = f_interp_age([lat,lon])[0]
            spreading_rates[ipt]          = f_interp_spreading_rate([lat,lon])[0]
            azimuth1, azimuth2, dist[ipt] = geod.inv(lon1, lat1, lon, lat)
            ipt = ipt + 1
        
        if not plot_wintercg_data:
            cmap        = mpl.cm.jet
            cmaplist    = [cmap(i) for i in range(cmap.N)]
            cmaplist[0] = (.5, .5, .5, 1.0) # force the first color entry to be grey
            cmap        = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
            bounds      = np.linspace(0, 4, 9)
            norm        = mpl.colors.BoundaryNorm(bounds, cmap.N)
            pts         = plt.scatter(dist/1e3,elev_etopo1,cmap=cmap,c=ages,norm=norm)
            plt.colorbar(pts)

        plt.xlim(xlim)
        plt.ylim(ylim)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True))
        ax.yaxis.get_major_formatter().set_powerlimits((0, 1))

        if plot_wintercg_data:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
            plt.gca().invert_yaxis()
    
        plt.text(0.8,0.9,outputname.split('.')[0],transform=ax.transAxes) ;
        plt.text(0.02,0.9,'MOR elev {0:7.1f}'.format(np.nanmax(elev_etopo1)),transform=ax.transAxes) ;
        plt.text(0.02,0.8,'Full spreading rate (median) {0:4.1f} cm/yr'.format(np.nanmedian(spreading_rates)),transform=ax.transAxes) ;
        #plt.savefig(outputname,dpi=300)
        
    plt.show()
    
