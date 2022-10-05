import numpy as np
from mpl_toolkits.basemap import Basemap
#import basemap as Basemap
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pickle
import os
import csv

from .external import *

def plot_MOR(what='elevation',path='./data/topography/',selection_name='MOR_pts_all',distance_between_pts_along_ridges = 25,savefig = False):
    '''  
        This routine plots elevation or spreading rate of a selection of points located along mid-ocean ridges
             what           = 'spreadingRate' or 'elevation'
             selection_name =  'MOR_pts_all' or 'MOR_pts_far_from_hs' or 'MOR_pts_close_to_hs'
    '''
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'

    define_rcParams()

    filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'

    if not (os.path.isfile(path+filename)):
        print("You did not select MOR points using `select_MOR_pts()` or selection_name is wrong ({})".format(path+filename))
        exit()
    else:
        with open(path+filename, 'rb') as filehandle:
            # read the data as binary data stream
            [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)
   
    plt.figure(figsize=(25/2.54,12.5/2.54))

    # miller projection
    map = Basemap(projection='mill',lon_0=0)
    map.drawcoastlines(linewidth=0.5)
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
    map.drawmapboundary()
    #xx, yy = np.meshgrid(x, y)
    #map.pcolormesh(xx, yy, active_MOR,latlon=True, cmap='jet',vmin=-4500,vmax=0)
    
    filename2='Morgan-Morgan_2007_hs.txt'
    data = np.loadtxt(path+filename2,dtype='float',delimiter=' ')
    lat_read = data[:,1]
    lon_read = data[:,0]
    #x,y = map(lon_read, lat_read)
    #map.plot(x, y, 'bo', markersize=24)
    radius = 1e3 #km
    for lon,lat in zip(lon_read, lat_read):
        equi(map,lon,lat,radius,lw=0.5,color='black')
    x,y = map(new_active_MOR_x, new_active_MOR_y)
    cmap = mpl.cm.jet
    if what == 'spreadingRate':
        norm = mpl.colors.Normalize(vmin=12,vmax=30)
        pts = plt.scatter(x, y, marker='o', cmap=cmap,c=new_active_MOR_spreading_rate,norm=norm)
    elif what == 'elevation':
        norm = mpl.colors.Normalize(vmin=-4500.,vmax=0.) 
        pts = plt.scatter(x, y, marker='o', cmap=cmap,c=new_active_MOR_elev,norm=norm)
    else:
        print("what = '{}' does not exist".format(what))
        pass
    
    lon_read = [] ; lat_read = []
    with open(path+'MOR_examples_selection-Fig3.txt', mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file,delimiter=' ')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                pass
            lon_read.append(float(row['lon'])) ; lat_read.append(float(row['lat']))
            line_count += 1
    x,y = map(lon_read, lat_read)
    map.plot(x, y, 'bo', markersize=1)
    plt.colorbar(pts)
   
    if savefig:
        if what == 'spreadingRate':
            plt.savefig('map_active_MOR_spreading_rate.pdf',dpi=300)
        elif what == 'elevation':
            plt.savefig('map_active_MOR_elevation.pdf',dpi=300)
        else:
            print("what = '{}' does not exist".format(what))
            pass

    plt.show()


def plot_global_maps_continents(x, y, data,
                                path                  = './data/topography/',
                                typedat               = 'elevation',
                                plotcircles_hs        = True,
                                orientationbar        = 'vertical',
                                savefig              = False,
                                vmin                 = None,
                                vmax                 = None):
    '''
        x,y : 1-D arrays
        data : 2-D arrays
        typedat = "distance_1st_hs", "seafloorage", "strainrate",
                  "elevation", "crustal_thickness", "lithos_thickness",
                  "age_lithos", "crustal_density", "mantle_density"
        orientationbar = 'horizontal', 'vertical'
    '''
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'

    define_rcParams()

    fig = plt.figure(figsize=(25/2.54,12.5/2.54))
    #fig = plt.figure(figsize=(50/2.54,25/2.54))
    # miller projection
    map = Basemap(projection='mill',lon_0=0)

    xpt,ypt = map(*np.meshgrid(x,y))

    if vmin is None:
        vmin = np.min(data)
    if vmax is None:
        vmax = np.max(data)

    if typedat == "elevation":
        cdict = gmtColormap('ETOPO1',GMTPath=path)
        cmap  = LinearSegmentedColormap('ETOPO1', cdict)
        vmin  = -7500  # -3500
        vmax  = 6000    # 2800
        pcm   = plt.pcolormesh(xpt,ypt,data,cmap=cmap,vmin=vmin,vmax=vmax,rasterized=True,shading='nearest')
        pcm.cmap.set_over('whitesmoke')
    elif ((typedat == "crustal_density") | (typedat == "mantle_density") | (typedat == "temperature") | (typedat == "pressure") | (typedat == "density")):
        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,vmin=vmin,vmax=vmax,rasterized=True,shading='nearest')
#    elif typedat == "mantle_density":
#        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,vmin=3350,vmax=3550,rasterized=True,shading='nearest')
    elif typedat == "lithos_thickness":
        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,vmin=40,vmax=250,rasterized=True,shading='nearest')
    elif typedat == "crustal_thickness":
        # GEMMA
        #if (mask>=1):
        #    pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,vmin=25,vmax=45,rasterized=True,shading='nearest')
        #else:
        if vmin is None or vmax is None:
            bounds = np.linspace(0,65,14)
        else:
            bounds = np.linspace(vmin,vmax,14)
        norm   = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pcm    = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,norm=norm,rasterized=True,shading='nearest')
    elif typedat == "strainrate":
        norm = mpl.colors.SymLogNorm(linthresh=1e-16, linscale=1e-16,
                                              vmin=1e-17, vmax=1e-14,base=np.e)
        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,norm=norm,rasterized=True,shading='nearest')
    elif typedat == "seafloorage":
        bounds = np.linspace(0,280,29)
        norm   = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet_r,norm=norm,rasterized=True,shading='nearest')
    elif typedat == "distance_1st_hs":
        bounds = np.linspace(0,4000,41)
        norm   = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pcm = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet_r,norm=norm,rasterized=True,shading='nearest')
    elif typedat == "age_lithos":
        bounds = np.linspace(0,4000,41)
        norm   = mpl.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
        pcm    = plt.pcolormesh(xpt,ypt,data,cmap=plt.cm.jet,norm=norm,rasterized=True,shading='nearest')


    if typedat == "strainrate":
        map.drawcoastlines(linewidth=0.25,)
    else:
        map.drawcoastlines(linewidth=0.5,)
    map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
    map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
    map.drawmapboundary()

    if plotcircles_hs:
        filename2='Morgan-Morgan_2007_hs.txt'
        hs = np.loadtxt(path+filename2,dtype='float',delimiter=' ')
        lat_read = hs[:,1]
        lon_read = hs[:,0]
        #x,y = map(lon_read, lat_read)
        #map.plot(x, y, 'bo', markersize=24)
        radius = 1e3 #km
        for lon,lat in zip(lon_read, lat_read):
            if typedat == "strainrate":
                equi(map,lon,lat,radius,lw=0.25,color='black',zorder=4)
            else:
                equi(map,lon,lat,radius,lw=0.5,color='black',zorder=4)
  
    if orientationbar=='horizontal':
        shrink=0.5
        aspect=25
    else:
        shrink=1
        aspect=20

    if typedat == "elevation":
        if vmax==2800:
            plt.colorbar(pcm,extend='max',orientation=orientationbar,shrink=shrink,aspect=aspect)
        else:
            plt.colorbar(pcm,orientation=orientationbar,shrink=shrink,aspect=aspect)
    elif typedat == "crustal_thickness":
        plt.colorbar(pcm,extend='both',orientation=orientationbar,shrink=shrink,aspect=aspect)
    elif typedat == "strainrate":
        plt.colorbar(pcm,extend='both',orientation=orientationbar,shrink=shrink,aspect=aspect)
    elif ((typedat == "crustal_density") | (typedat == "mantle_density") | (typedat == "temperature") | (typedat == "pressure") | (typedat == "density")):
        plt.colorbar(pcm,extend='both',orientation=orientationbar,shrink=shrink,aspect=aspect)
    else:
        plt.colorbar(pcm,orientation=orientationbar,shrink=shrink,aspect=aspect)

#    ax=plt.gca()
#    cax,kw = mpl.colorbar.make_axes(ax,location=locationbar)
#    #mpl.colorbar.Colorbar(cax,pcm,extend='both')
#    plt.colorbar(mappable=pcm, cax=cax, extend='both')

    if savefig:
        plt.savefig('map_'+typedat+'.pdf',dpi=300)
    else:
        plt.show()

    return data[~data.mask]




def plot_histo_elev_MOR(path='./data/topography/',selection_name='MOR_pts_all',distance_between_pts_along_ridges = 25,
                        GaussianModel=True,sigmamodel=250,meanmodel=-2750,
                        savefig=False):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'

    define_rcParams()

    filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'

    with open(path+filename, 'rb') as filehandle:
        # read the data as binary data stream
        [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)

    data             = new_active_MOR_elev
    title            = 'Active ridges ('+filename+')'
    plot_histo(data,title,'figure.pdf',xlabel='Elevation (m)',unit='m',GaussianModel=True,sigmamodel=sigmamodel,meanmodel=meanmodel,savefig=savefig,text="right")
    plt.show()

def plot_correlation_elev_spreading_rate(path='./data/topography/',selection_name='MOR_pts_all',distance_between_pts_along_ridges = 25,
                                         savefig=False):
    if os.path.exists('./data_figshare/topography/'):
        path='./data_figshare/topography/'

    define_rcParams()

    filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'

    with open(path+filename, 'rb') as filehandle:
        # read the data as binary data stream
        [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)

    y = new_active_MOR_elev
    x = new_active_MOR_spreading_rate*2/10

    title    = 'MOR elevation nostrain vs spreading rate'
    plot_correlation(x,y,title,'figure.pdf',nbins = 20,
                     ylabel='Elevation (m)',
                     xlabel='Full-spreading rate (cm/yr)',
                     xlim=[0,20],
                     ylim=[-4000,0],
                     savefig=savefig,
                     ticks=[2,1,1000,250])
    plt.show()


