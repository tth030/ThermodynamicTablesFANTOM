#!/usr/bin/env python

import numpy as np
import netCDF4
import os
import pickle
from statsmodels.stats.weightstats import DescrStatsW

def load_wintercg_crustal_thickness(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/',gemma=False):
    if os.path.exists('./data_figshare/wintercg'):
        path='./data_figshare/wintercg/'
    # Crustal thickness
    if gemma:
        filename = 'crust_thickness_GEMMA_0.5.grd'
    else:
        filename = 'Global_Moho_WINTERC-G_surf.grd'
    if (os.path.isfile(path+filename)):
        nc_crust = netCDF4.Dataset(path+filename)
    x=nc_crust.variables['lon'][:]
    y=nc_crust.variables['lat'][:]
    crustal_thickness=nc_crust.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        crustal_thickness=crustal_thickness.to_numpy()
    if gemma:
        print("Crustal thickness (Gemma) km min/max {} {}".format(np.nanmin(crustal_thickness),np.nanmax(crustal_thickness)))
    else:
        print("Crustal thickness (WINTERC-G) km min/max {} {}".format(np.nanmin(crustal_thickness),np.nanmax(crustal_thickness)))
    return x,y,crustal_thickness

def load_wintercg_lab(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/'):
    if os.path.exists('./data_figshare/wintercg'):
        path='./data_figshare/wintercg/'
    # Lithospheric thickness
    filename = 'Global_LAB_WINTERC-G_surf.grd'
    if (os.path.isfile(path+filename)):
        nc_lab = netCDF4.Dataset(path+filename)
    x=nc_lab.variables['lon'][:]
    y=nc_lab.variables['lat'][:]
    lab=nc_lab.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        lab=lab.to_numpy()
    print("LAB (WINTERC-G) km min/max {} {}".format(np.nanmin(lab),np.nanmax(lab)))
    return x,y,lab

def load_wintercg_crustal_density(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/'):
    if os.path.exists('./data_figshare/wintercg'):
        path='./data_figshare/wintercg/'
    # average crustal densities
    filename = 'rho_c_out_180-180_surf.grd'
    if (os.path.isfile(path+filename)):
        nc_rhoc = netCDF4.Dataset(path+filename)
    x=nc_rhoc.variables['lon'][:]
    y=nc_rhoc.variables['lat'][:]
    rhoc=nc_rhoc.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        rhoc=rhoc.to_numpy()
    print("Average crustal density (WINTERC-G) kg/m3 min/max {} {}".format(np.nanmin(rhoc),np.nanmax(rhoc)))
    return x,y,rhoc

def load_wintercg_density(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/density/',z=None,rhoref=False,alpha=3.2e-5,beta=0.77e-11,verbose=True,useCstAlphaBeta=False):
    if os.path.exists('./data_figshare/wintercg/density'):
        path='./data_figshare/wintercg/density/'

    alldepth = np.arange(-401,6,2)
    idepth   = np.argmin(np.abs(alldepth - z))
    # Density
    filename = 'd'+str(alldepth[idepth])+'.grd'
    if (os.path.isfile(path+filename)):
        nc_density = netCDF4.Dataset(path+filename)
    #print(nc_density)
    x=nc_density.variables['lon'][:]
    y=nc_density.variables['lat'][:]
    density=nc_density.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        density=density.to_numpy()
    if rhoref:
        x,y,temperature = load_wintercg_temperature(z=z,verbose=verbose)
        x,y,pressure    = load_wintercg_pressure(z=z,verbose=verbose)
        if not useCstAlphaBeta:
            Ax = 0;    Ay = 0;   Az = alpha;
            Bx = 1100; By = 6e9; Bz = alpha;
            Cx = 1304; Cy = 2e9; Cz = alpha+0.6e-5;
            a=(By-Ay)*(Cz-Az)-(Cy-Ay)*(Bz-Az)
            b=(Bz-Az)*(Cx-Ax)-(Cz-Az)*(Bx-Ax)
            c=(Bx-Ax)*(Cy-Ay)-(Cx-Ax)*(By-Ay)
            d=-1*(a*Ax+b*Ay+c*Az)
            #print(a,b,c,d)
            #d = 179968000
            #c = -5624e9
            #b = -0.0066
            #a = 36e3
            alpha = -1*(a*temperature+b*pressure+d)/c
    
            Ax = 0;    Ay = 0;     Az = beta;
            Bx = 1500; By = 6e9;   Bz = beta;
            Cx = 1304; Cy = 2e9;   Cz = beta+1.3e-12;
            a=(By-Ay)*(Cz-Az)-(Cy-Ay)*(Bz-Az)
            b=(Bz-Az)*(Cx-Ax)-(Cz-Az)*(Bx-Ax)
            c=(Bx-Ax)*(Cy-Ay)-(Cx-Ax)*(By-Ay)
            d=-1*(a*Ax+b*Ay+c*Az)
            #print(a,b,c,d)
            #d = 40.14
            #c = -5222e9
            #b = -1.95e-9
            #a = 0.00825
            beta     = -1*(a*temperature+b*pressure+d)/c
        density  = density/((1-temperature*alpha)*(1+pressure*beta))
    else:
        if verbose:
            print("Density (WINTERC-G) kg/m3 min/max {} {} at z= {} (closest from model z={})".format(np.nanmin(density),np.nanmax(density),z,alldepth[idepth]))
    return x,y,density

def load_wintercg_temperature(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/temperature/',z=None,verbose=True):
    if os.path.exists('./data_figshare/wintercg/temperature'):
        path='./data_figshare/wintercg/temperature/'

    alldepth = np.arange(-401,6,2)
    idepth   = np.argmin(np.abs(alldepth - z))
    # Density
    filename = 't'+str(alldepth[idepth])+'.grd'
    if (os.path.isfile(path+filename)):
        nc_temperature = netCDF4.Dataset(path+filename)
    x=nc_temperature.variables['lon'][:]
    y=nc_temperature.variables['lat'][:]
    temperature=nc_temperature.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        temperature=temperature.to_numpy()
    if verbose:
        print("Temperature (WINTERC-G) degC min/max {} {} at z= {} (closest from model z={})".format(np.nanmin(temperature),np.nanmax(temperature),z,alldepth[idepth]))
    return x,y,temperature

def load_wintercg_pressure(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/pressure/',z=None,verbose=True):
    if os.path.exists('./data_figshare/wintercg/pressure'):
        path='./data_figshare/wintercg/pressure/'

    alldepth = np.arange(-401,6,2)
    idepth   = np.argmin(np.abs(alldepth - z))
    # Density
    filename = 'p'+str(alldepth[idepth])+'.grd'
    if (os.path.isfile(path+filename)):
        nc_pressure = netCDF4.Dataset(path+filename)
    x=nc_pressure.variables['lon'][:]
    y=nc_pressure.variables['lat'][:]
    pressure=nc_pressure.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        pressure=pressure.to_numpy()
    if verbose:
        print("Pressure (WINTERC-G) Pa min/max {} {} at z= {} (closest from model z={})".format(np.nanmin(pressure),np.nanmax(pressure),z,alldepth[idepth]))
    return x,y,pressure

def build_wintercg_pressure(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/pressure/'):
    if os.path.exists('./data_figshare/wintercg/pressure'):
         path='./data_figshare/wintercg/pressure/'
    #load densities from top to bottom to compute lithostatic pressure
    alldepth = np.arange(-401,6,2)
    dh = 2000 # m
    g  = 10   # m/s2
    for depth in np.flip(alldepth):
        x,y,density = load_wintercg_density(z=depth)
        if (depth<5):
            x,y,pressureabove = load_wintercg_pressure(z=depth+2)
            x,y,densityabove  = load_wintercg_density(z=depth+2)
            pressureabove     = pressureabove + np.asarray(densityabove)*(dh/2)*g
        else:
            pressureabove = 0.
        pressure = np.asarray(density)*(dh/2)*g + pressureabove
        #pressure = np.array([[i for i in k] for k in pressure], dtype=np.float32)

        #print(pressure)

        # write file
        try: ncfile.close()  # just to be safe, make sure dataset is not already open.
        except: pass
        ncfile          = netCDF4.Dataset(path+'p'+str(depth)+'.grd',mode='w',format='NETCDF4')
        lat_dim         = ncfile.createDimension('lat', len(y)) # latitude axis
        lon_dim         = ncfile.createDimension('lon', len(x)) # longitude axis
        ncfile.title    = 'lithostatic pressure at depth '+str(depth)
        lon             = ncfile.createVariable('lon', np.float64, ('lon',))
        lon.units       = 'degrees_east'
        lon.long_name   = 'longitude'
        lat             = ncfile.createVariable('lat', np.float64, ('lat',))
        lat.units       = 'degrees_north'
        lat.long_name   = 'latitude'
        z               = ncfile.createVariable('z', np.float32, ('lat','lon'))
        z.units         = 'Pa'
        z.standard_name = 'pressure'
        z.long_name     = 'pressure'
        nlats           = len(lat_dim)
        nlons           = len(lon_dim)
        lat[:]          = np.arange(-90,90.5,0.5)
        lon[:]          = np.arange(-180,180.5,0.5)
        z[:,:]          = pressure
        #print(ncfile)
        ncfile.close()
        #break
        
    print("All .xyz are built")
    return 0

def load_wintercg_strain_rate(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/'):
    if os.path.exists('./data_figshare/wintercg/'):
        path='./data_figshare/wintercg/'
    # Strain rate
    nc_strain_rate = netCDF4.Dataset(path+'GSRM_strain_0.5.grd')
    #print(nc_strain_rate.variables.keys())
    x=nc_strain_rate.variables['lon'][:]
    y=nc_strain_rate.variables['lat'][:]
    strain_rate=nc_strain_rate.variables['z'][:]
    strain_rate=strain_rate*1e-9/(365*24*3600)
    strain_rate=np.where(np.isnan(strain_rate),1e-21,strain_rate)
    print("Strain rate (Kreemer et al, 2014) s-1 min/max {} {}".format(np.nanmin(strain_rate),np.nanmax(strain_rate)))
    return x,y,strain_rate

def load_wintercg_seafloor_ages(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/'):
    if os.path.exists('./data_figshare/wintercg/'):
        path='./data_figshare/wintercg/'
    # Seafloor age
    nc_age = netCDF4.Dataset(path+'age.3.2_0.5.nc')
    #print(nc_age.variables.keys())
    x=nc_age.variables['lon'][:]
    y=nc_age.variables['lat'][:]
    age=nc_age.variables['z'][:]
    age=age/100
    print("Seafloor age (Muller et al, 2008) Myrs min/max {} {}".format(np.nanmin(age),np.nanmax(age)))
    return x,y,age

def load_wintercg_hotspots(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/'):
    if os.path.exists('./data_figshare/wintercg/'):
        path='./data_figshare/wintercg/'
    # Hot spots
    filename = 'dist_closest_hs_0.5.grd'
    if (os.path.isfile(path+filename)):
        nc_dist_closest_hs = netCDF4.Dataset(path+filename)
    else:
        name               = filename.split('.grd')[0]
        nc_dist_closest_hs = xr.open_mfdataset([path+name+'_1.grd',path+name+'_2.grd',path+name+'_3.grd',path+name+'_4.grd',path+name+'_5.grd',path+name+'_6.grd'],
                                      concat_dim=['x'], combine='nested',engine='netcdf4')
    #print(nc_dist_closest_hs.variables.keys())
    x=nc_dist_closest_hs.variables['lon'][:]
    y=nc_dist_closest_hs.variables['lat'][:]
    dist_closest_hs=nc_dist_closest_hs.variables['z'][:]
    if not os.path.isfile(path+filename):
        x=x.to_numpy()
        y=y.to_numpy()
        dist_closest_hs=dist_closest_hs.to_numpy()
    print("Distances closest hot spot (Morgan and Morgan, 2007) m min/max {} {}".format(np.nanmin(dist_closest_hs),np.nanmax(dist_closest_hs)))
    
    return x,y,dist_closest_hs


def load_wintercg_temp_profiles(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/',
                                mask=None,coord=None,weights=None,
                                ridges=False, selection_name='MOR_pts_all', distance_between_pts_along_ridges = 25, distance_to_ridges=0):
    if os.path.exists('./data_figshare/wintercg/'):
         path='./data_figshare/wintercg/'

    profiles = []
    profile_avg = []
    profile_med = []
    profile_std = []
    alldepth = np.arange(-401,6,2)
    x,y,temperature = load_wintercg_temperature(z=5)
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    if mask is None:
        mask=(temperature!=temperature)

    if ridges:
        path_ridges = './data/topography/'
        if os.path.exists('./data_figshare/topography/'):
             path_ridges='./data_figshare/topography/'
        filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'
        if not (os.path.isfile(path_ridges+filename)):
            print("You did not select MOR points using `select_MOR_pts()` or selection_name is wrong ({})".format(path_ridges+filename))
            exit()
        else:
            with open(path_ridges+filename, 'rb') as filehandle:
                # read the data as binary data stream
                [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)
            coord      = None
            npts       = len(new_active_MOR_x)
            #print("Ridges npts = {}".format(npts))
            xi         = [] ; yi = []
            mask_saved = mask
            mask[:,:]  = True
            for ipt in np.arange(0,npts):
                xipt = np.argmin(np.abs(x - new_active_MOR_x[ipt]))
                yipt = np.argmin(np.abs(y - new_active_MOR_y[ipt]))
                xi.append(xipt)
                yi.append(yipt)
                mask[yipt,xipt] = False
            mask = np.where((~mask_saved),mask,mask_saved)
            npts = np.count_nonzero(mask==False)
                
    if coord:
        npts       = len(coord)
        mask_saved = mask
        mask[:,:]  = True
        xi = [] ; yi = []
        for ipt in np.arange(0,npts):
            xipt = np.argmin(np.abs(x - coord[ipt][0]))
            yipt = np.argmin(np.abs(y - coord[ipt][1]))
            xi.append(xipt)
            yi.append(yipt)
            mask[yipt,xipt] = False
        mask = np.where((~mask_saved),mask,mask_saved)
        npts = np.count_nonzero(mask==False)
    else:
        npts = np.count_nonzero(mask==False)
        xi = np.arange(0,len(x))
        yi = np.arange(0,len(y))

    temp = []
    for depth in np.flip(alldepth):
        x,y,T = load_wintercg_temperature(z=depth,verbose=False)
        temp.append(T)

    print("npts {}".format(npts))

    ipt = 0
    for j in yi:
        for i in xi: 
            prof = []
            if not mask[j,i]:
                for k in np.arange(0,len(alldepth)):
                    temperature = temp[k]
                    # mask
                    data = temperature
                    data = np.ma.masked_array(data, mask)
                    prof.append(data[j,i])
                    if ipt==0:
                        mean   = np.nanmean(data)
                        median = np.ma.median(data)
                        sigma  = np.nanstd(data)

                        if weights is not None:
                            ma      = np.ma.MaskedArray(data, mask=np.isnan(data))
                            meanW   = np.ma.average(ma.reshape(-1),weights=weights)
                            dsw     = DescrStatsW(ma.reshape(-1),weights=weights)
                            stdW    = dsw.std  # weighted std
                            mean    = meanW
                            sigma   = stdW

                        profile_avg.append(mean)
                        profile_med.append(median)
                        profile_std.append(sigma)
                profiles.append(prof)
                ipt = ipt + 1
                if (ipt%1000==0):
                    print(ipt)
    return profiles, profile_avg, profile_med, profile_std, alldepth
                    
def load_wintercg_density_profiles(path='/Users/thomastheunissen/Documents/Recherche/Norway/data/WINTERC-G/',
                                   mask=None,coord=None,weights=None,
                                   rhoref=False,useCstAlphaBeta=False,alpha=3.2e-5,beta=0.77e-11,
                                   lab=0,extractLABthickness=False,
                                   ridges=False, selection_name='MOR_pts_all', distance_between_pts_along_ridges = 25, distance_to_ridges=0):
    if os.path.exists('./data_figshare/wintercg/'):
         path='./data_figshare/wintercg/'
    profiles = []
    profile_avg = []
    profile_med = []
    profile_std = []
    alldepth = np.arange(-401,6,2)
    x,y,density = load_wintercg_density(z=5,rhoref=rhoref,verbose=False,useCstAlphaBeta=useCstAlphaBeta,alpha=alpha,beta=beta)
    if mask is None:
        mask=(density!=density)

    if (lab>0):
        x,y,lithos_thickness  = load_wintercg_lab()
        x,y,crustal_thickness = load_wintercg_crustal_thickness(gemma=False)
        lithos_thickness      = -1*lithos_thickness
        crustal_thickness     = -1*crustal_thickness

    if ridges:
        path_ridges = './data/topography/'
        if os.path.exists('./data_figshare/topography/'):
             path_ridges='./data_figshare/topography/'
        filename = './'+selection_name+'_'+str(distance_between_pts_along_ridges)+'km.dat'
        if not (os.path.isfile(path_ridges+filename)):
            print("You did not select MOR points using `select_MOR_pts()` or selection_name is wrong ({})".format(path_ridges+filename))
            exit()
        else:
            with open(path_ridges+filename, 'rb') as filehandle:
                # read the data as binary data stream
                [new_active_MOR_x,new_active_MOR_y,new_active_MOR_elev,new_active_MOR_spreading_rate] = pickle.load(filehandle)
            coord      = None
            npts       = len(new_active_MOR_x)
            #print("Ridges npts = {}".format(npts))
            xi         = [] ; yi = []
            mask_saved = mask
            mask[:,:]  = True
            for ipt in np.arange(0,npts):
                xipt = np.argmin(np.abs(x - new_active_MOR_x[ipt]))
                yipt = np.argmin(np.abs(y - new_active_MOR_y[ipt]))
                xi.append(xipt)
                yi.append(yipt)
                mask[yipt,xipt] = False
            mask = np.where((~mask_saved),mask,mask_saved)
            npts = np.count_nonzero(mask==False)
 
    if coord:
        npts       = len(coord)
        mask_saved = mask
        mask[:,:]  = True
        xi = [] ; yi = []
        for ipt in np.arange(0,npts):
            xipt = np.argmin(np.abs(x - coord[ipt][0]))
            yipt = np.argmin(np.abs(y - coord[ipt][1]))
            xi.append(xipt)
            yi.append(yipt)
            mask[yipt,xipt] = False
        mask = np.where((~mask_saved),mask,mask_saved)
        npts = np.count_nonzero(mask==False)
    else:
        npts = np.count_nonzero(mask==False)
        xi = np.arange(0,len(x))
        yi = np.arange(0,len(y))


    allrho = []
    for depth in np.flip(alldepth):
        x,y,rho = load_wintercg_density(z=depth,verbose=False,rhoref=rhoref,useCstAlphaBeta=useCstAlphaBeta,alpha=alpha,beta=beta)
        allrho.append(rho)
    allrho = np.flip(allrho)

    print("npts {}".format(npts))

    ipt = 0
    diff_average = []
    allLABthickness = []
    for j in yi:
        for i in xi: 
            prof = []
            if not mask[j,i]:
                allLABthickness.append(lithos_thickness[j,i])
                for k in np.arange(0,len(alldepth)):
                    density = allrho[k]
                    # mask
                    data = density
                    data = np.ma.masked_array(data, mask)
                    if lab==1:
                        if alldepth[k]<=lithos_thickness[j,i]:
                            prof.append(data[j,i])
                        else:
                            prof.append(np.nan)

                    if lab==2:
                        if ((alldepth[k]>=lithos_thickness[j,i]) and (alldepth[k]<crustal_thickness[j,i]-2.5)):
                            prof.append(data[j,i])
                        else:
                            prof.append(np.nan)

                    if lab==3:
                        if (alldepth[k]>crustal_thickness[j,i]):
                            prof.append(data[j,i])
                        else:
                            prof.append(np.nan)

                    if lab==0:
                        prof.append(data[j,i])

                    if ipt==0:
                        mean   = np.nanmean(data)
                        median = np.ma.median(data)
                        sigma  = np.nanstd(data)

                        if weights is not None:
                            ma      = np.ma.MaskedArray(data, mask=np.isnan(data))
                            meanW   = np.ma.average(ma.reshape(-1), weights=weights)
                            dsw     = DescrStatsW(ma.reshape(-1), weights=weights)
                            stdW    = dsw.std  # weighted std
                            mean    = meanW
                            sigma   = stdW

                        profile_avg.append(mean)
                        profile_med.append(median)
                        profile_std.append(sigma)
                    
                    diff_average.append(prof[k]-profile_avg[k])
                profiles.append(np.flip(prof))
                ipt = ipt + 1
                if (ipt%1000==0):
                    print(ipt)

    avg_diff = np.nanmean(diff_average)
    med_diff = np.nanmedian(diff_average)
    std_diff = np.nanstd(diff_average)
    print("Average difference density with respect to the average is {} +/- {} median = {}".format(avg_diff,std_diff,med_diff))

    if extractLABthickness:
        return profiles, np.flip(profile_avg), np.flip(profile_med), np.flip(profile_std), alldepth, allLABthickness
    else:
        return profiles, np.flip(profile_avg), np.flip(profile_med), np.flip(profile_std), alldepth
     
