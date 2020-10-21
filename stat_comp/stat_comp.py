# stat_comp.py
'''
 Selection of argo profiles to compare to a simulation
 In this script the argo database is browsed in order to find the profiles that will be compare to the outputs of one simulation. 
 The comparison will not be a simple colocation of the profile inside the model grid but we want to make a statistical comparison of the observed profile with a significant number of profiles close to it in the model (close in terms of space and time, for instance in a 0.5° radius around the profile location and 10 days before and after it has been sampled)
 Therefore the selected profiles must have to be relevant according to some criteria :
   - they must be inside the domain of simulation (the mask file of the configuration file must be provided)  
   - they must be sampled inside the period of simulation (the period shortened by a certain amount of days must be provided)
   - they must go as deep as a given depth (according to the desired depth for the comparison profiles)
   - they must be in a location where there are enough model profiles (for instance if profile is too close to an island or the coast)
'''
# Imports of librairies
import numpy as np

import time
import sys
import os
import glob
import io

import datetime
from datetime import date

from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import dask

import json
import re

import seawater
import cmocean

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import warnings
warnings.filterwarnings('ignore')

# Define the criteria for one profile

# Localisation of a given profile location inside model domain
def loc(latargo,lonargo,nprof,meshfile,namlatmod,namlonmod,nammaskmod,sosie_exec):
    # open the maskfile and get lat lon and mask
    ds=xr.open_dataset(meshfile)
    lat=ds[namlatmod]
    lon=ds[namlonmod]
    latmin,latmax,lonmin,lonmax=(lat.min(),lat.max(),lon.min(),lon.max())
    
    if (lonmin.values<lonargo)&(lonargo<lonmax.values)&(latmin.values<latargo)&(latargo<latmax.values):
        if not os.path.exists('txt'):
            os.makedirs('txt')    
        with open('txt/prof0.txt','w') as txt_file:
            txt_file.write('Profile_'+str(nprof)+' '+str(lonargo)+' '+str(latargo))
        if os.path.exists('ij_found.out'):
            get_ipython().system('rm ij_found.out')
            
        os.system(str(sosie_exec)+' -i '+str(meshfile)+' -p txt/prof0.txt  > txt/output')
        if os.path.exists('ij_found.out'):
            with open('ij_found.out','r') as txt_file:
                last_line = txt_file.readlines()[-1]
                if last_line[1] == '#':
                    print('Profile no '+str(nprof)+' is not in the domain, do not keep')
                    i0,j0=-1,-1
                else:
                    print('Profile no '+str(nprof)+' is in the domain, go proceed')
                    i0=int(last_line.split()[1])
                    j0=int(last_line.split()[2])
        else:
            print('Problem with ij_found program, try debugging '+str(sosie_exec)+' -i '+str(meshfile)+' -p txt/prof0.txt')
    else:
        print('Profile no '+str(nprof)+' is not in the domain, do not keep')
        i0,j0=-1,-1
       
    return i0,j0


# Check if profile location is on land
def check_prof_in_ocean(i0,j0,meshfile,nammaskmod):
    print('check if profile is in the ocean : ')
    dsN=xr.open_dataset(meshfile)
    tmaskN=dsN[nammaskmod] 
    if tmaskN[0,0,int(j0)-1,int(i0)-1] == 1:
        check=0
    else:
        check=1
    return check

# Check if too close to boundaries of model domain
def check_close_to_boundaries(nprof,maskfile,nammaskmod,i0,j0,radius_max):
    print('check if profile is not too close to model boundaries : ')
    dsm=xr.open_dataset(maskfile)
    tmask=dsm[nammaskmod][0,0]
    ly,lx=tmask.shape
    gdpts=np.int(np.round(radius_max*60))
    if (j0-gdpts<0) or (j0+gdpts>ly-1) or (i0-gdpts<0) or (i0+gdpts>lx-1):
        check=1
    else:
        check=0
    return check
            
# Check depth of profile
def check_prof_depth(nprof,presargo,latargo,depthmin):
    print('check if profile has a good depth : ')
    #convert pressure to depth
    depthargo=seawater.dpth(presargo,latargo)
    #look for the last level
    indzprof=np.max(np.where(np.isnan(depthargo)==False))
    dmax=depthargo[indzprof]
    print('profile max depth is '+str(dmax)+' m')
    if dmax >= depthmin:
        check=0
    else:
        check=1
    return check


# Check if there are enough model profiles around the obs profile
def map_profile_from_jsonfile(lonprof,latprof,radius,lonmodmin, lonmodmax,latmodmin, latmodmax,nprof,plotdir):
    # Produce a map with the profile location and the area in which the model profiles are sampled
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.patches as mpatches
    fig=plt.figure(figsize=(20,15))
    ax = plt.subplot(111,projection=ccrs.PlateCarree(central_longitude=0))
    ax.set_extent((lonmodmin, lonmodmax, latmodmin, latmodmax))
    ax.coastlines(resolution="10m")
    gl = ax.gridlines(draw_labels=True, linestyle=':', color='black',
                      alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.tick_params('both',labelsize=22)
    ax.add_patch(mpatches.Circle(xy=[lonprof,latprof], radius=radius, color='green', alpha=0.3, transform=ccrs.PlateCarree(), zorder=30))
    plt.scatter(lonprof,latprof, c='g', linewidth='0', s=18);
    plt.savefig(plotdir+'/debug_map_'+str(nprof)+'.png')

def check_number_profile(nprof,i0,j0,depthmin,coordfile,maskfile,zgrfile,namlatmod,namlonmod,namdepmod,nammaskmod,lonargo,latargo,radius_max,number_of_model_profiles,period,lonmodmin,lonmodmax,latmodmin,latmodmax,dmap=0):
    print('check if there are enough model profiles : ')
    if radius_max > 0:
        #open mask file and read lat, lon and mask
        ds=xr.open_dataset(coordfile)
        gdpts=np.int(np.round(radius_max*60))      
        lat=ds[namlatmod][j0-gdpts:j0+gdpts+1,i0-gdpts:i0+gdpts+1]
        lon=ds[namlonmod][j0-gdpts:j0+gdpts+1,i0-gdpts:i0+gdpts+1]
        dsz=xr.open_dataset(zgrfile)
        depth=dsz[namdepmod][0]
        dsm=xr.open_dataset(maskfile)
        tmask=dsm[nammaskmod][0,:,j0-gdpts:j0+gdpts+1,i0-gdpts:i0+gdpts+1]
        # Stack the variables
        lon_stacked = lon.stack(profile=('x', 'y'))
        lat_stacked = lat.stack(profile=('x', 'y'))
        mask_stacked = tmask.stack(profile=('x', 'y'))
        #Get the depth at every grid point
        d,ly,lx=tmask.shape
        depthmod2d=np.zeros([lx,ly])
        for j in np.arange(ly):
            for i in np.arange(lx):
                depthmod2d[j,i]=depth[np.min(np.where(tmask[:,j,i].values<1))].values
        xr_depthmod2d=xr.DataArray(depthmod2d, dims=("y", "x"))    
        depth_stacked = xr_depthmod2d.stack(profile=('x', 'y'))
        #find the profiles filling criteria
        distance_threshold = radius_max
        square_distance_to_observation = (lon_stacked - lonargo)**2 + (lat_stacked-latargo)**2
        is_close_to_observation = (square_distance_to_observation < distance_threshold) & (depth_stacked > depthmin)
        nb_profiles=np.sum(1*is_close_to_observation)*24*(period*2+1)
        print('There is a total of '+str(nb_profiles.values)+' model oceanic profiles with enough depth')
        if dmap == 1:
            map_profile_from_jsonfile(lonprof,latprof,radius,lonmodmin, lonmodmax,latmodmin, latmodmax,nprof,plotdir)
    else:
        nb_profiles=1
        print('There is a total of '+str(nb_profiles)+' model oceanic profiles with enough depth')
    if nb_profiles >= number_of_model_profiles:
        check=0
    else:
        check=1
    return check



# Make the selection of profiles
def selection(config='MEDWEST60',case='BLBT02',member='',dirmod='/mnt/alberta/equipes/IGE/meom/workdir/lerouste/MEDWEST60/',coordfile='/mnt/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_coordinates_v3.nc4',maskfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mask.nc4',hgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_hgr.nc4',zgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_zgr.nc4',batfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4',namlatmod='nav_lat',namlonmod='nav_lon',namdepmod='gdept_1d',nammaskmod='tmask',namtempmod='votemper',
namsaltmod='vosaline',ymin=2009,mmin=7,dmin=1,ymax=2010,mmax=6,dmax=30,depthmin=0,radius_max=0.25,period=0,number_of_model_profiles=100000,plotdir='plots',ncdir='/gpfswork/rech/egi/rote001/ARGO',dmap=1,sosie_exec='/mnt/meom/workdir/alberta/DEV/sosie/bin/ij_from_lon_lat.x',srcargo='erddap'):
    # determining the dates
    from argopy import DataFetcher as ArgoDataFetcher

    datemin=datetime.date(ymin,mmin,dmin)+datetime.timedelta(days=period)
    datemax=datetime.date(ymax,mmax,dmax)-datetime.timedelta(days=period)
    # determining the area
    ds=xr.open_dataset(coordfile)
    lat=ds[namlatmod]
    lon=ds[namlonmod]
    latmodmin,latmodmax,lonmodmin,lonmodmax=(lat.min(),lat.max(),lon.min(),lon.max())
    # use argopy to get the selection of profiles
    ds_points=ArgoDataFetcher(src=srcargo).region([lonmodmin.values,lonmodmax.values,latmodmin.values,latmodmax.values,0,12000,str(datemin),str(datemax)]).to_xarray()
    ds_profiles=ds_points.argo.point2profile()
    index_i_model=[]
    index_j_model=[]
    # get rid of profiles not following criteria
    for nprof in range(len(ds_profiles.N_PROF)):
        print('Processing profile no '+str(nprof))
        latargo=ds_profiles.LATITUDE[nprof].values
        lonargo=ds_profiles.LONGITUDE[nprof].values
        presargo=ds_profiles.PRES[nprof,:].values
        i0,j0=loc(latargo,lonargo,nprof,hgrfile,namlatmod,namlonmod,nammaskmod,sosie_exec)
        if (i0,j0) == (-1,-1):
            print('profile is not in the domain at all')
            continue
        check=check_prof_in_ocean(i0,j0,maskfile,nammaskmod)
        if check == 1:
            print('no, profile is on the land')
            continue
        print('yes, profile is in the ocean')
        check=check_close_to_boundaries(nprof,maskfile,nammaskmod,i0,j0,radius_max)
        if check == 1:
            print('no, profile is too close to model boundaries')
            continue
        print('yes, profile is not too close to model boundaries')
        check=check_prof_depth(nprof,presargo,latargo,depthmin)
        if check == 1:
            print('no, profile is not deep enough')
            continue
        print('yes, profile is deep enough')
        check=check_number_profile(nprof,i0,j0,depthmin,coordfile,maskfile,zgrfile,namlatmod,namlonmod,namdepmod,nammaskmod,lonargo,latargo,radius_max,number_of_model_profiles,period,lonmodmin, lonmodmax,latmodmin, latmodmax,dmap=0)
        if check == 1:
            print('no, there are not enough model profiles')
            continue
        print('yes, there are enough model profiles')
        index_i_model.append(i0)
        index_j_model.append(j0)
        ds_one=ds_profiles.isel(N_PROF=nprof)
        ds_prof=ds_one.expand_dims({'N_PROF':1})
        try:
            ds_profiles_out=xr.concat([ds_profiles_out,ds_prof],dim='N_PROF')
        except NameError:
            ds_profiles_out=ds_prof  
    indexI_da=xr.DataArray(index_i_model,dims='N_PROF')
    indexJ_da=xr.DataArray(index_j_model,dims='N_PROF')
    ds_profiles_out['index_i_model']=indexI_da
    ds_profiles_out['index_j_model']=indexJ_da
    if not os.path.exists(ncdir):
        os.makedirs(ncdir)
    netcdf_file=ncdir+'/ARGO_profiles_selection_for_'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.nc'
    ds_profiles_out.to_netcdf(path=netcdf_file,mode='w')
    return ds_profiles_out

                        
# Colocation of profile in model with hourly outputs of gridT and gridS

def colocation_profiles_argo(ds_profiles,config='MEDWEST60',case='BLBT02',member='',dirmod='/mnt/alberta/equipes/IGE/meom/workdir/lerouste/MEDWEST60/',coordfile='/mnt/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_coordinates_v3.nc4',maskfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mask.nc4',hgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_hgr.nc4',zgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_zgr.nc4',batfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4',namlatmod='nav_lat',namlonmod='nav_lon',namdepmod='gdept_1d',nammaskmod='tmask',namtempmod='votemper',
namsaltmod='vosaline',ymin=2009,mmin=7,dmin=1,ymax=2010,mmax=6,dmax=30,depthmin=0,radius_max=0.25,period=0,number_of_model_profiles=100000,plotdir='plots',ncdir='/gpfswork/rech/egi/rote001/ARGO',dmap=1,sosie_exec='/mnt/meom/workdir/alberta/DEV/sosie/bin/ij_from_lon_lat.x'):
    nb_profilesargo=len(ds_profiles.N_PROF)
    profile_temp_model=[]
    profile_salt_model=[]
    for prof in np.arange(nb_profilesargo):
        i0=ds_profiles.index_i_model[prof]
        j0=ds_profiles.index_j_model[prof]
        date_argo=ds_profiles.TIME[prof]
        t_argo=pd.to_datetime(date_argo.values)
        year_argo=t_argo.year
        month_argo=t_argo.month
        mm_argo = "{:02d}".format(month_argo)
        day_argo=t_argo.day
        dd_argo = "{:02d}".format(day_argo)
        hour_argo=t_argo.hour
        model_file_nameT=str(config)+'-'+str(case)+'_y'+str(year_argo)+'m'+str(mm_argo)+'d'+str(dd_argo)+'.1h_gridT.nc'
        model_file_nameS=str(config)+'-'+str(case)+'_y'+str(year_argo)+'m'+str(mm_argo)+'d'+str(dd_argo)+'.1h_gridS.nc'
        if os.path.exists(dirmod+model_file_nameT):
            print('the model file exists for this profile, proceed')
            dsT=xr.open_dataset(dirmod+model_file_nameT)
            dsS=xr.open_dataset(dirmod+model_file_nameS)
            model_profile_T=dsT[namtempmod][hour_argo,:,j0,i0].values
            model_profile_S=dsS[namsaltmod][hour_argo,:,j0,i0].values
            profile_temp_model.append(model_profile_T)
            profile_salt_model.append(model_profile_S)
            ds_one=ds_profiles.isel(N_PROF=prof)
            ds_prof=ds_one.expand_dims({'N_PROF':1})
            try:
                ds_profiles_out=xr.concat([ds_profiles_out,ds_prof],dim='N_PROF')
            except NameError:
                ds_profiles_out=ds_prof  
        else:
            print('the model file does not exist for this profile')
            continue
    dsz=xr.open_dataset(zgrfile)
    ds_profiles_out['DEPTH_MOD']=dsz[namdepmod][0]
    ds_profiles_out=ds_profiles_out.rename({'z':'DEPTH_MOD'})
    ds_profiles_out=ds_profiles_out.rename_vars({'DEPTH_MOD':'model_levels'})
    profile_temp_model_da=xr.DataArray(profile_temp_model,dims={'N_PROF','DEPTH_MOD'})
    profile_salt_model_da=xr.DataArray(profile_salt_model,dims={'N_PROF','DEPTH_MOD'})
    ds_profiles_out['profile_temp_model']=profile_temp_model_da
    ds_profiles_out['profile_salt_model']=profile_salt_model_da
    if not os.path.exists(ncdir):
        os.makedirs(ncdir)
    netcdf_file=ncdir+'/ARGO_profiles_selection_and_model_colocation_for_'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.nc'
    ds_profiles_out.to_netcdf(path=netcdf_file,mode='w')
    return ds_profiles_out
            
            


# Plot the locations of all profiles

def plot_profiles_argo(ds_profiles,config='MEDWEST60',case='BLBT02',member='',dirmod='/mnt/alberta/equipes/IGE/meom/workdir/lerouste/MEDWEST60/',coordfile='/mnt/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_coordinates_v3.nc4',maskfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mask.nc4',hgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_hgr.nc4',zgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_zgr.nc4',batfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4',namlatmod='nav_lat',namlonmod='nav_lon',namdepmod='gdept_1d',nammaskmod='tmask',namtempmod='votemper',
namsaltmod='vosaline',ymin=2009,mmin=7,dmin=1,ymax=2010,mmax=6,dmax=30,depthmin=0,radius_max=0.25,period=0,number_of_model_profiles=100000,plotdir='plots',ncdir='/gpfswork/rech/egi/rote001/ARGO',dmap=1,sosie_exec='/mnt/meom/workdir/alberta/DEV/sosie/bin/ij_from_lon_lat.x'):

    nb_profilesargo=len(ds_profiles.N_PROF)
    all_lat=np.zeros((nb_profilesargo))
    all_lon=np.zeros((nb_profilesargo))

    for prof in np.arange(nb_profilesargo):
        lat_prof = ds_profiles.LATITUDE[prof]
        lon_prof = ds_profiles.LONGITUDE[prof]
        all_lat[prof]=lat_prof
        all_lon[prof]=lon_prof
                        
    # location and name of the maskfile of the model configuration
    ds=xr.open_dataset(coordfile)
    lat=ds.nav_lat
    lon=ds.nav_lon
    dsm=xr.open_dataset(maskfile)
    tmask=dsm[nammaskmod]
    dsb=xr.open_dataset(batfile)
    bathy=dsb.Bathymetry
    bathy_mask=np.ma.masked_where(tmask[0,0]==0.,bathy)
    latmin,latmax,lonmin,lonmax=(lat.min(),lat.max(),lon.min(),lon.max())
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)


    fig, axs = plt.subplots(1,2, figsize=(15, 7.5), gridspec_kw={'width_ratios': [4, 1]}, subplot_kw={'projection': ccrs.PlateCarree()})
    axs = axs.ravel()
    axs[0].set_extent((lonmin, lonmax, latmin, latmax))
    pcolor=axs[0].pcolormesh(lon,lat,bathy_mask,transform=ccrs.PlateCarree(),
                             cmap=cmocean.cm.deep,vmin=0,vmax=4000)
    axs[0].coastlines(resolution="10m")
    gl = axs[0].gridlines(draw_labels=True, linestyle=':', color='black',
                          alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    axs[0].tick_params('both',labelsize=22)

    cbar = plt.colorbar(pcolor,orientation='vertical',shrink=0.75,label='m',ax=axs[0])
    axs[0].scatter(all_lon, all_lat, c='r', linewidth=0, s=18);
    axs[0].set_title('There are '+str(len(all_lon))+' argo profiles', size=20);

    textstr = '\n'.join((
                ' simulation = MEDWEST60-BLBT02',
                ' dates = '+str(datemin)+' '+str(datemax),
                ' radius max = '+str(radius_max)+'°',
                ' period = '+str(period)+'d',
                ' depth min = '+str(depthmin)+'m',
                ' nb_profiles = '+str(number_of_model_profiles)))        
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axs[1].text(0.05, 0.95, textstr, transform=axs[1].transAxes, fontsize=14,verticalalignment='top', bbox=props)
    axs[1].axis('off')

    plt.savefig(plotdir+'/map-profiles-'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.png')

def process_one_profile(prof,infos,dirmod,config,case,meshfile,dirargo,radius_max,depthmin,period,number_of_model_profiles,debug_plot=False):

        # Process one profile
        list_profiles = infos.keys()
        reference =  str(list(list_profiles)[prof])
        print('Processing profile ', reference)

        # Get all infos from json file
        lat_prof = infos[list(list_profiles)[prof]]['latitude']
        lon_prof = infos[list(list_profiles)[prof]]['longitude']
        date_prof = infos[list(list_profiles)[prof]]['date']
        file_prof = infos[list(list_profiles)[prof]]['file']
        prof_prof = infos[list(list_profiles)[prof]]['profile no']
        i0 = infos[list(list_profiles)[prof]]['i0']
        j0 = infos[list(list_profiles)[prof]]['j0']

        # List of all model files involved
        date_profmin=datetime.date(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))-datetime.timedelta(days=int(period))
        date_profmax=datetime.date(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))+datetime.timedelta(days=int(period))
        def date_range(start, end):
            r = (end+datetime.timedelta(days=1)-start).days
            return [start+datetime.timedelta(days=i) for i in range(r)]
        dateList = date_range(date_profmin, date_profmax) 
        list_filesmod_T=[]
        list_filesmod_S=[]
        for date in dateList:
            year=date.year
            month=date.month
            day=date.day
            mm="{:02d}".format(month) #month on 2 digits
            dd="{:02d}".format(day) # day on 2 digits
            list_filesmod_T.append(dirmod+'/'+str(config)+'-'+str(case)+'_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridT.nc')
            list_filesmod_S.append(dirmod+'/'+str(config)+'-'+str(case)+'_y'+str(year)+'m'+str(mm)+'d'+str(dd)+'.1h_gridS.nc')
        print(list_filesmod_T)  

        # Read model files
        dsT=xr.open_mfdataset(list_filesmod_T)
        dsS=xr.open_mfdataset(list_filesmod_S)
        gdpts=np.int(np.round(radius_max*60))

        lonmod=dsT.nav_lon[0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
        latmod=dsT.nav_lat[0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
        tempmod=dsT.votemper[:,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
        saltmod=dsS.vosaline[:,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
        depthmod=dsT.deptht
        ds=xr.open_dataset(meshfile)
        maskmod=ds.tmask[0,:,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]
        maskmod0=ds.tmask[0,0,j0-gdpts:j0+gdpts,i0-gdpts:i0+gdpts]

        #Get the depth at every grid point
        d,ly,lx=maskmod.shape
        depthmod2d=np.zeros([lx,ly])
        for j in np.arange(ly):
            for i in np.arange(lx):
                depthmod2d[j,i]=depthmod[np.min(np.where(maskmod[:,j,i].values<1))].values

        # Read the profile file
        tfileargo=dirargo+file_prof
        dsargo=xr.open_dataset(tfileargo)
        latargo=dsargo['LATITUDE'][prof_prof]
        lonargo=dsargo['LONGITUDE'][prof_prof]
        dayargo=dsargo['JULD'][prof_prof]
        tempargo=dsargo['TEMP_ADJUSTED'][prof_prof]
        saltargo=dsargo['PSAL_ADJUSTED'][prof_prof]
        presargo=dsargo['PRES_ADJUSTED'][prof_prof]
        depthargo=seawater.dpth(presargo,latargo)

        # Find the last level and reduce the profiles
        indzprof=np.min(np.where(np.isnan(depthargo)==True))
        dmax=depthargo[indzprof-1]
        obsred_dep=np.zeros(int(indzprof))
        obsred_temp=np.zeros(int(indzprof))
        obsred_salt=np.zeros(int(indzprof))
        for z in np.arange(int(indzprof)):
            obsred_dep[int(z)]=depthargo[int(z)]
            obsred_temp[int(z)]=tempargo[int(z)]
            obsred_salt[int(z)]=saltargo[int(z)]

        # Remove the NaN in the profiles
        indtempnan=np.where(np.isnan(obsred_temp)==True)
        if len(indtempnan[0]) > 0:
            obsred_dep = np.delete(obsred_dep, indtempnan[0])    
            obsred_salt = np.delete(obsred_salt, indtempnan[0])    
            obsred_temp = np.delete(obsred_temp, indtempnan[0])    
        indsaltnan=np.where(np.isnan(obsred_salt)==True)
        if len(indsaltnan[0]) > 0:
            obsred_dep = np.delete(obsred_dep, indsaltnan[0])    
            obsred_salt = np.delete(obsred_salt, indsaltnan[0])    
            obsred_temp = np.delete(obsred_temp, indsaltnan[0])    

        # Find the model profiles
        lon_stacked = lonmod.stack(profile=('x', 'y'))
        lat_stacked = latmod.stack(profile=('x', 'y'))
        mask_stacked = maskmod0.stack(profile=('x', 'y'))
        xr_depthmod2d=xr.DataArray(depthmod2d, dims=('y','x'))
        depth_stacked = xr_depthmod2d.stack(profile=('x', 'y'))

        distance_threshold = radius_max
        square_distance_to_observation = (lon_stacked - lon_prof)**2 + (lat_stacked-lat_prof)**2
        square_distance_to_observation_mask = np.ma.masked_where(mask_stacked==0.,square_distance_to_observation) 
        square_distance_to_observation_sorted = np.sort(square_distance_to_observation_mask)
        nb_profiles_per_timestep=number_of_model_profiles/(24*period*2+24)
        new_threshold=square_distance_to_observation_sorted[int(np.round(nb_profiles_per_timestep)+1)]
        is_closer_to_observation = (square_distance_to_observation < new_threshold) & (depth_stacked > depthmin)

        model_temperature_stacked = tempmod.stack(profile=('x', 'y'))
        model_salinity_stacked = saltmod.stack(profile=('x', 'y'))

        model_temperature_near_observation = model_temperature_stacked.where(is_closer_to_observation,drop=True)
        model_salinity_near_observation = model_salinity_stacked.where(is_closer_to_observation, drop=True)
        lat_near_observation = lat_stacked.where(is_closer_to_observation, drop=True)
        lon_near_observation = lon_stacked.where(is_closer_to_observation, drop=True)

        # Compute statistics on the model profiles
        temp_model_mean = np.mean(model_temperature_near_observation,axis=(0,2))
        temp_percentile_10= np.percentile(model_temperature_near_observation,10,axis=(0,2))
        temp_percentile_90= np.percentile(model_temperature_near_observation,90,axis=(0,2))
        salt_model_mean = np.mean(model_salinity_near_observation,axis=(0,2))
        salt_percentile_10= np.percentile(model_salinity_near_observation,10,axis=(0,2))
        salt_percentile_90= np.percentile(model_salinity_near_observation,90,axis=(0,2))

        # Interpolate on obs vertical grid
        temp_model_mean_depobs=np.interp(obsred_dep,depthmod,temp_model_mean)
        temp_model_percentile_10_depobs=np.interp(obsred_dep,depthmod,temp_percentile_10)
        temp_model_percentile_90_depobs=np.interp(obsred_dep,depthmod,temp_percentile_90)
        salt_model_mean_depobs=np.interp(obsred_dep,depthmod,salt_model_mean)
        salt_model_percentile_10_depobs=np.interp(obsred_dep,depthmod,salt_percentile_10)
        salt_model_percentile_90_depobs=np.interp(obsred_dep,depthmod,salt_percentile_90)

        # Make a debug plot
        if debug_plot == True:
            fig, axs = plt.subplots(1,2, figsize=(10, 6))
            axs = axs.ravel()
            title = 'Temperature and Salinity Profiles for profile '+reference
            plt.suptitle(title,size = 25,y=1.05)
            axs[0].plot(temp_model_mean_depobs,obsred_dep,'b.-', label='temp model')
            axs[0].plot(obsred_temp,obsred_dep,'k.-', label='temp argo')
            axs[0].set_ylabel('Depth [m]', size=14)
            axs[0].set_ylim(2000, 0)
            axs[0].grid(True, which='both')
            axs[0].xaxis.tick_top()
            axs[0].xaxis.set_label_position('top') 
            axs[0].plot(temp_model_percentile_10_depobs,obsred_dep,'b-', label='percent10')
            axs[0].plot(temp_model_percentile_90_depobs,obsred_dep,'b-', label='percent90')
            axs[0].fill_betweenx(obsred_dep, temp_model_percentile_10_depobs, x2=temp_model_percentile_90_depobs, alpha=0.2, facecolor='b')

            axs[1].plot(salt_model_mean_depobs,obsred_dep,'b.-', label='salt model')
            axs[1].plot(obsred_salt,obsred_dep,'k.-', label='salt argo')
            axs[1].set_ylabel('Depth [m]', size=14)
            axs[1].set_ylim(2000, 0)
            axs[1].grid(True, which='both')
            axs[1].xaxis.tick_top()
            axs[1].xaxis.set_label_position('top') 
            axs[1].plot(salt_model_percentile_10_depobs,obsred_dep,'b-', label='percent10')
            axs[1].plot(salt_model_percentile_90_depobs,obsred_dep,'b-', label='percent90')
            axs[1].fill_betweenx(obsred_dep, salt_model_percentile_10_depobs, x2=salt_model_percentile_90_depobs, alpha=0.2, facecolor='b')
            fig.tight_layout()
            plt.savefig(dirplot+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'_prof'+str(prof)+'.png')

        # Write netcdf file
        match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
        debut_ref=match.group(1)
        fin_ref=match.group(2)
        dirname=dirargo+'profiles_files/'+str(config)+'-'+str(case)
        if not os.path.exists(dirname):
            os.makedirs(dirname)    

        outname=dirname+'/'+str(debut_ref)+str(fin_ref)+'_'+str(config)+'-'+str(case)+'_'+str(depthmin)+'m_TS.nc'
        print('output file is '+outname)
        dsout=Dataset(outname,'w')

        today=date.today()
        dsout.description = 'This file contains one profile of temperature and salinity from argo dataset and the mean and 10 and 90 percentile of '+str(config)+'-'+str(case)+' data within a '+str(radius_max)+'deg circle around the location of the profile and '+str(period)+' days before and after it has been sampled. This file has been created '+str(today.day)+'/'+str(today.month)+'/'+str(today.year)

        depth=dsout.createDimension('depth',len(obsred_dep))
        x=dsout.createDimension('x',1)
        y=dsout.createDimension('y',1)

        lat = dsout.createVariable('latitude_profileargo', 'f8', ('y','x'))
        lat.standart_name="latitude_profileargo"
        lat.long_name = "Latitude of selected argo profile" 
        lat.units = "degrees_north"

        lon = dsout.createVariable('longitude_profileargo', 'f8', ('y','x'))
        lon.standart_name="longitude_profileargo"
        lon.long_name = "Longitude of selected argo profile" 
        lon.units = "degrees_east"

        time = dsout.createVariable('time_profileargo', 'f8', ('y','x'))
        time.standart_name="time_profileargo"
        time.timeg_name = "Time in seconds from 1-1-1958 of selected argo profile" 
        time.units = "s"

        depth_argo = dsout.createVariable('depth_argo', 'f8', ('depth'),fill_value=0.)
        depth_argo.units = "m" 
        depth_argo.valid_min = 0.
        depth_argo.valid_max = 8000.
        depth_argo.long_name = "Depth" 

        temp_argo = dsout.createVariable('temp_profileargo', 'f8', ('depth'),fill_value=0.)
        temp_argo.units = "degC" 
        temp_argo.valid_min = -10.
        temp_argo.valid_max = 40.
        temp_argo.long_name = "Temperature profile of the selected argo profile" 

        salt_argo = dsout.createVariable('salt_profileargo', 'f8', ('depth'),fill_value=0.)
        salt_argo.units = "PSU" 
        salt_argo.valid_min = 20.
        salt_argo.valid_max = 40.
        salt_argo.long_name = "Salinity profile of the selected argo profile" 

        mean_temp_model = dsout.createVariable('mean_temp_model', 'f8', ('depth'),fill_value=0.)
        mean_temp_model.units = "degC" 
        mean_temp_model.valid_min = -10.
        mean_temp_model.valid_max = 40.
        mean_temp_model.long_name = "Mean Temperature profile of the model" 

        mean_salt_model = dsout.createVariable('mean_salt_model', 'f8', ('depth'),fill_value=0.)
        mean_salt_model.units = "PSU" 
        mean_salt_model.valid_min = 20.
        mean_salt_model.valid_max = 40.
        mean_salt_model.long_name = "Mean Salinity profile of the model" 

        percent10_temp_model = dsout.createVariable('percent10_temp_model', 'f8', ('depth'),fill_value=0.)
        percent10_temp_model.units = "degC" 
        percent10_temp_model.valid_min = -10.
        percent10_temp_model.valid_max = 40.
        percent10_temp_model.long_name = "Percent 10 Temperature profile of the model" 

        percent10_salt_model = dsout.createVariable('percent10_salt_model', 'f8', ('depth'),fill_value=0.)
        percent10_salt_model.units = "PSU" 
        percent10_salt_model.valid_min = 20.
        percent10_salt_model.valid_max = 40.
        percent10_salt_model.long_name = "Percent 10 Salinity profile of the model" 

        percent90_temp_model = dsout.createVariable('percent90_temp_model', 'f8', ('depth'),fill_value=0.)
        percent90_temp_model.units = "degC" 
        percent90_temp_model.valid_min = -90.
        percent90_temp_model.valid_max = 40.
        percent90_temp_model.long_name = "Percent 90 Temperature profile of the model" 

        percent90_salt_model = dsout.createVariable('percent90_salt_model', 'f8', ('depth'),fill_value=0.)
        percent90_salt_model.units = "PSU" 
        percent90_salt_model.valid_min = 20.
        percent90_salt_model.valid_max = 40.
        percent90_salt_model.long_name = "Percent 90 Salinity profile of the model" 


        lat[:]=lat_prof
        lon[:]=lon_prof
        time[:]=(datetime.datetime(int(date_prof[0:4]),int(date_prof[5:7]),int(date_prof[8:10]))-datetime.datetime(1958,1,1,0,0)).total_seconds()
        depth_argo[:]=obsred_dep
        temp_argo[:]=obsred_temp
        salt_argo[:]=obsred_salt
        mean_temp_model[:]=temp_model_mean_depobs
        mean_salt_model[:]=salt_model_mean_depobs
        percent10_temp_model[:]=temp_model_percentile_10_depobs
        percent10_salt_model[:]=salt_model_percentile_10_depobs
        percent90_temp_model[:]=temp_model_percentile_90_depobs
        percent90_salt_model[:]=salt_model_percentile_90_depobs
        dsout.close()  # close the new file


def process_profiles(config='MEDWEST60',case='BLBT02',member='',dirmod='/mnt/alberta/equipes/IGE/meom/workdir/lerouste/MEDWEST60/',maskfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mask.nc4',hgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_hgr.nc4',zgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_zgr.nc4',batfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4',namlatmod='nav_lat',namlonmod='nav_lon',namdepmod='gdept_1d',nammaskmod='tmask',ymin=2009,mmin=7,dmin=1,ymax=2010,mmax=6,dmax=30,depthmin=0,radius_max=0.25,period=0,number_of_model_profiles=100000,plotdir='plots',jsondir='txt',dmap=1):
    # Loop over all the profiles listed in the json file                   
    datemin=datetime.date(ymin,mmin,dmin)
    datemax=datetime.date(ymax,mmax,dmax)
    jsonfile=jsondir+'/'+str(config)+'-'+str(case)+'_'+str(datemin)+'-'+str(datemax)+'_'+str(depthmin)+'m_'+str(radius_max)+'x'+str(period)+'d_'+str(number_of_model_profiles)+'.json'

    # Read the jsonfile
    sourcefile=open(jsonfile,'rU')
    infos=json.load(sourcefile)
    nb_profilesargo=len(infos)

    # Loop on the number of profiles
    print("Nb de profiles : "+str(nb_profilesargo))
    print(time.strftime('%d/%m/%y %H:%M',time.localtime()))

    for prof in range(nb_profilesargo):
        list_profiles = infos.keys()
        reference =  str(list(list_profiles)[prof])
        match=re.search(r'([\w.-]+).nc([\w.-]+)', reference)
        debut_ref=match.group(1)
        fin_ref=match.group(2)
        dirname=dirargo+'profiles_files/'+str(config)+'-'+str(case)
        outname=dirname+'/'+str(debut_ref)+str(fin_ref)+'_'+str(config)+'-'+str(case)+'_'+str(depthmin)+'m_TS.nc'
        if not os.path.exists(outname):
            process_one_profile(prof,infos,dirmod,config,case,meshfile,dirargo,radius_max,depthmin,period,number_of_model_profiles)
                        
    print(time.strftime('%d/%m/%y %H:%M',time.localtime()))






                
                        
                        


