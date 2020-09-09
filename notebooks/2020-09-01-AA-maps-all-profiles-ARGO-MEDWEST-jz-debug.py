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

from argopy import DataFetcher as ArgoDataFetcher

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

import sys
sys.path.insert(0,'/gpfswork/rech/egi/rote001/git/validation-TSprofiles/stat_comp')
sys.path.insert(0,'/gpfswork/rech/egi/rote001/git/validation-TSprofiles/param')

import stat_comp as sc
import param_all_profiles_MEDWEST60_1an_jz as param

config=param.config
case=param.case
member=param.member
dirmod=param.dirmod
coordfile=param.coordfile
maskfile=param.maskfile
hgrfile=param.hgrfile
zgrfile=param.zgrfile
batfile=param.batfile
namlatmod=param.namlatmod
namlonmod=param.namlonmod
namdepmod=param.namdepmod
nammaskmod=param.nammaskmod
ymin=param.ymin
mmin=param.mmin
dmin=param.dmin
ymax=param.ymax
mmax=param.mmax
dmax=param.dmax
depthmin=param.depthmin
radius_max=param.radius_max
period=param.period
number_of_model_profiles=param.number_of_model_profiles
plotdir=param.plotdir
ncdir=param.ncdir
dmap=param.dmap
sosie_exec=param.sosie_exec

    datemin=datetime.date(ymin,mmin,dmin)+datetime.timedelta(days=period)
    datemax=datetime.date(ymax,mmax,dmax)-datetime.timedelta(days=period)
    # determining the area
    ds=xr.open_dataset(coordfile)
    lat=ds[namlatmod]
    lon=ds[namlonmod]
    latmodmin,latmodmax,lonmodmin,lonmodmax=(lat.min(),lat.max(),lon.min(),lon.max())
    # use argopy to get the selection of profiles
    ds_points=ArgoDataFetcher(src='erddap').region([lonmodmin.values,lonmodmax.values,latmodmin.values,latmodmax.values,0,12000,str(datemin),str(datemax)]).to_xarray()

    ds_profiles=ds_points.argo.point2profile()
    index_i_model=[]
    index_j_model=[]
    # get rid of profiles not following criteria
    for nprof in range(2)):
        print('Processing profile no '+str(nprof))
        latargo=ds_profiles.LATITUDE[nprof].values
        lonargo=ds_profiles.LONGITUDE[nprof].values
        presargo=ds_profiles.PRES[nprof,:].values
        i0,j0=sc.loc(latargo,lonargo,nprof,hgrfile,namlatmod,namlonmod,nammaskmod,sosie_exec)
        if (i0,j0) == (-1,-1):
            print('profile is not in the domain at all')
        check=sc.check_prof_in_ocean(i0,j0,maskfile,nammaskmod)
        if check == 1:
            print('no, profile is on the land')
        print('yes, profile is in the ocean')
        check=sc.check_close_to_boundaries(nprof,maskfile,nammaskmod,i0,j0,radius_max)
        if check == 1:
            print('no, profile is too close to model boundaries')
        print('yes, profile is not too close to model boundaries')
        check=sc.check_prof_depth(nprof,presargo,latargo,depthmin)
        if check == 1:
            print('no, profile is not deep enough')
        print('yes, profile is deep enough')
        check=sc.check_number_profile(nprof,i0,j0,depthmin,coordfile,maskfile,zgrfile,namlatmod,namlonmod,namdepmod,nammaskmod,lonargo,latargo,radius_max,number_of_model_profiles,period,lonmodmin, lonmodmax,latmodmin, latmodmax,dmap=0)
        if check == 1:
            print('no, there are not enough model profiles')
        print('yes, there are enough model profiles')
        index_i_model.append(i0)
        index_j_model.append(j0)
        ds_one=ds_profiles.sel(N_PROF=nprof)
        ds_prof=ds_one.expand_dims({'N_PROF':1})
        try:
            ds_profiles_out=xr.concat([ds_profiles_out,ds_prof],dim='N_PROF')
        except NameError:
            ds_profiles_out=ds_prof












    nprof=0
        print('Processing profile no '+str(nprof))
        latargo=ds_profiles.LATITUDE[nprof].values
        lonargo=ds_profiles.LONGITUDE[nprof].values
        presargo=ds_profiles.PRES[nprof,:].values
        i0,j0=sc.loc(latargo,lonargo,nprof,hgrfile,namlatmod,namlonmod,nammaskmod,sosie_exec)
        if (i0,j0) == (-1,-1):
            print('profile is not in the domain at all')
        check=sc.check_prof_in_ocean(i0,j0,maskfile,nammaskmod)
        if check == 1:
            print('no, profile is on the land')
        print('yes, profile is in the ocean')
        check=sc.check_close_to_boundaries(nprof,maskfile,nammaskmod,i0,j0,radius_max)
        if check == 1:
            print('no, profile is too close to model boundaries')
        print('yes, profile is not too close to model boundaries')
        check=sc.check_prof_depth(nprof,presargo,latargo,depthmin)
        if check == 1:
            print('no, profile is not deep enough')
        print('yes, profile is deep enough')
        check=sc.check_number_profile(nprof,i0,j0,depthmin,coordfile,maskfile,zgrfile,namlatmod,namlonmod,namdepmod,nammaskmod,lonargo,latargo,radius_max,number_of_model_profiles,period,lonmodmin, lonmodmax,latmodmin, latmodmax,dmap=0)
        if check == 1:
            print('no, there are not enough model profiles')
        print('yes, there are enough model profiles')
        k+=1
        index_i_model.append(i0)
        index_j_model.append(j0)
        ds_one=ds_profiles.sel(N_PROF=nprof)
        ds_prof=ds_one.expand_dims({'N_PROF':k})
        if ds_profiles_out == None:
            ds_profiles_out=ds_prof
        else:
            ds_profiles_out=xr.concat([ds_profiles_out,ds_prof],dim='N_PROF')

