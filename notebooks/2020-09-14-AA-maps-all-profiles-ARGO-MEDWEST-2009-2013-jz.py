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

import stat_comp_debug as sc
import param_all_profiles_MEDWEST60_2009_jz as param

params = {
    'config':param.config,
    'case':param.case,
    'member':param.member,
    'dirmod':param.dirmod,
    'coordfile':param.coordfile,
    'maskfile':param.maskfile,
    'hgrfile':param.hgrfile,
    'zgrfile':param.zgrfile,
    'batfile':param.batfile,
    'namlatmod':param.namlatmod,
    'namlonmod':param.namlonmod,
    'namdepmod':param.namdepmod,
    'nammaskmod':param.nammaskmod,
    'ymin':param.ymin,
    'mmin':param.mmin,
    'dmin':param.dmin,
    'ymax':param.ymax,
    'mmax':param.mmax,
    'dmax':param.dmax,
    'depthmin':param.depthmin,
    'radius_max':param.radius_max,
    'period':param.period,
    'number_of_model_profiles':param.number_of_model_profiles,
    'plotdir':param.plotdir,
    'ncdir':param.ncdir,
    'dmap':param.dmap,
    'sosie_exec':param.sosie_exec}

datemin=datetime.date(params['ymin'],params['mmin'],params['dmin'])+datetime.timedelta(days=params['period'])
datemax=datetime.date(params['ymax'],params['mmax'],params['dmax'])-datetime.timedelta(days=params['period'])

netcdf_file=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin)+'-'+str(datemax)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'

ds_profiles=xr.open_dataset(netcdf_file)

sc.plot_profiles_argo(ds_profiles,**params)


datemin09=datetime.date(2009,1,1)
datemax09=datetime.date(2009,12,31)
netcdf_file09=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin09)+'-'+str(datemax09)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'
ds_profiles09=xr.open_dataset(netcdf_file09)

datemin10=datetime.date(2010,1,1)
datemax10=datetime.date(2010,12,31)
netcdf_file10=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin10)+'-'+str(datemax10)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'
ds_profiles10=xr.open_dataset(netcdf_file10)

datemin11=datetime.date(2011,1,1)
datemax11=datetime.date(2011,12,31)
netcdf_file11=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin11)+'-'+str(datemax11)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'
ds_profiles11=xr.open_dataset(netcdf_file11)

datemin12=datetime.date(2012,1,1)
datemax12=datetime.date(2012,12,31)
netcdf_file12=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin12)+'-'+str(datemax12)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'
ds_profiles12=xr.open_dataset(netcdf_file12)

datemin13=datetime.date(2013,1,1)
datemax13=datetime.date(2013,12,31)
netcdf_file13=params['ncdir']+'/ARGO_profiles_selection_for_'+str(params['config'])+'-'+str(params['case'])+'_'+str(datemin13)+'-'+str(datemax13)+'_'+str(params['depthmin'])+'m_'+str(params['radius_max'])+'x'+str(params['period'])+'d_'+str(params['number_of_model_profiles'])+'.nc'
ds_profiles13=xr.open_dataset(netcdf_file13)

ds_profiles_all=xr.concat([ds_profiles09,ds_profiles10,ds_profiles11,ds_profiles12,ds_profiles13],dim='N_PROF')
sc.plot_profiles_argo(ds_profiles_all,**params)

