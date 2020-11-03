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

import sys
sys.path.insert(0,'/mnt/meom/workdir/alberta/DEV/git/validation-TSprofiles/stat_comp')
sys.path.insert(0,'/mnt/meom/workdir/alberta/DEV/git/validation-TSprofiles/param')

from argopy import DataFetcher as ArgoDataFetcher
import stat_comp as sc
import param_all_profiles_NANFL60_2009_cal1 as param

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
srcargo=param.srcargo 

datemin=datetime.date(ymin,mmin,dmin)+datetime.timedelta(days=period)
datemax=datetime.date(ymax,mmax,dmax)-datetime.timedelta(days=period)
# determining the area
ds=xr.open_dataset(coordfile)
lat=ds[namlatmod]
lon=ds[namlonmod]
latmodmin,latmodmax,lonmodmin,lonmodmax=(lat.values.min(),lat.values.max(),lon.values.min(),lon.max())

print(

ds_points=ArgoDataFetcher(src=srcargo,parallel=True).region([np.float(lonmodmin),np.float(lonmodmax),np.float(latmodmin),np.float(latmodmax),0,10000,str(datemin),str(datemax)]).to_xarray()


