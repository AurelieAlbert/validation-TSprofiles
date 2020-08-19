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
import param_1000m as param

    config=param.config
    case=param.case
    member=param.member
    dirmod=param.dirmod
    meshfile=param.meshfile
    batfile=param.batfile
    namlatmod=param.namlatmod
    namlonmod=param.namlonmod
    namdepmod=param.namdepmod
    nammaskmod=param.nammaskmod
    ymin=param.ymin;mmin=param.mmin;dmin=param.dmin
    ymax=param.ymax;mmax=param.mmax;dmax=param.dmax
    depthmin=param.depthmin
    radius_max=param.radius_max
    period=param.period
    number_of_model_profiles=param.number_of_model_profiles
    plotdir=param.plotdir
    jsondir=param.jsondir
    dmap=param.dmap
    datemin=datetime.date(ymin,mmin,dmin)+datetime.timedelta(days=period)
    datemax=datetime.date(ymax,mmax,dmax)-datetime.timedelta(days=period)
    # determining the area
    ds=xr.open_dataset(meshfile)
    lat=ds[namlatmod]
    lon=ds[namlonmod]
    latmodmin,latmodmax,lonmodmin,lonmodmax=(lat.min(),lat.max(),lon.min(),lon.max())
    # use argopy to get the selection of profiles
    ds_points=ArgoDataFetcher().region([lonmodmin,lonmodmax,latmodmin,latmodmax,0,12000,str(datemin),str(datemax)]).to_xarray()
    ds_profiles=ds_points.argo.point2profile()
    #HTTPError on jean zay => test that part on pc, download meshfile

