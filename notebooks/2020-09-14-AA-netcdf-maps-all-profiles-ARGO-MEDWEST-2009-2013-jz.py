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
import param_all_profiles_MEDWEST60_2012_2013_jz as param

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

ds_profiles=sc.selection(**params)
sc.plot_profiles_argo(ds_profiles,**params)


