#!/usr/bin/env python
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
sys.path.insert(0,'/mnt/meom/workdir/alberta/DEV/git/validation-TSprofiles/stat_comp')
sys.path.insert(0,'/mnt/meom/workdir/alberta/DEV/git/validation-TSprofiles/param')

import stat_comp as sc
import param_all_profiles_NANFL60_2009_cal1 as param

## parser and main
def script_parser():
	"""Customized parser.
	"""
	from optparse import OptionParser
	usage = "usage: %prog  --jsonfile name --dir dir"
	parser = OptionParser(usage=usage)
	parser.add_option('--datemin', help="datemin", dest="datemin", type="string", nargs=1)
	parser.add_option('--datemax', help="datemax", dest="datemax", type="string", nargs=1)
	return parser

def main():
	parser = script_parser()
	(options, args) = parser.parse_args()
	optdic=vars(options)

	datemin = optdic['datemin']
	datemax = optdic['datemax']


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
    		'ymin':np.int(datemin[6:]),
    		'mmin':np.int(datemin[3:5]),
    		'dmin':np.int(datemin[0:2]),
    		'ymax':np.int(datemax[6:]),
    		'mmax':np.int(datemax[3:5]),
    		'dmax':np.int(datemax[0:2]),
    		'depthmin':param.depthmin,
    		'radius_max':param.radius_max,
    		'period':param.period,
    		'number_of_model_profiles':param.number_of_model_profiles,
    		'plotdir':param.plotdir,
    		'ncdir':param.ncdir,
    		'dmap':param.dmap,
    		'sosie_exec':param.sosie_exec,
    		'srcargo':param.srcargo}

	ds_profiles=sc.selection(**params)
	#sc.plot_profiles_argo(ds_profiles,**params)

if __name__ == '__main__':
	sys.exit(main() or 0)
