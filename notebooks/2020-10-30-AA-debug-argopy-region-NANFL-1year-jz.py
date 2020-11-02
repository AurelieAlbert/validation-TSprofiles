# Imports of librairies

import sys
sys.path.insert(0,'/gpfswork/rech/egi/rote001/git/argopy')

from argopy import DataFetcher as ArgoDataFetcher

ds_points=ArgoDataFetcher(src='erddap', parallel=True, progress=True).region([-55,-50,30,40,0,10000,'2009-01-01','2009-12-31']).to_xarray()
