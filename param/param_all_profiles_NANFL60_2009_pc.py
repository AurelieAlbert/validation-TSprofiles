# param_all_profiles_MEDWEST60.py
'''
All the parameters
'''

# simulation to be compared and where to find it
config='NATL60'
case='BLBT02'
member=''
dirmod='/media/extra/DATA/eNATL60'
coordfile='/media/extra/DATA/eNATL60/coordinates_eNATL60.nc'
maskfile='/media/extra/DATA/eNATL60/mesh_mask_eNATL60_3.6.nc'
hgrfile='/media/extra/DATA/eNATL60/mesh_hgr_eNATL60_3.6.nc'
zgrfile='/media/extra/DATA/eNATL60/mesh_zgr_eNATL60_3.6.nc'
batfile='/media/extra/DATA/eNATL60/eNATL60_BATHY_GEBCO_2014_2D_msk_v3.1_lb.nc4'
namlatmod='nav_lat'
namlonmod='nav_lon'
namdepmod='gdept_1d'
nammaskmod='tmask'
namtempmod='votemper'
namsaltmod='vosaline'
# period considered for the comparison
ymin=2009;mmin=1;dmin=1
ymax=2009;mmax=12;dmax=31
# depth of the desired comparison profile in m
depthmin=0
# radius of the circle around the profile location in which we take the modeled profiles, in °  
radius_max=0
# period of time around the profile sampling date in which we take the modeled profiles, in days
period=0
# minimum amount of model profiles to be considered to make a significant statistical comparison, for instance in a 1° square and 30-days window we have 2.6 millions modeled profiles, in a 0.5°x10 days 216 000
number_of_model_profiles=1
# where to store all the plots and json files
plotdir='/home/alberta/Work/git/validation-TSprofiles/plots'
ncdir='/media/extra/DATA/ARGO'
# where to find the sosie executable
sosie_exec='/home/brodeau/DEV/sosie/bin/ij_from_lon_lat.x'
#source for argo fetcher
srcargo='erddap' #'erddap' or 'argovis' 
# wether we want to plot the area of comparison for each profile
dmap=0
