# param_all_profiles_MEDWEST60.py
'''
All the parameters
'''

# simulation to be compared and where to find it
config='MEDWEST60'
case='BLBT02'
member=''
dirmod='/mnt/alberta/equipes/IGE/meom/workdir/lerouste/MEDWEST60/'
coordfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_coordinates_v3.nc4'
maskfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mask.nc4'
hgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_hgr.nc4'
zgrfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_mesh_zgr.nc4'
batfile='/mnt/alberta/equipes/IGE/meom/MODEL_SET/MEDWEST60/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'
namlatmod='nav_lat'
namlonmod='nav_lon'
namdepmod='gdept_1d'
nammaskmod='tmask'
# period considered for the comparison
ymin=2009;mmin=7;dmin=1
ymax=2010;mmax=6;dmax=30
# depth of the desired comparison profile in m
depthmin=0
# radius of the circle around the profile location in which we take the modeled profiles, in °  
radius_max=0.25
# period of time around the profile sampling date in which we take the modeled profiles, in days
period=0
# minimum amount of model profiles to be considered to make a significant statistical comparison, for instance in a 1° square and 30-days window we have 2.6 millions modeled profiles, in a 0.5°x10 days 216 000
number_of_model_profiles=100000
# where to store all the plots and json files
plotdir='/home/alberta/Work/git/validation-TSprofiles/plots'
jsondir='txt'
# where to find the sosie executable
sosie_exec='/mnt/alberta/equipes/IGE/meom/workdir/alberta/DEV/sosie/bin/ij_from_lon_lat.x'
# wether we want to plot the area of comparison for each profile
dmap=0
