# validation-TSprofiles

In this repo, a statistical comparison of modelled temperature and salinity profiles in the ocean is done between the outputs of a model and the EN4 database.

Install the python libraries :
  - conda env create -f environment.yml
  - source activate stat-comp
  - python -m ipykernel install --user --name stat-comp --display stat-comp

Install the sosie module :
  - git clone https://github.com/brodeau/sosie.git
  - compile with make.macro adapted to the machine : ```ln -sf macro/make.macro_machine make.macro; make bin/ij_from_lon_lat.x```
  - report the path to the bin/ij_from_lon_lat.x executable

  - first, we download the EN4 database relevant to the period of simulation
  - two options are then possible :
     - we run the comparison in notebooks, allowing the parallelization of computation with dask and xarray
     - we run the comparison by submitting scripts to a job scheduler, the parallelization is done by using mpi
  - both options rely on the stat-comp library
