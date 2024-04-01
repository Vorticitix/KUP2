#!/bin/bash -l
# load modules
module unload PrgEnv-gnu
module unload cray-netcdf
module unload cray-hdf5
module load PrgEnv-intel
module load cray-netcdf-hdf5parallel

. config_${exp_name}.sh
# set path to files
# path to the compiled version of WRF
export path_compiled="/project/s1171/odoensen"
# path to template namelist.wps_control and namelist.input.in.control
export path_namelist="/project/s1171/odoensen/WRF_run_utilities"
#path to grib files (should contain the year in their file name)
export path_grib="$SCRATCH/CESM-input/${exp_name}"
#export path_grib="/project/s905/santosgo/DATA/ERA5_WW/"
