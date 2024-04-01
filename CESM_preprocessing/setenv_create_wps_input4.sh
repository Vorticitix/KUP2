#!/bin/sh

#-- Libraries Upload
#module unload cray-netcdf-hdf5parallel
#module load CDO
#module load NCL
#module load NCO

module load NCO/4.8.1-foss-2020a


set -ex

export CESM_PATH="/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/"  #set path to the folder of CESM input
export FILE_CAM="TRANS.3501BP.cam.h1."                        #set first part of filename for CAM
export FILE_CLM="TRANS.3501BP.clm2.h1."                       #set first part of filename for CLM
export scratchdir="/storage/climatestor/PleioCEP/doensen/CESM-input"
# To do:
#pull pressure interpolation ncl scripts from github to $scratchdir
#adjust ncl files

# adjust the arguments according to your needs
#$1 = start_year
#$2 = end_year
#$3 = off_set needed to transform unrealistic years into regular years, e.g., 1990
#$4 = exp_name used to create the folder for the corresponding experiment
bash run_cesm_preprocessing.sh 3443 3472 -1502 Test
