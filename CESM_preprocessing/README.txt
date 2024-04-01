###########################################################
# 
#
###########################################################

E. Russo and M. Messmer, 24th Nov 2020

This folder will be copied to scratch before running WRF. This folder must first
be cloned into your project (i.e., cd $PROJECT; git clone
http://gitlab.climate.unibe.ch/messmerm/KUP-ESM-AD.git; enter user name and
password)

Please make sure to adjust the files in jour project folder according to your 
needs.

#### Preprocessing and WPS ####

## for leap years, i.e., reanalysis data (ERA-5)


## for no-leap years, i.e. climate model simulations (CESM)

# setenv_create_wps_input.sh
Can be run to initiate all the environmental variables and to run 
run_cesm_preprocessing.sh

# run_cesm_preprocessing.sh 
This script prepares the CESM output to be digested by the WPS program
of WRF. Only needed in combination with a climate model output.
Can be used with the simulations from Jonathan, but not Sandro Blumer.

# set_env_wps.sh
This scripts sets the environment to run the WRF preprocessing system. You will
have to adapte the paths according to your needs

# config.sh_orig
This file exports  all the user dependent values (i.e., similar to a namelist),
so please make sure to adjust it according to your needs.
This file is copied to config_${exp_name}.sh, which will be copied into the according
scratch folder, so that various different simulations can be run at the same time

# run_wps.sh
This script copies the compiled WPS folder to scratch and overwrites the 
GEOGRID.TBL, METGRID.TBL and Vtable used for leap or noleap simulations from the
gitlab folder on your $PROJECT, so that this is not forgotten. Then it starts
to run geogrid.exe in case it is the first loop and then it starts ungrib and 
metgrid for the current chunk.

# wrf.start.sh, wrf.restart.sh
The main difference between the wrf.start.sh and wrf.restart.sh is to set the
boolean researt in the namelist.input.in.control to true after the first
iteration.
These two scripts calculate the end time of the chunk according to user 
specified duration in the config.sh. It creates the run file for real.exe (which
preciously starts the run_wps.sh in the wrf.restart.sh) run_real.sh.
Additionally the submit_wrf.sh and run_wrf.sh are created.

submit_wrf.sh starts run_wrf.sh in dependence of the real.exe job and submits
the postprocessing and compressing scripts written by post_greasy.sh.

run_wrf.sh initiates the next round of chunks, by running wrf.restart.sh, and it
starts the current wrf.exe file. 

the script run_real is submitted and the post_greasy.sh is started as well as
the submit_wrf.sh script

# post_greasy.sh
It creates the three files prepare_post.sh, wrf_post.sbatch, wrf_comp.sbatch. 
The prepare_post.sh file prepares the greasy tasks for the variable extraction
and the compression of wrfout files.
The wrf_post and wrf_comp submit the greasy jobs. 

# main.sh
This script must be copied to $SCRATCH. It can then be run to start the whole 
simulation chain. If the config.sh_orig and the set_env_wps.sh are setup properly
you can start your simulation by submitting ./main.sh in the terminal.

#!!!! Attention !!!!
This folder structure and the fact that WPS and WRF folders are created with all
the rsl.out and rsl.error files the file limit imposed on scratch can be reached
quite quickly. So you have to keep and eye on it and delete rsl files from WRF
parts that have been successfully finished and you can also delete WPS folders. 
These can be created quickly again if needed.
Maybe we can think of including some lines in the script to delete these fodlers
and files.

WRF and postporcessing are still missing
