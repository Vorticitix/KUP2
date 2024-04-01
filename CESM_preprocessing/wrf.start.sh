#!/bin/bash
set -ex
# loding the config file again, to export all the needed user specific variables
. config_${exp_name}.sh

part=$1
nd=$nr_days

# calculating the start and end date of the current simulation chunk
start_date=$(date -d "${start_date:0:10} 00:00:00" +"%Y-%m-%d %H:%M:%S")
end_date_sim=$(date -d "${start_date:0:10} +${nd}days 00:00:00" +"%Y-%m-%d %H:%M:%S")
if [[ $model_calendar == "noleap" ]]; then
  # split date string into three parts separated by "-", i.e., year , month and day
  # assign the three splits to new variables
  year_s=${start_date:0:4}
  mon_s=${start_date:5:2}
  day_s=${start_date:8:2}
  # use only month and day and the year 1950, because this year is for sure a noleap year. So we are substituting every year with this noleap year to avoid leap years.
  date_new=1950-${mon_s}-${day_s}
  diff=$((year_s - 1950))
  end_date_int=$(date -d "${date_new} +${nr_days}days 00:00:00" +"%Y-%m-%d %H:%M:%S")
  year_e=${end_date_int:0:4}
  mon_e=${end_date_int:5:2}
  day_e=${end_date_int:8:2}
  end_date_sim=$((year_e + diff))-${mon_e}-${day_e}
fi

year_s=${start_date:0:4}
mon_s=${start_date:5:2}
day_s=${start_date:8:2}
year_e=${end_date_sim:0:4}
mon_e=${end_date_sim:5:2}
day_e=${end_date_sim:8:2}

cd $SCRATCH/${exp_name}/jobs

#write file to submit real.exe
cat>run_real.sh<<EOF
#!/bin/bash
#SBATCH --output=O.REAL.${part}
#SBATCH --error=E.REAL.${part}
#SBATCH --job-name=real-p${part}-${exp_name}
#SBATCH --nodes=${nnodes_real}
#SBATCH --ntasks-per-node=${nmpi_tasks}
#SBATCH --partition=${partition}
#SBATCH --constrain=mc
#SBATCH --time=${real_time}
#SBATCH --mail-type=${mail_type}
#SBATCH --mail-user=${email}
#SBATCH --account=${account}

set -ex
# load the needed modules
. load_modules.sh
# check if folder is available and if executables are present
mkdir -p $SCRATCH/${exp_name}/WRF
cp -rf $SCRATCH/WRF_executables/${exp} $SCRATCH/${exp_name}/WRF/${start_date:0:10}
cd $SCRATCH/${exp_name}/WRF/${start_date:0:10}/
cp -rf $SCRATCH/${exp_name}/jobs/namelist.input.in.control namelist.input

#set restart intervall in a way that we only have an output at the end of the simulation (can be adapted if needed)
restart_interval=$((nd*24*60))
sed -i 's/SIMDAY01/'${nd}'/g'    namelist.input
sed -i 's/YEAR1/'${year_s}'/g'   namelist.input 
sed -i 's/MONTH1/'${mon_s}'/g'   namelist.input 
sed -i 's/DAY1/'${day_s}'/g'     namelist.input 
sed -i 's/HOUR1/'00'/g'          namelist.input
sed -i 's/YEAR2/'${year_e}'/g'   namelist.input 
sed -i 's/MONTH2/'${mon_e}'/g'   namelist.input 
sed -i 's/DAY2/'${day_e}'/g'     namelist.input 
sed -i 's/BOOLEAN/'false'/g'     namelist.input
sed -i 's/RESTART_INT/'\${restart_interval}'/g'     namelist.input

ln -sf $SCRATCH/${exp_name}/WPS/${start_date:0:10}/met_em.d0*.* .

srun -n ${nmpi_tasks} -C mc ./real.exe

EOF


#create a script that submits wrf.exe in dependence of real.exe
cat>submit_wrf.sh<<EOF
#!/bin/bash
set -ex
cd $SCRATCH/${exp_name}/jobs

# export the job_id of the real.exe submission (variable kkk), so that a dependency can be created
. export_real.sh
sbatch --dependency=afterok:\$kkk run_wrf.sh

kk=\`squeue -hn WRF-p${part}-${exp_name} -u ${user_name} | awk '{print substr(\$0,1,8);exit}'\`

# save the wrf job_id into a file
rm -rf export_job.sh
cat>export_job.sh<<EOF3
#!/bin/bash
export kk=\$kk
EOF3

# get varible list from config.sh
list_var=(`echo "${var_list[@]}"`)
#submit the wrf variable extraction as soon as wrf.exe is successfully done
sbatch --dependency=afterok:\$kk wrf_post_${part}.sbatch
p1=\`squeue -hn Post-${exp_name}-part${part} -u ${user_name} -t Pending | awk '{print substr(\$0,1,8);exit}'\`
sbatch --dependency=afterok:\$p1 wrf_comp_${part}.sbatch
p2=\`squeue -hn Comp-${exp_name}-part${part} -u ${user_name} -t Pending | awk '{print substr(\$0,1,8);exit}'\`
module unload xalt
sbatch --dependency=afterok:\$p2 --job-name=copy_wrfout $SCRATCH/${exp_name}/jobs/copy.sbatch $SCRATCH/${exp_name}/WRF/${start_date:0:10}/wrfout*.gz  $path_project_wrfout
sbatch --dependency=afterok:\$p2 --job-name=copy_wrfrst $SCRATCH/${exp_name}/jobs/copy.sbatch $SCRATCH/${exp_name}/WRF/${start_date:0:10}/wrfrst*.gz  $path_project_wrfrst
for var_name in "\${list_var[@]}";do
  module unload xalt
  sbatch --dependency=afterok:\$p1 --job-name=copy_\${var_name} $SCRATCH/${exp_name}/jobs/copy.sbatch $SCRATCH/${exp_name}/postproc/\${var_name}/\${var_name}_D*_${start_date:0:10}.nc $path_project_post
done

EOF

#write file to submit wrf.exe
cat>run_wrf.sh<<EOF
#!/bin/bash
#SBATCH --output=O.WRF.${part}
#SBATCH --error=E.WRF.${part}
#SBATCH --job-name=WRF-p${part}-${exp_name}
#SBATCH --nodes=${nnodes_wrf}       
#SBATCH --ntasks-per-node=${nmpi_tasks}
#SBATCH --partition=${partition}
#SBATCH --constrain=mc
#SBATCH --time=${wrf_time}
#SBATCH --mail-type=${mail_type}
#SBATCH --mail-user=${email}
#SBATCH --account=${account}

set -ex
# if you want  to load and unload other modules, please adapt the load_modules.sh
. load_modules.sh

cd $SCRATCH/${exp_name}/WRF/${start_date:0:10}

rm -rf ../../jobs/simday.sh
cat>../../jobs/simday.sh<<EOF2
#!/bin/bash
export sday=${nd}
EOF2

# initialize next chunk
cd $SCRATCH/${exp_name}/jobs
./wrf.restart.sh $((part+1)) ${end_date_sim:0:10}

cd $SCRATCH/${exp_name}/WRF/${start_date:0:10}
tasks=$((nnodes_wrf * nmpi_tasks))
# execute current job
srun -n \${tasks} -C mc ./wrf.exe  

EOF

mkdir -p $SCRATCH/WRF_executables
# check if executables are present, if not copy again from project and submit real in dependence to it, otherwise just submit real.exe
if [ ! -d $SCRATCH/WRF_executables/${exp} ] || [ ! -f $SCRATCH/WRF_executables/${exp}/wrf.exe ] || [ ! -f $SCRATCH/WRF_executables/${exp}/real.exe ]; then
  cp -rf ${path_wrf_executables} $SCRATCH/WRF_executables/${exp}
  rm -rf $SCRATCH/WRF_executables/${exp}/*.exe
  cp -rf ${path_wrf_executables}../main/*.exe $SCRATCH/WRF_executables/${exp}/.
#  ff=`squeue -hn copy_exefiles -u ${user_name} | awk '{print substr($0,1,8);exit}'`
#  sbatch --dependency=afterok:$ff run_real.sh
  sbatch run_real.sh
else 
  sbatch run_real.sh
fi
kkk=`squeue -hn real-p${part}-${exp_name} -u ${user_name} | awk '{print substr($0,1,8);exit}'`
# save real.exe job_id to file
rm -rf export_real.sh
cat>export_real.sh<<EOF3
#!/bin/bash
export kkk=$kkk
EOF3

# start postprocessing and compressing script 
./post_greasy.sh ${start_date:0:10} ${part}

chmod +x submit_wrf.sh
./submit_wrf.sh

