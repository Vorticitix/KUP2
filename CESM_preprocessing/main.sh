#!/bin/bash
set -ex
module unload xalt
cat>copy.sbatch<<EOF
#!/bin/bash -l 
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --account=s881
#SBATCH --constraint=mc

module unload xalt
command="cp -r"
echo -e "\$SLURM_JOB_NAME started on \$(date):\n \$command \$1 \$2\n"
srun -n \$SLURM_NTASKS \$command \$1 \$2
echo -e "\$SLURM_JOB_NAME finished on \$(date)\n"

if [ -n "\$3" ]; then
 # submit job with dependency
  sbatch --dependency=afterok:\$SLURM_JOB_ID \$3
fi
EOF

cp -rf $PROJECT/KUP-ESM-AD/WRF/scripts $SCRATCH/.
cd $SCRATCH/scripts
#you will have to modify the config file according to your needs.
. config.sh_orig
cp config.sh_orig config_${exp_name}.sh

#you will have to modify the set_env_wps according to your needs.
. set_env_wps.sh

part=1
date=$start_date
end_date_sim=$(date -d "${start_date} +${nr_days}days 00:00:00" +"%Y-%m-%d %H:%M:%S")
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
 ./run_wps.sh $date ${end_date_sim:0:10} $exp_name ${part}
 ./wrf.start.sh $part
#./wrf.restart.sh 2 2007-11-02
