#!/bin/bash
#-- Bern 15.10.2019, Bern, 1.12.2020
#-- Author: Emmanuele Russo, M.Messmer
#-- Script for postprocessing the WRF outputs for present days. 
#   In this case we consider T2M from 4 simulations. Monthly values of Daily means are calculated.
#   We use the CSCS multitasks software greasy, for running one job for each cpu of a node.

#-- Variables and Paths Declaration
. config_${exp_name}.sh
#period of the current simulation corresponding to specific directory name
start_date=$1
part=$2
indir=$SCRATCH/${exp_name}/WRF/${start_date}
outdir=$SCRATCH/${exp_name}/postproc
list_var=(${var_list[@]})

# number of domains of the current simulation
ndom=$ndoms 

cat>prepare_post_${part}.sh<<EOF
#!/bin/bash
set -ex

. config_${exp_name}.sh
list_var=(${var_list[@]})

mkdir -p $SCRATCH/${exp_name}/scripts
scriptdir=$SCRATCH/${exp_name}/scripts
mkdir -p $outdir

#-- Remove existing tasks
rm -f \$scriptdir/task*.sh
rm -f \$scriptdir/tasks*.txt

#-- Enter Main Directory of experimens
cd $SCRATCH/${exp_name}/WRF/${start_date}

#-- Produce the tasks scripts with each CDO task that we want to perform
itask=0 # Number of current task to execute

for var_name in "\${list_var[@]}";do

   mkdir -p ${outdir}/${name_var}
   test_num=\$(ls -1 wrfout* | wc -l)
   if ((\$test_num==${ndoms})); then
       
      dom=1
      for name_out in wrfout*; do
         itask=\$(( itask+1 ))              
         if [[ \$itask -gt 36 ]]; then
            echo "Too many tasks (itask="\$itask") are requested. Only 36 tasks are available in a single node."
            exit 1
         fi

#-- Create Task script to be run by greasy
cat>\${scriptdir}/task.\${itask}.sh<<EOF2
#!/bin/bash
module load CDO
mkdir -p ${outdir}/\${var_name}
cdo selvar,\${var_name} ${indir}/\${name_out} ${outdir}/\${var_name}/\${var_name}_D0\${dom}_${start_date}.nc 
EOF2

         chmod +x \${scriptdir}/task.\${itask}.sh
         echo \${scriptdir}/task.\${itask}.sh >> \${scriptdir}/tasks.txt

         dom=\$((dom+1))

      done
   fi
done

# prepare greasy to compress wrfout files
ntask=0 # Number of current task to execute

for (( name_dom=1; name_dom<=${ndom}; name_dom++ ));do

   #check if number of output files is equal to 4
   test_num=\$(ls -1 wrfout*00 | wc -l)
   if ((\$test_num==${ndom})); then

      # the 00 at the end is to avoid repeating the operation is gzip files are already present
      for name_out in wrfout_d0\${name_dom}*00; do
         ntask=\$(( ntask+1 ))
         if [[ \$ntask -gt 36 ]]; then
            echo "Too many tasks (ntask="\$ntask") are requested. Only 36 tasks are available in a single node."
            exit 1
         fi

#-- Create Task script to be run by greasy
cat>\${scriptdir}/task_compress_\${ntask}.sh<<EOF3
#!/bin/bash
pigz ${indir}/\${name_out}
EOF3

         chmod +x \$scriptdir/task_compress_\${ntask}.sh

         echo \$scriptdir/task_compress_\${ntask}.sh >> \${scriptdir}/tasks_compress.txt

      done
   fi
##   test_num=\$(ls -1 wrfrst*00 | wc -l)
##   if ((\$test_num==${ndom})); then
##
##      # the 00 at the end is to avoid repeating the operation is gzip files are already present
##      for name_out in wrfrst_d0\${name_dom}*00; do
##         ntask=\$(( ntask+1 ))
##         if [[ \$ntask -gt 36 ]]; then
##            echo "Too many tasks (ntask="\$ntask") are requested. Only 36 tasks are available in a single node."
##            exit 1
##         fi
##
###-- Create Task script to be run by greasy
##cat>\${scriptdir}/task_compress_\${ntask}.sh<<EOF4
###!/bin/bash
##pigz ${indir}/\${name_out}
##EOF4
##
##         chmod +x \$scriptdir/task_compress_\${ntask}.sh
##
##         echo \$scriptdir/task_compress_\${ntask}.sh >> \${scriptdir}/tasks_compress.txt
##
##      done
##   fi
done
cd -
rm -rf export_tasks.sh
echo \$itask \$ntask
cat>export_tasks.sh<<EOF5
#!/bin/bash
export itask=\$itask
export ntask=\$ntask
export scriptdir=\$scriptdir
EOF5
EOF

cat>wrf_post_${part}.sbatch<<EOF
#!/bin/bash -l
#SBATCH --out=O.POST.${part}
#SBATCH --error=E.POST.${part}
#SBATCH --job-name=Post-${exp_name}-part${part}
#SBATCH --time=${post_time}
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
##SBATCH --ntasks-per-node=\${itask}
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$email
#SBATCH --account=$account

set -ex
./prepare_post_${part}.sh
. export_tasks.sh

module load GREASY

export CRAY_CUDA_MPS=1
export CUDA_VISIBLE_DEVICES=0
export GPU_DEVICE_ORDINAL=0
 
greasy \${scriptdir}/tasks.txt

rm -rf prepare_post_$((part-1)).sh
rm -rf wrf_post_$((part-1)).sbatch
EOF

cat>wrf_comp_${part}.sbatch<<EOF
#!/bin/bash -l
#SBATCH --out=O.COMP.${part}
#SBATCH --error=E.COMP.${part}
#SBATCH --job-name=Comp-${exp_name}-part${part}
#SBATCH --time=${post_time}
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
##SBATCH --ntasks-per-node=\${ntask}
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --mail-type=$mail_type
#SBATCH --mail-user=$email
#SBATCH --account=$account

set -ex
. export_tasks.sh
module load GREASY

export CRAY_CUDA_MPS=1
export CUDA_VISIBLE_DEVICES=0
export GPU_DEVICE_ORDINAL=0

greasy \${scriptdir}/tasks_compress.txt

rm -rf wrf_comp_$((part-1)).sbatch
EOF

chmod +x prepare_post_${part}.sh

