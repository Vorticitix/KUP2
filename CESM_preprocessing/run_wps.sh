#script to start the whole modelling process from WPS, real and WRF.
#!/bin/bash -l
set -ex
module unload xalt

if [[ $# -ne 4 ]]; then
  echo "./start_simulations.sh start_date end_date experiment_name part, e.g., ./start_simulations 2007-11-01 2007-11-08 Peru 1"
  exit
fi

start_date=$1
end_date=$2
exp_name=$3
part=$4

# split the start_date and end_date into single strings for year, month and day
year_s=${start_date:0:4}
mon_s=${start_date:5:2}
day_s=${start_date:8:2}
year_e=${end_date:0:4}
mon_e=${end_date:5:2}
day_e=${end_date:8:2}

    . set_env_wps.sh

    if [ ! -d "$SCRATCH/${exp_name}/WPS" ]; then 
      mkdir -p $SCRATCH/${exp_name}
      mkdir -p $SCRATCH/${exp_name}/WPS
    fi
    if [ ! -d "$SCRATCH/${exp_name}/WPS/${start_date}" ]; then
      cp -r ${path_compiled}/WPSv${wrf_version}_${model_calendar} $SCRATCH/${exp_name}/WPS/${start_date}
    fi
    # get GEOGRID.TBL, METGRID.TBL and Vtable directly from gitlab folder in $PROJECT into WPS
    rm -rf $SCRATCH/${exp_name}/WPS/${start_date}/geogrid/GEOGRID.TBL
    cp $PROJECT/KUP-ESM-AD/WRF/Utilities/GEOGRID.TBL $SCRATCH/${exp_name}/WPS/${start_date}/geogrid/.
    if [[ $model_calendar == "noleap" ]]; then
      rm -rf $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL
      cp $PROJECT/KUP-ESM-AD/WRF/Utilities/METGRID.TBL.CESM $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/.
      mv -f $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL.CESM $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL
      cp $PROJECT/KUP-ESM-AD/WRF/Utilities/Vtable.CESM $SCRATCH/${exp_name}/WPS/${start_date}/ungrib/Variable_Tables/.
    elif [[ $model_calendar == "leap" ]]; then
      rm -rf $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL
      cp $PROJECT/KUP-ESM-AD/WRF/Utilities/METGRID.TBL.ECMWF $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/.
      mv -f $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL.ECMWF $SCRATCH/${exp_name}/WPS/${start_date}/metgrid/METGRID.TBL
    fi

    cd $SCRATCH/${exp_name}/WPS/${start_date}
    # you need to have your own namelist, as it is domain dependent. Check gitlab for example namelist including placeholders for date.
    cp ${path_namelist}/namelist.wps_control_${exp_name} namelist.wps_control
    cp namelist.wps_control namelist.wps
    sed -i 's/YEAR_S/'${year_s}'/g'   namelist.wps
    sed -i 's/MON_S/'${mon_s}'/g'     namelist.wps
    sed -i 's/DAY_S/'${day_s}'/g'     namelist.wps
    sed -i 's/YEAR_E/'${year_e}'/g'   namelist.wps
    sed -i 's/MON_E/'${mon_e}'/g'     namelist.wps
    sed -i 's/DAY_E/'${day_e}'/g'     namelist.wps

    if [ $model_calendar == "leap" ]; then
      ln -sf ungrib/Variable_Tables/Vtable.ERA-interim.pl Vtable
    elif [ $model_calendar == "noleap" ]; then
      ln -sf ungrib/Variable_Tables/Vtable.CESM Vtable
    fi
    #tbd split yearly grib files into months first to increase performance
    if [ $year_s -ne $year_e ]; then
      nr_grib=$(ls -1 ${path_grib}/grib-${year_s}/*${year_s}*.[23]d.grib | wc -l)
      if [ $nr_grib -ne $(ls -1 ../*${year_s}*.[23]d.grib | wc -l) ]; then
        cp ${path_grib}/grib-${year_s}/*${year_s}*grib ../.
      fi 
      nr_grib=$(ls -1 ${path_grib}/grib-${year_e}/*${year_e}*.[23]d.grib | wc -l)
      if [ $nr_grib -ne $(ls -1 ../*${year_e}*.[23]d.grib | wc -l) ]; then
        cp ${path_grib}/grib-${year_e}/*${year_e}*grib ../.
      fi
      if [ ! -e ../${year_s}_${year_e}_2d.grib ] || [ ! -e ../${year_s}_${year_e}_3d.grib ]; then
        cdo mergetime ../*${year_s}*2d.grib ../*${year_e}*2d.grib ../${year_s}_${year_e}_2d.grib
        cdo mergetime ../*${year_s}*3d.grib ../*${year_e}*3d.grib ../${year_s}_${year_e}_3d.grib
      fi
      rm -rf date_to_run.2d.grib
      rm -rf date_to_run.3d.grib
      cdo seldate,${year_s}-${mon_s}-${day_s},${year_e}-${mon_e}-${day_e} ../${year_s}_${year_e}_2d.grib date_to_run.2d.grib
      cdo seldate,${year_s}-${mon_s}-${day_s},${year_e}-${mon_e}-${day_e} ../${year_s}_${year_e}_3d.grib date_to_run.3d.grib
    else
      nr_grib=$(ls -1 ${path_grib}/grib-${year_s}/*${year_s}*.[23]d.grib | wc -l)
      if [ $nr_grib -ne $(ls -1 ../*${year_s}*.[23]d.grib | wc -l) ]; then
        cp ${path_grib}/grib-${year_s}/*${year_s}*grib ../.
      fi
      if [ ! -e ../${year_s}_2d.grib ] || [ ! -e ../${year_s}_3d.grib ]; then
        cdo mergetime ../*${year_s}*2d.grib ../${year_s}_2d.grib
        cdo mergetime ../*${year_s}*3d.grib ../${year_s}_3d.grib
      fi
      rm -rf date_to_run.2d.grib
      rm -rf date_to_run.3d.grib
      cdo seldate,${year_s}-${mon_s}-${day_s},${year_e}-${mon_e}-${day_e} ../${year_s}_2d.grib date_to_run.2d.grib
      cdo seldate,${year_s}-${mon_s}-${day_s},${year_e}-${mon_e}-${day_e} ../${year_s}_3d.grib date_to_run.3d.grib
    fi  
    ./link_grib.csh date_to_run*.grib

    if [[ ! -e ../geo_em.d01.nc ]]; then
      ./geogrid.exe
      cp geo_em* ../.
    fi
    ./ungrib.exe

    if [[ -e ungrib.log ]]; then
      success=$(tail ungrib.log | grep -i Success | wc -l)
      if [[ $success -lt 1 ]]; then
         mailx -s "WPS $year_s $month_s $day_s broken, ungrib didn't run successfully" ${email} < /dev/null
      else
  
        ln -sf ../geo_em* .
        ./metgrid.exe

        if [[ -e metgrid.log ]]; then
          success=$(tail metgrid.log | grep -i Success | wc -l)
          if [[ $success -lt 1 ]]; then
             mailx -s "WPS $year_s $month_s $day_s broken, metgrid didn't run successfully" ${email} < /dev/null
          fi
        fi
      fi
    fi

    mkdir -p $SCRATCH/${exp_name}/jobs
    cd $SCRATCH/${exp_name}/jobs
    if [ ${part} -eq 1 ]; then
      cp $SCRATCH/scripts/* .
    fi
    if [ ! -e namelist.input.in.control ]; then
      cp ${path_namelist}/namelist.input.in.control_${exp_name}  $SCRATCH/${exp_name}/jobs/namelist.input.in.control
    fi
