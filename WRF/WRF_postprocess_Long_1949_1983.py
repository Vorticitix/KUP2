"""
WRF Postprossing script using WRF python
"""

#%%
#Import modules
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import wrf
from netCDF4 import Dataset
import netCDF4
from glob import glob
import sys
from pathlib import Path
import time
#Set path where WRF output is saved. And set path where to save the postprocessed data
path = '/storage/climatestor/PleioCEP/data/WRF_MedCyclones/Long_1949_1983/'
path_to_save = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/'
files=sorted(glob(path+'*/wrfout*'))
#Set variables to preprocess
vars_2D = ['slp','T2','td2','wspd_wdir10','pw']
#vars_2D = ['wspd_wdir10']
vars_3D = ['QVAPOR','z','ua','va','tk']
fluxs = ['ACGRDFLX','ACHFX','ACLHF']

vars_2D_total = ['slp','T2','td2','wspd_wdir10','pw','ACGRDFLX','ACHFX','ACLHF','prec']


#F_ = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/F.nc')['F']

#%%
for i,file in enumerate(files):
    print(file)
    start = time.time()
    exist_date = file[-19:-9]
    exist_file = 'WRF_{}_3D.nc'.format(exist_date)
    my_file = Path(path_to_save+'Long_1949_1983/'+exist_file)
    if my_file.is_file():
        print(exist_file+' exists')
        continue
    try:
        wrfout = Dataset(file)
    except:
        continue
    # cache = wrf.extract_vars(wrfout, wrf.ALL_TIMES, ("P", "PSFC", "PB", "PH", "PHB",
    #                                        "T", "QVAPOR", "HGT", "U", "V",
    #                                        "W"))
    #wrf.omp_set_num_threads(wrf.omp_get_num_procs())
    # F = F_.squeeze().expand_dims({'Time':len(wrfout['U'])}).transpose('Time','south_north','west_east')
    # wrfout.createVariable('F',F.dtype,('Time','south_north','west_east'))

    # U = wrf.getvar(wrfout,'U',wrf.ALL_TIMES)
    # V = wrf.getvar(wrfout,'V',wrf.ALL_TIMES)
    # T = wrf.getvar(wrfout,'T',wrf.ALL_TIMES)
    # P = wrf.getvar(wrfout,'P',wrf.ALL_TIMES)
    # PB = wrf.getvar(wrfout,'PB',wrf.ALL_TIMES)
    # PRES = P + PB
    # MAPFAC_U = wrf.getvar(wrfout,'MAPFAC_U',wrf.ALL_TIMES)
    # MAPFAC_V = wrf.getvar(wrfout,'MAPFAC_V',wrf.ALL_TIMES)
    # MAPFAC_M = wrf.getvar(wrfout,'MAPFAC_M',wrf.ALL_TIMES)
    # F = F_.squeeze().expand_dims({'Time':len(U)}).transpose('Time','south_north','west_east')
    # pvo = wrf.pvo(ustag=U, vstag=V, theta=T,
    #               pres=PRES, msfu=MAPFAC_U,
    #               msfv=MAPFAC_V, msfm=MAPFAC_M,
    #               dx=20000,dy=20000,cor=F)
    # pvo_interp = wrf.vinterp(wrfout,
    #                          field=pvo,
    #            vert_coord='theta',
    #            interp_levels=[310,320,330,340,350],
    #            extrapolate=True,
    #            field_type="theta",
    #            timeidx=wrf.ALL_TIMES)
    total_2D = []
    # #Extract the convective and non-convective precipitation part and sum them together
    new_dic = {}
    rainc = wrf.getvar(wrfout,'RAINC',wrf.ALL_TIMES)
    rainnc = wrf.getvar(wrfout,'RAINNC',wrf.ALL_TIMES)
    rain_total = rainc + rainnc
    #Precipitation must be de-accumaletd. Subtract the precipitation at time t+1 from the precipitation at time t.
    #This yields hourly precipitation. To compute the hourly precipitation for the last time step, we must open the first time step of the second
    prec = rain_total.copy(deep=True)
    for j in range(len(prec)):
        if j==0:
            if i == 0:
                prec[j,:,:] = rain_total[j,:,:]
            else:
                wrfout2 = Dataset(files[i-1])
                rainc2 = wrf.getvar(wrfout2,'RAINC',wrf.ALL_TIMES)
                rainnc2 = wrf.getvar(wrfout2,'RAINNC',wrf.ALL_TIMES)
                rain_total2 = rainc2 + rainnc2
                prec[j,:,:] = rain_total[j,:,:] - rain_total2[j-1,:,:]
        else:
                prec[j,:,:] = rain_total[j,:,:] - rain_total[j-1,:,:]
    #If path for precipitation does not yet exist, create it
    #Path(path_to_save+'Long_1949_1983/prec/').mkdir(parents=True,exist_ok=True)
    #Define date for file to save. Then save to netcdf file
    init_time = prec.Time[0].dt.strftime("%Y-%m-%d")
    print(str(init_time.values))
    print('prec')
    ds = prec.to_dataset(name='prec')
    total_2D.append(ds)
    
    #ds.to_netcdf(path_to_save+'Long_1949_1983/prec/prec_{}.nc'.format(str(init_time.values)))

    #Energy fluxes need to be deaccumaletd as well. Subtract the energy flux at time t+1 from the energy flux at time t.
    #Divide by 3600 to convert from to W/m^2
    for flux in fluxs:
        #Path(path_to_save+'Long_1949_1983/{}/'.format(flux)).mkdir(parents=True,exist_ok=True)
        da = wrf.getvar(wrfout,flux,wrf.ALL_TIMES)
        da.attrs['projection'] = str(da.attrs['projection'])
        init_time = da.Time[0].dt.strftime("%Y-%m-%d")
        dafl = da.copy(deep=True)
        for j in range(len(prec)):
            if j==0:
                if i == 0:
                    dafl[j,:,:] = da[j,:,:]
                else:
                    wrfout2 = Dataset(files[i-1])
                    da2 = wrf.getvar(wrfout2,flux,wrf.ALL_TIMES)
                    dafl[j,:,:] = da[j,:,:] - da2[j-1,:,:]
            else:
                dafl[j,:,:] = da[j,:,:] - da[j-1,:,:]
        print(flux)
        dafl = dafl/3600
        ds = dafl.to_dataset(name=flux)
        total_2D.append(ds)
        #ds.to_netcdf(path_to_save+'Long_1949_1983/{}/{}_{}.nc'.format(flux,flux,str(init_time.values)))

    # Loop over 3D variables
    dic = {}
    for j in range(len(wrfout['Times'])):
        print(j)
        # cache = wrf.extract_vars(wrfout, timeidx=j,varnames= ("P", "PSFC", "PB", "PH", "PHB",
        #                                    "T", "QVAPOR", "HGT", "U", "V",
        #                                    "W"))
        for var in vars_3D:
            print(var)
            if j == 0:
                dic[var] = []
            #If path with designated variable does not exist, create it
            #Path(path_to_save+'Long_1949_1983/{}/'.format(var)).mkdir(parents=True,exist_ok=True)
            # extract variable using wrf.getvar
            da = wrf.getvar(wrfout,var,timeidx=j)
            #Interpolate variable to pressure levels: 1000,850,700,500,300,200 hPa
            interp = wrf.vinterp(wrfout,
                field=da,
                vert_coord='pres',
                interp_levels=[1000,850,500,300,200],
                extrapolate=True,
                field_type="pres",
                log_p=True)
            #Convert cartopy object to string so it can be saved to netcdf file
            interp.attrs['projection'] = str(da.attrs['projection'])
            #Define initial time of the simulation and save to netcdf file
            #init_time = interp.Time.dt.strftime("%Y-%m-%d")
            dic[var].append(interp)
            #convert to dataset and save to netcdf file
            # ds = interp.to_dataset(name=var)
            # ds.to_netcdf(path_to_save+'Long_1949_1983/{}/{}_{}.nc'.format(var,var,str(init_time.values)))



    # Loop over 2D variables
        for var in vars_2D:
            print(var)
            if j==0:
                dic[var] = []
            #If path with designated variable does not exist, create it
            
            # extract variable using wrf.getvar
            da = wrf.getvar(wrfout,var,timeidx=j)
            #Convert cartopy object to string so it can be saved to netcdf file
            da.attrs['projection'] = str(da.attrs['projection'])
            dic[var].append(da)
            #Define initial time of the simulation and save to netcdf file
            
            #convert to dataset and save to netcdf file
            # ds = da.to_dataset(name=var)
            # ds.to_netcdf(path_to_save+'Long_1949_1983/{}/{}_{}.nc'.format(var,var,str(init_time.values)))
    
    for key in vars_2D:
        print(key)
        da = xr.concat(dic[key],dim='Time')
        ds = da.to_dataset(name=key)
        init_time = da.Time[0].dt.strftime("%Y-%m-%d")
        total_2D.append(ds)
    ds_2d = xr.merge(total_2D)
    ds_2d.to_netcdf(path_to_save+'Long_1949_1983/WRF_{}_2D.nc'.format(str(init_time.values)))
        
    total_3D = []
    for key in vars_3D:
        print(key)
        da = xr.concat(dic[key],dim='Time')
        ds = da.to_dataset(name=key)
        init_time = da.Time[0].dt.strftime("%Y-%m-%d")
        #Path(path_to_save+'Long_1949_1983/{}/'.format(key)).mkdir(parents=True,exist_ok=True)
        total_3D.append(ds)
        #ds.to_netcdf(path_to_save+'Long_1949_1983/{}/{}_{}.nc'.format(key,key,str(init_time.values)))
    ds_3d = xr.merge(total_3D)
    ds_3d.to_netcdf(path_to_save+'Long_1949_1983/WRF_{}_3D.nc'.format(str(init_time.values)))
    
    end = time.time()
    print(end-start)

# %%
