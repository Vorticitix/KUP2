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
#Set path where WRF output is saved. And set path where to save the postprocessed data
path = '/storage/climatestor/PleioCEP/data/WRF_MedCyclones/Long/'
path_to_save = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/'
files=sorted(glob(path+'*/wrfout*'))
#Set variables to preprocess
vars_2D = ['slp','T2','td2','wspd_wdir10','pw']
vars_3D = ['QVAPOR','z','ua','va','tk','avo']
fluxs = ['ACGRDFLX','ACHFX','ACLHF']
#F_ = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/F.nc')['F']

#%%
for i,file in enumerate(files):
    print(file)
    wrfout = Dataset(file)
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
    
    #Extract the convective and non-convective precipitation part and sum them together
    rainc = wrf.getvar(wrfout,'RAINC',wrf.ALL_TIMES)
    rainnc = wrf.getvar(wrfout,'RAINNC',wrf.ALL_TIMES)
    rain_total = rainc + rainnc
    #Precipitation must be de-accumaletd. Subtract the precipitation at time t+1 from the precipitation at time t.
    #This yields hourly precipitation. To compute the hourly precipitation for the last time step, we must open the first time step of the second 
    prec = rain_total.copy(deep=True)
    for j in range(len(prec)):
        if (j == len(prec)-1)&(i<len(files)-1):
            wrfout2 = Dataset(files[i+1])
            rainc2 = wrf.getvar(wrfout2,'RAINC',wrf.ALL_TIMES)
            rainnc2 = wrf.getvar(wrfout2,'RAINNC',wrf.ALL_TIMES)
            rain_total2 = rainc2 + rainnc2
            prec[j,:,:] = rain_total2[0,:,:] - rain_total[j,:,:]
        else:
            prec[j,:,:] = rain_total[j+1,:,:] - rain_total[j,:,:]
    Path(path_to_save+'Long/prec/').mkdir(parents=True,exist_ok=True)
    init_time = prec.Time[0].dt.strftime("%Y-%m-%d")
    print(str(init_time.values))
    print('prec')
    ds = prec.to_dataset(name='prec')
    ds.to_netcdf(path_to_save+'Long/prec/prec_{}.nc'.format(str(init_time.values)))
    for flux in fluxs:
        Path(path_to_save+'Long/{}/'.format(flux)).mkdir(parents=True,exist_ok=True)
        da = wrf.getvar(wrfout,flux,wrf.ALL_TIMES)
        da.attrs['projection'] = str(da.attrs['projection'])
        init_time = da.Time[0].dt.strftime("%Y-%m-%d")
        for j in range(len(da)):
            if j == len(da)-1:
                wrfout2 = Dataset(files[i+1])
                da2 = wrf.getvar(wrfout2,flux,wrf.ALL_TIMES)
                da[j,:,:] = da2[0,:,:] - da[j,:,:]
            else:
                da[j,:,:] = da[j+1,:,:] - da[j,:,:]
        print(flux)
        da = da/3600
        ds = da.to_dataset(name=flux)
        ds.to_netcdf(path_to_save+'Long/{}/{}_{}.nc'.format(flux,flux,str(init_time.values)))
    
    # Loop over 3D variables
    for var in vars_3D:
        #If path with designated variable does not exist, create it
        Path(path_to_save+'Long/{}/'.format(var)).mkdir(parents=True,exist_ok=True)
        # extract variable using wrf.getvar
        da = wrf.getvar(wrfout,var,wrf.ALL_TIMES)
        #Interpolate variable to pressure levels: 1000,850,700,500,300,200 hPa
        interp = wrf.vinterp(wrfout,
               field=da,
               vert_coord='pres',
               interp_levels=[1000,850,700,500,300,200],
               extrapolate=True,
               field_type="pres",
               log_p=True,
               timeidx=wrf.ALL_TIMES)
        #Convert cartopy object to string so it can be saved to netcdf file
        interp.attrs['projection'] = str(da.attrs['projection'])
        #Define initial time of the simulation and save to netcdf file
        init_time = interp.Time[0].dt.strftime("%Y-%m-%d")
        print(var)
        #convert to dataset and save to netcdf file
        ds = interp.to_dataset(name=var)
        ds.to_netcdf(path_to_save+'Long/{}/{}_{}.nc'.format(var,var,str(init_time.values)))
        sys.exit()
        
    # Loop over 2D variables
    for var in vars_2D:
        #If path with designated variable does not exist, create it
        Path(path_to_save+'Long/{}/'.format(var)).mkdir(parents=True,exist_ok=True)
        # extract variable using wrf.getvar
        da = wrf.getvar(wrfout,var,wrf.ALL_TIMES)
        #Convert cartopy object to string so it can be saved to netcdf file
        da.attrs['projection'] = str(da.attrs['projection'])
        #Define initial time of the simulation and save to netcdf file
        init_time = da.Time[0].dt.strftime("%Y-%m-%d")
        print(var)
        #convert to dataset and save to netcdf file
        ds = da.to_dataset(name=var)
        ds.to_netcdf(path_to_save+'Long/{}/{}_{}.nc'.format(var,var,str(init_time.values)))
    sys.exit()
    if i==5:
        sys.exit()
# %%
