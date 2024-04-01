#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 17:17:24 2022

@author: doensen
"""

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import Normalize, BoundaryNorm
import matplotlib.ticker as mticker
from glob import glob
import sys
path = '/storage/climatestor/PleioCEP/doensen/data/count_cyclone/'
#%%
files = sorted(glob(path+'cyclone_count_cesm_????_????_RCP85.nc'))
averages = []
for file in files:
    ds = xr.open_dataset(file)
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    ds = ds.where((ds.dayofyear>=335)|(ds.dayofyear<=59),drop=True)
    ds_yy_mean = ds.sum(dim='dayofyear')
    ds_yy_mean['density'] = ds_yy_mean['density']/90
    weights = np.cos(np.deg2rad(ds.lat))
    weighted = ds_yy_mean.sel(lat=slice(70,30),lon=slice(-65,40)).weighted(weights)
    fldmean = weighted.mean(('lon','lat'))
    averages.append(fldmean)
    print(int(ds.year.min()))
all_fldmean = xr.concat(averages,dim='year')
all_fldmean.to_netcdf(path+'cyclone_density_fldmean_med_all_atl_RCP85.nc')

#%%
ds_all = xr.open_dataset(path+'cyclone_density_fldmean_med_all.nc')
ds_rcp85 = xr.open_dataset(path+'cyclone_density_fldmean_med_all_RCP85.nc').sel(year=slice(3515,3600))
ds_com = xr.concat([ds_all,ds_rcp85],dim='year')
ds_com.to_netcdf(path+'cyclone_density_med.nc')
#%%
fig,ax = plt.subplots()
ds_com = ds_com.sel(year=slice(3000,3600))
ax.plot(ds_com['year']-1502,ds_com.density)
ax.plot(ds_com['year']-1502,ds_com.density.rolling(year=30,center=True).mean())
ax.grid()


    