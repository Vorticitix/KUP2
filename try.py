#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:05:07 2022

@author: doensen
"""

import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import sys
from scipy.ndimage.filters import gaussian_filter
import scipy
path_cyc = '/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/'
path_T850 = '/storage/climatestor/PleioCEP/doensen/data/extracted/'

def average_2D_matrix(old_array):
    array=old_array.copy()
    #get shape of 2D
    y,x = array.shape
    for i in range(1,x-1):
        for j in range(1,y-1):
            if np.isnan(array[j,i])==False:
                array[j,i]=np.nanmean(array[j-1:j+2,i-1:i+2])
            else:
                array[j,i]=np.nan
    array[0,:]=np.nan;array[y-1,:]=np.nan;array[:,0]=np.nan;array[:,x-1]=np.nan
    return array
#%%
cyc = xr.open_dataset(path_cyc+'all_stats_monmean.nc').squeeze()
cyc = cyc.assign_coords(lon=(((cyc.lon + 180) % 360) - 180)).sortby('lon')
cyc = cyc.where((cyc.time.dt.year>=2000)&(cyc.time.dt.year<=3352),drop=True).squeeze()
cyc = cyc.where((cyc.time.dt.month>=12)|(cyc.time.dt.month<=2),drop=True)
cyc = cyc.drop_vars('time_bnds')
lats_cyc = cyc.lat.values; lons_cyc = cyc.lon.values
T = xr.open_dataset(path_T850+'T850_anom_all.nc').squeeze()
T = T.assign_coords(lon=(((T.lon + 180) % 360) - 180)).sortby('lon')
T = T.where((T.time.dt.year>=2000)&(T.time.dt.year<=3352),drop=True)
T = T.where((T.time.dt.month>=12)|(T.time.dt.month<=2),drop=True)
T = T.sel(lat=lats_cyc,lon=lons_cyc, method='nearest')

cyc_time=np.array([x.dt.year.values + x.dt.month.values/12 - 1/12 for x in cyc.time])
T_time = np.array([x.dt.year.values + x.dt.month.values/12 - 1/12 for x in T.time])
cyc = cyc.assign_coords(time=cyc_time)
T = T.assign_coords(time=T_time)

T = T.reindex(time=cyc.time)
#%%

cyc_rol = cyc.rolling(time=90,min_periods=6).mean(skipna=True)
T_rol = T.rolling(time=90,min_periods=6).mean(skipna=True)

#%%
names = ['Mean Sea Level Pressure',
         'Gradient',
         '100 hPa Geopotential Height',
         'Radius',
         'Depth',
         'Mean Cyclone Related Precipitation',
         'Max Cyclone Related Precipitation',
         'Mean Cyclone Related Wind Speed',
         'Max Cyclone Related Wind Speed']
for i,key in enumerate(list(cyc_rol.keys())):
    lons = cyc[key].lon.values; lats = cyc[key].lat.values
    corr = xr.corr(cyc_rol[key],T_rol.T,dim='time')
    corr = average_2D_matrix(corr.values)
    fig,ax = plt.subplots(figsize=(14,8),subplot_kw={'projection': ccrs.PlateCarree()})
    
    im=ax.pcolormesh(lons,lats,corr,transform=ccrs.PlateCarree(),
                cmap='seismic',vmin=-.5,vmax=.5)
    ax.coastlines()
    fig.colorbar(im,orientation='horizontal')
    fig.suptitle('Correlation between Monthly Average {} and Monthly 850 hPa Temperature Anomaly DJF'\
                 .format(names[i]),fontsize=16)
    fig.subplots_adjust(left=.05,right=.95,top=.95,bottom=.035,hspace=.075)
    fig.savefig(path_cyc+'figs/{}_corr.png'.format(key),dpi=300)
#T_rol = T.rolling(time=360).mean()




