#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Thu Nov 17 15:34:17 2022

@author: doensen
"""

import pandas as pd
import xarray as xr
import numpy as np
import cftime
import sys
from matplotlib import pyplot as plt
from copy import deepcopy
from geopy.distance import great_circle

#%%
path = "/storage/climatestor/PleioCEP/doensen/data/"
file = "fort_36_total_whem_RCP85.txt"

#Open cyclone tracks from CESM
df = pd.read_csv(path+file,header=0,index_col=0)
df = df.where(((df['year']>=3507)&(df['year']<=3601))).dropna()


#%%

#extract values from dataframe into numpy array to severly speed up the process.
iics = df['iic'].values
lats = df['lat'].values
lons = df['lon'].values
years = df['year'].values
months = df['month'].values
#day of month was not present in the dataframe yet so extract it manually
days = np.array([x.zfill(8)[6:8] for x in df['date'].astype(int).astype(str).values]).astype(int)
hours = df['idumi'].values
gzs = df['gz'].values
czs = df['cz'].values
zrads = df['zrad'].values
zdeps = df['zdep'].values
slps = df['slp'].values
precmeans = df['precmean'].values
precmaxs = df['precmax'].values
wsmeans = df['wsmean'].values
wsmaxs = df['wsmax'].values
#Concatenate all arrays into a numpy matrix
total = np.stack([iics,lats,lons,years,months,days,hours,
                  gzs,czs,zrads,zdeps,slps,precmeans,precmaxs,
                  wsmeans,wsmaxs],axis=-1)
#filter erroneous dates
wrong_index = np.where(years<5)[0]
total = np.delete(total,wrong_index,axis=0)
#sort martix by year
total = total[total[:,3].argsort()]
#calcluate index where year variable switches to the next one
jump_years = np.where(np.diff(total[:,3],prepend=total[0,3]))[0]
#calculate matrix length
length = np.arange(len(iics))


#%%
# calculate unique lats and lons
lats_uq = np.unique(lats); lons_uq = np.unique(lons)
#cet variable names and atrributes
new_keys = ['GZ','CZ','ZRAD','ZDEP','PRECMEAN','PRECMAX','WSMEAN','WSMAX','PSL']
new_attrs = [('cyclone gradient','gpm/1000 km'),
             ('1000 hPa geopotential height','gpm'),
             ('cyclone radius','1000 km'),
             ('cyclone depth','gpm'),
             ('mean cylone related precipitation','mm'),
             ('max cyclone related precipitation','mm'),
             ('mean cylone related wind speed','m/s'),
             ('max cylone related wind speed','m/s'),
             ('sea level pressure','hPa')]
#open netcdf file that will act as template to save cyclone data
# file_framework = 'TRANS.3501BP.cam.h1.0005-01-01.sel_alt.nc'
# ds = xr.open_dataset(path+'extracted/'+file_framework)

# for key in list(ds.keys()):
#     if key != 'PSL':
#         ds = ds.drop_vars(key)
        

        
# ds = ds.sel(lat=lats_uq,lon=lons_uq,method='nearest')
# #set all values to nan so it can be replaced later with real values
# ds = ds.where(ds==-99999999999)
# ds['PSL'].attrs['units'] = 'hPa'
# for j,key in enumerate(new_keys):
    
   
#     #ds[key].values = ds['PSL'].values
#     ds[key].attrs['long_name'] = new_attrs[j][0]
#     ds[key].attrs['units'] = new_attrs[j][1]

#loop over matrix rows
for i in length:
    print(i)
    #select variables from matrix for every row
    lat = total[i,1]; lon = total[i,2]
    year = total[i,3]; month = total[i,4]; day = total[i,5]; hour = total[i,6]//10000
    gz = total[i,7]; cz = total[i,8]; zrad = total[i,9]; zdep = total[i,10]
    slp =total[i,11]; precmean  = total[i,12]; precmax = total[i,13]
    wsmean = total[i,14]; wsmax =total[i,15]
    # open new dataset if new year is entered
    if (i==0)|(np.isin(i,jump_years)):
        if i!=0:
            ds.to_netcdf(path+'cyclone_spatial_stats/TRANS.3501BP.RCP85.cam.h1.{:04d}-01-01.cyclone_XL.nc'\
                         .format(int(total[i-1,3])))
        print(year)
        #check whether jump index actually corresponds to a new year
        assert(total[i,3]!=total[i-1,3])
        assert(total[i,3]>=total[i+1,3])
        file_framework = 'TRANS.3501BP.RCP85.cam.h1.{:04d}-01-01.sel_alt.nc'.format(int(year))
        ds = xr.open_dataset(path+'extracted/'+file_framework)
    
        #drop old avriables in CESM output. Only keep sea level pressure    
        for key in list(ds.keys()):
            if key != 'PSL':
                ds = ds.drop_vars(key)
                
    
        #select Europe and Atlantic region only        
        ds = ds.sel(lat=lats_uq,lon=lons_uq,method='nearest')
        #set all values to nan so it can be replaced later with real values
        ds = ds.where(ds==-99999999999)
        #set new variables with corresponding attributes
        ds['PSL'].attrs['units'] = 'hPa'
        for i,key in enumerate(new_keys):
            ds[key] = deepcopy(ds['PSL'])
            ds[key].attrs['long_name'] = new_attrs[i][0]
            ds[key].attrs['units'] = new_attrs[i][1]

    #select indexes that correspond to indexes in the netcdf file
    date = cftime.DatetimeNoLeap(year=year,month=month,day=day,hour=hour,
                                 has_year_zero=True)
    time_idx = np.where(ds.time.values==date)[0][0]
    lat_idx = np.where(np.abs(ds.lat.values-lat)<0.00001)[0][0]
    lon_idx = np.where(np.abs(ds.lon.values-lon)<0.00001)[0][0]
    #assign variables to 
    ds['PSL'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = slp
    ds['GZ'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = gz
    ds['CZ'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = cz
    ds['ZRAD'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = zrad
    ds['ZDEP'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = zdep
    ds['PRECMEAN'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = precmean
    ds['PRECMAX'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = precmax
    ds['WSMEAN'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = wsmean
    ds['WSMAX'].values[time_idx,lat_idx-1:lat_idx+2,lon_idx-1:lon_idx+2] = wsmax    
    


# %%
