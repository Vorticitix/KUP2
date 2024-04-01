#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from my_tools import *
import xarray as xr
from geopy.distance import great_circle
import matplotlib.ticker as mticker
from glob import glob
from cftime import DatetimeNoLeap
import sys
import time 

def haversine(ref_lon, ref_lat, lon, lat):
    """
    Vectorized Haversine function to calculate distance between two sets of lon-lat coordinates.
    """
    ref_lon, ref_lat, lon, lat = map(np.radians, [ref_lon, ref_lat, lon, lat])
    dlon = lon - ref_lon
    dlat = lat - ref_lat
    a = np.sin(dlat/2.0)**2 + np.cos(ref_lat) * np.cos(lat) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of Earth in kilometers
    return c * r
# %%
print('hallo')
#set path and filename to analyze 
path = "/storage/climatestor/PleioCEP/doensen/data/"
files = sorted(glob(path+'cyclone_1805_2404/fort_30_????_????'))
filez_to_write = 'count_cyclone/cyclone_count_cesm_{:04d}_{:04d}.nc'

#List all necessary cyclone files

#

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep']
for file in files:
    t0=time.time()
    print(file)
#Open 10 year cyclone script
    df = pd.read_csv(file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
    #Convert coordinate indexes to real  cartesian coordinates
    df['lat']=ilat_to_lat(df['ilat'])
    df['lon']=ilon_to_lon(df['ilon'])
    df['lon']= (df['lon'] + 180) % 360 - 180
    #df = df.where((df['lon']>=200)|(df['lon']<=60)).dropna().reset_index()
    #Exclude all dates that are lower than year 5
    df = df.where(df['date']>50000).dropna().reset_index()
    #Extract year, month and dayofyear from cyclone data
    df['month'] = [int(x[4:6]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['year'] = [int(x[:4]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['day'] = [int(x[6:8]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['hour'] = [x for x in (df['idumi'].astype(int)//10000).astype(str).str.zfill(2)]
    df = df.where(~((df['month']==2)&(df['day']==29))).dropna()
    df['zrad'] = df['zrad'].values*1.5
    #print(df)
    #df = df.where((df['year']>=start_year)&(df['year']<=start_year + dy -1)).dropna()
    # df['dayofyear'] = [pd.to_datetime('2022-{}-{}'.format(x[4:6],x[6:8])).dayofyear\
    #                    for x in df['date'].astype(int).astype(str).str.zfill(8)]
    # df['dayofyear'] = df['dayofyear'] + df['idumi']/240000
    df['datetime'] = [df['year'].iloc[i].astype(str).zfill(4)+df['month'].iloc[i].astype(str).zfill(2)\
                      +df['day'].iloc[i].astype(str).zfill(2)+df['hour'].iloc[i] for i in range(len(df['year']))]

    #Create arrays of all available cartesian coordinates and unique years
    lats_uq = df['lat'].sort_values(ascending=False).dropna().unique()
    lons_uq_sorted = df['lon'].sort_values().dropna().unique()
    lons_uq = np.append(lons_uq_sorted[lons_uq_sorted>60]
                        ,lons_uq_sorted[lons_uq_sorted<=60])
    lats_idx = np.arange(len(lats_uq))
    lons_idx = np.arange(len(lons_uq))
    range_years =np.unique(df['year'].values)
    #Create empty array with shape of 
    zeros = np.full((len(range_years)*1460,len(lats_uq),len(lons_uq)),False)
    #   zeros_frac = np.zeros((len(lats_uq),len(lons_uq),len(range_years),1460))
    #zeros_2 = np.zeros((len(lats_uq),len(lons_uq),12))\
    
    lons_grid,lats_grid=np.meshgrid(lons_uq,lats_uq)
    dates = np.array([])
    for year in range_years:
        start = DatetimeNoLeap(year,1,1,0)
        end = DatetimeNoLeap(year,12,31,18)
        date_range = xr.cftime_range(start,end,freq='6H',calendar='noleap')
        dates = np.append(dates,date_range.values)
    #range_dayofyears = np.unique(df['dayofyear'])
    for k,datetime in enumerate(np.unique(df['datetime'])):
        df_datetime = df.where(df['datetime']==datetime).dropna(how='all')
        # year = df_datetime['year'].iloc[0]
        # month = df_datetime['month'].iloc[0]
        # day = df_datetime['day'].iloc[0]
        # hour = int(df_datetime['hour'].iloc[0])
        # date = DatetimeNoLeap(year,month,day,hour)
        print(k) 
        for i,latlon in df_datetime[['lat','lon','zrad']].iterrows():
    
            lat_fix,lon_fix,zrad = latlon.values
            dist = haversine(lon_fix, lat_fix, lons_grid, lats_grid)
            dist = dist < zrad*1000
            zeros[k,:,:] =  zeros[k,:,:] + dist
            
            # lat_idx = np.where(lat_fix==lats_uq)[0][0]
            # lon_idx = np.where(lon_fix==lons_uq)[0][0]
            # lon_idx_min = lon_idx - 14 ; lon_idx_max = lon_idx + 16
            # lat_idx_min = lat_idx - 6 ; lat_idx_max = lat_idx + 7
            # if lon_idx_min<0:
            #     lon_idx_min=len(lons_uq) - lon_idx_min
            # elif lon_idx_max>len(lons_uq):
            #     lon_idx_max=lon_idx_max - len(lons_uq)
            # if lat_idx_min<0:
            #     lat_idx_min = 0
            # elif lat_idx_max>len(lats_uq):
            #     lat_idx_max=len(lats_uq)
            # if lon_idx_min > lon_idx_max:
            #     dist_1_left = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
            #         for lon in lons_uq[lon_idx_min:len(lons_uq)]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
            #     dist_1_right = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
            #         for lon in lons_uq[0:lon_idx_max]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
            #     dist_1_left_bool = (dist_1_left<=zrad*1000).astype(int)
            #     dist_1_right_bool = (dist_1_right<=zrad*1000).astype(int)
            #     zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:len(lons_uq),j,k] = np.add(dist_1_left_bool,zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:len(lons_uq),j,k])
            #     zeros_1[lat_idx_min:lat_idx_max,0:lon_idx_max,j,k] = np.add(dist_1_right_bool,zeros_1[lat_idx_min:lat_idx_max,0:lon_idx_max,j,k])
            # else:
            #     dist_1 = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
            #         for lon in lons_uq[lon_idx_min:lon_idx_max]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
            
    
            #     dist_1_bool = (dist_1<=zrad*1000).astype(int)
    
            #     zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,k] = np.add(dist_1_bool,zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,k])
    
        #zeros_1 = np.clip(zeros_1,0,1)
        #zeros_frac[:,:,:,k] = zeros_1[:,:,:,k]/4

    
    ds_count = xr.Dataset(data_vars=dict(density=(['date','lat','lon'],zeros)),
                  coords=dict(lon=(['lon'],lons_uq),
                              lat=(['lat'],lats_uq),
                              time=(['date'],dates)))
    
    ds_count.to_netcdf(path+filez_to_write.format(int(range_years[0]),int(range_years[-1])))
    t1 = time.time()
    print('{} minutes'.format(t1-t0))

# %%
