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
# %%
print('hallo')
#set path and filename to analyze 
path = "/storage/climatestor/PleioCEP/doensen/data/"
#List all necessary cyclone files
filez = sorted(glob(path+"cyclone_RCP85_3512_3601/fort_36_????_????"))
print(filez)
#output filename
filez_to_write = 'count_cyclone/cyclone_count_RCP85_{}_{}.nc'
#
start_year = 1
dy = 10
for j,file in enumerate(filez):
    header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
             'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmean',
		 'precmax','preccmean','preclmean','wsmean','wsmax']

    #Open 10 year cyclone script
    df = pd.read_csv(file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
    #print(df)    
    #Convert coordinate indexes to real  cartesian coordinates
    df['lat']=ilat_to_lat(df['ilat'])
    df['lon']=ilon_to_lon(df['ilon'])
    #df = df.where((df['lon']>=200)|(df['lon']<=60)).dropna().reset_index()
    #Exclude all dates that are lower than year 5
    df = df.where(df['date']>50000).dropna().reset_index()
    #Extract year, month and dayofyear from cyclone data
    df['month'] = [int(x[4:6]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['year'] = [int(x[:4]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['day'] = [int(x[6:8]) for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df = df.where(~((df['month']==2)&(df['day']==29))).dropna()
    #print(df)
    #df = df.where((df['year']>=start_year)&(df['year']<=start_year + dy -1)).dropna()
    df['dayofyear'] = [pd.to_datetime('2022-{}-{}'.format(x[4:6],x[6:8])).dayofyear\
                       for x in df['date'].astype(int).astype(str).str.zfill(8)]
    df['dayofyear'] = df['dayofyear'] + df['idumi']/240000
    #print(df)
    #Create arrays of all available cartesian coordinates and unique years
    lats_uq = df['lat'].sort_values(ascending=False).dropna().unique()
    lons_uq_sorted = df['lon'].sort_values().dropna().unique()
    lons_uq = np.append(lons_uq_sorted[lons_uq_sorted>60]
                        ,lons_uq_sorted[lons_uq_sorted<=60])
    lats_idx = np.arange(len(lats_uq))
    lons_idx = np.arange(len(lons_uq))
    range_years =np.unique(df['year'].values)
    #Create empty array with shape of 
    zeros_1 = np.zeros((len(lats_uq),len(lons_uq),len(range_years),1460))
    #   zeros_frac = np.zeros((len(lats_uq),len(lons_uq),len(range_years),1460))
    #zeros_2 = np.zeros((len(lats_uq),len(lons_uq),12))\
    count=0
    for j,year in enumerate(range_years):
        print(year)
        df_year = df.where(df['year']==year).dropna()
        range_dayofyears = np.unique(df_year['dayofyear'])
        for k,dayofyear in enumerate(range_dayofyears):
            print(k)
            df_dayofyear = df_year.where(df['dayofyear']==dayofyear).dropna()
            for i,latlon in df_dayofyear[['lat','lon','zrad']].iterrows():
                count+=1
                lat_fix,lon_fix,zrad = latlon.values
                zrad = zrad*2
                lat_idx = np.where(lat_fix==lats_uq)[0][0]
                lon_idx = np.where(lon_fix==lons_uq)[0][0]
                lon_idx_min = lon_idx - 14 ; lon_idx_max = lon_idx + 16
                lat_idx_min = lat_idx - 6 ; lat_idx_max = lat_idx + 7
                if lon_idx_min<0:
                    lon_idx_min=len(lons_uq) - lon_idx_min
                elif lon_idx_max>len(lons_uq):
                    lon_idx_max=lon_idx_max - len(lons_uq)
                if lat_idx_min<0:
                    lat_idx_min = 0
                elif lat_idx_max>len(lats_uq):
                    lat_idx_max=len(lats_uq)
                if lon_idx_min > lon_idx_max:
                    dist_1_left = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
                        for lon in lons_uq[lon_idx_min:len(lons_uq)]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
                    dist_1_right = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
                        for lon in lons_uq[0:lon_idx_max]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
                    dist_1_left_bool = (dist_1_left<=zrad*1000).astype(int)
                    dist_1_right_bool = (dist_1_right<=zrad*1000).astype(int)
                    zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:len(lons_uq),j,k] = np.add(dist_1_left_bool,zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:len(lons_uq),j,k])
                    zeros_1[lat_idx_min:lat_idx_max,0:lon_idx_max,j,k] = np.add(dist_1_right_bool,zeros_1[lat_idx_min:lat_idx_max,0:lon_idx_max,j,k])
                else:
                    dist_1 = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
                        for lon in lons_uq[lon_idx_min:lon_idx_max]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
                

                    dist_1_bool = (dist_1<=zrad*1000).astype(int)

                    zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,k] = np.add(dist_1_bool,zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,k])

            zeros_1 = np.clip(zeros_1,0,1)
            #zeros_frac[:,:,:,k] = zeros_1[:,:,:,k]/4

    ds_count = xr.Dataset(data_vars=dict(density=(['lat','lon','year','dayofyear'],zeros_1)),
                  coords=dict(lon=(['lon'],lons_uq),
                              lat=(['lat'],lats_uq),
                              dayofyear=(['dayofyear'],np.arange(1,1461)/4 + 1),
                              year=(['year'],range_years)))
    ds_count.to_netcdf(path+filez_to_write.format(int(range_years[0]),int(range_years[-1])))

# %%
