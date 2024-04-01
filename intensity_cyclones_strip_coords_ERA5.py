#%%
import pandas as pd
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import cftime
import sys 
from datetime import timedelta
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/'
locations = {
               'emed':[(17.6,40),(28,42)],
               'cmed':[(2.6,17.5),(28,47)],
               #'wmed':[(-10,2.5),(28,40)],
               #'med':[(-10,40),(28,47)]
              }
events = ['precmax','wsmax','compoundmax']

#%%
df_iic = pd.read_csv(path+'fort_36_total_whem_ERA5.txt')
df_iic = df_iic[['iic','lat','lon','date','idumi','slp']]
df_iic['lon'] = (df_iic['lon'] + 180) % 360 - 180

#%%


new_time = xr.cftime_range('3352-01-01','3353-01-01', freq='1MS', calendar='noleap')

ERA5_all = xr.open_dataset(path+'cyclone_era5_1981_2010/ERA5_merge_all.nc')
#ERA5_sfc = xr.open_dataset(path+'cyclone_era5_1981_2010/ERA5_merge_prec_all.nc')
print('NetCDF files opened')
ERA5 = ERA5_all.sel(level=850)
ERA5 = ERA5.assign(WS=np.sqrt(ERA5['u']**2+ERA5['v']**2))

#%%
ERA5.coords['lon'] = (ERA5.coords['lon'] + 180) % 360 - 180
ERA5 = ERA5.sortby(ERA5.lon).squeeze()
#Loop over regions
for region in locations:
    print(region)
    #Loop over eventtypes
    for event in events:

        check_rows=[]
        check_latlon = []
        print(event)
        #Open corresponding ranking for extreme event type and location
        df_ext = pd.read_csv(path+'cyclone_intensity/{}_ext_{}_ERA5_JJA.csv'.format(event,region))
        #Transform latitudes
        df_ext['lon'] = (df_ext['lon'] + 180) % 360 - 180
        
        for i,iic in enumerate(df_ext['iic'].values[0:20]):
            print(i)
            #Select cyclones from ranking. From strongest to less strong
            idx = np.where(df_ext['iic']==iic)[0]
            #Locate cyclone in full dataset
            idx_iic = np.where(df_iic['iic']==iic)[0]
            #Select cyclones from both datasets
            df_sel = df_ext.iloc[idx]
            df_full_track = df_iic.iloc[idx_iic]
            #For calculating minimum only consider minima located in the region of interest
            #However, the full track outside the region is considered as well
            df_full_track_reg = df_full_track.where((df_full_track['lon']>=locations[region][0][0])&
                                                     (df_full_track['lon']<=locations[region][0][1])&
                                                     (df_full_track['lat']>=locations[region][1][0])&
                                                     (df_full_track['lat']<=locations[region][1][1]))
            if len(df_full_track_reg)==0:
                continue
            #Compute slp minimum and center datframe around this minimum.
            idx_slp_min = df_full_track_reg['slp'].reset_index().idxmin().values[1]
            min_bnd = -idx_slp_min
            max_bnd = len(df_full_track['slp'])-idx_slp_min
            df_full_track.index = np.arange(min_bnd,max_bnd)
            lat_series = df_full_track.reindex(np.linspace(-8,8,17))['lat'].round(1)
            lon_series = df_full_track.reindex(np.linspace(-8,8,17))['lon']
            latlon_series = pd.Series(list(zip(lat_series, lon_series)), 
                                      index=list(zip(lat_series.index,lat_series.index)))
            check_rows.append(df_full_track.reindex(np.linspace(-8,8,17))['slp'])
            check_latlon.append(latlon_series)

            # if df_full_track_reg.isnull().any().any():
            #     sys.exit()
            #Convert date of minimum, convert this to cftime object and use this to open CESM data
            date = str(df_sel['datetime'].astype(str).str.zfill(10).values[0])
            cfdate = cftime.DatetimeNoLeap(year=int(date[:4]),month=int(date[4:6]),day=int(date[6:8]),hour=int(date[8:10]))
            npdt = np.datetime64('{:04d}-{:02d}-{:02d}T{:02d}:00'.format(cfdate.year,cfdate.month,cfdate.day,cfdate.hour))
            
            #new_time_clim = np.arange(np.datetime64('{:04d}-01-01'.format(cfdate.year)),'{:04d}-01-01'.format(cfdate.year+1), freq='6H', calendar='noleap')[:-1]
            #cfdates = [cfdate + timedelta(hour=x) for x in np.linspace(-48,48,17)]
            #file_nc = path+'era5/ERA5_{:04d}_sfc.nc'.format(cfdate.year)
            # ds = ERA5
            # ds = xr.open_dataset(file_nc).squeeze()
            # idx_ds = np.where(ds['time'].dt.year==cfdate.year)[0]
            # ds = ds.isel(time=idx_ds)
            # #ds = ds.reindex(time=new_time_clim)
            # ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180
            # ds = ds.sortby(ds.longitude).squeeze()
            # ds = ds.assign(PRECT=ds['mcpr']+ds['mlspr'])
            
            # file_nc = path+'era5/ERA5_{:04d}.nc'.format(cfdate.year)
            # ds_ws = xr.open_dataset(file_nc).squeeze()
            # idx_ds_ws = np.where(ds_ws['time'].dt.year==cfdate.year)[0]
            # ds_ws = ds_ws.isel(time=idx_ds_ws)
            # #ds_ws = ds_ws.reindex(time=new_time_clim)
            # ds_ws.coords['longitude'] = (ds_ws.coords['longitude'] + 180) % 360 - 180
            # ds_ws = ds_ws.sortby(ds_ws.longitude).squeeze()
            # ws = np.sqrt(ds_ws.sel(level=850)['u']**2+ds_ws.sel(level=850)['v']**2).to_dataset(name='WS')
            
            # ds = xr.merge([ds,ws],compat='override')
            
            #Open Cesm data 8 timesteps before and 8 timesteps after slp minimum
            idx_time = np.where(ERA5['time'].values==npdt)[0]
            time_bnd = 8
            #If these timesteps are in multiple years open the corresponding file as well
            #ds_all = ds.copy(deep=True)
            #ds_all = ds_all.rename({'longitude':'lon','latitude':'lat'})
            #Set spatial data range    
            lon_bnd = 10
            lat_bnd = 10
            #Select all time steps from data

            ds_sub = ERA5.isel(time=np.arange(idx_time-time_bnd,idx_time+time_bnd+1))
            #Set actual times to indices around the slp minimum
            new_time = np.linspace(-time_bnd,time_bnd,2*time_bnd+1)
            ds_sub['time'] = new_time
            #Loop over time steps
            ds_sub_lst = []
            for time in ds_sub['time'].values:
                #Select time step
                ds_sub_time = ds_sub.sel(time=time)
                if np.isin(time,df_full_track.index):
                    #if time step present in cyclone track, center grid around the cyclone center
                    #Cyclone center will have the coordinate (0,0)
                    lat_iic = df_full_track.loc[time]['lat']
                    lon_iic = df_full_track.loc[time]['lon']
                    idx_lat = np.argmin(np.abs(ERA5['lat'].values-lat_iic))
                    idx_lon = np.argmin(np.abs(ERA5['lon'].values-lon_iic))
                    lat_weights = ds_sub_time.isel(lat=slice(idx_lat-lat_bnd,idx_lat+lat_bnd+1))['lat']
                    ds_sub_time = ds_sub_time.isel(lat=np.arange(idx_lat-lat_bnd,idx_lat+lat_bnd+1),
                                                   lon=np.arange(idx_lon-lon_bnd,idx_lon+lon_bnd+1)) 
                    new_lat = ds_sub_time['lat'].values - lat_iic
                    new_lon = ds_sub_time['lon'].values - lon_iic
                    #Account for numerical erros
                    ds_sub_time['lat'] = new_lat.round(7)
                    ds_sub_time['lon'] = new_lon.round(7)
                    ds_sub_lst.append(ds_sub_time)
                else:
                    #if time step not in cyclone track, set all values to NaN
                    #However, grid still needs to have a matching grid
                    gen_lat = 48
                    gen_lon = 72
                    ds_sub_time = ds_sub_time.isel(lat=np.arange(gen_lat-lat_bnd,gen_lat+lat_bnd+1),
                                                   lon=np.arange(gen_lon-lon_bnd,gen_lon+lon_bnd+1))
                    new_lat = ds_sub_time['lat'].values - ds_sub_time['lat'][lat_bnd].values
                    new_lon = ds_sub_time['lon'].values - ds_sub_time['lon'][lon_bnd].values
                    
                    ds_sub_time['lat'] = new_lat.round(7)
                    ds_sub_time['lon'] = new_lon.round(7)
                    
                    for var in list(ds_sub.keys()):
                        ds_sub_time[var][:]=np.nan
                    #print(ds_sub_time['lat'].values)
                    ds_sub_lst.append(ds_sub_time)
            #Concatenate transofrmed time steps
            ds_sub_new = xr.concat(ds_sub_lst,dim='time') 
            #Save file to corresponding event and region directory
            ds_sub_new.to_netcdf(path+'cyclone_intensity/ERA5/{}/{}_rank_{:04d}_{}_all_JJA.nc'.format(region,event,i,int(iic)))
        # check_latlon = pd.concat(check_latlon,axis=1).T.reset_index(drop=True)
        # check_latlon.to_csv(path+'cyclone_intensity/{}/{}_check_latlon.csv'.format(region,event))
        # check_df = pd.concat(check_rows,axis=1).T.reset_index(drop=True)
        # check_df.to_csv(path+'cyclone_intensity/{}/{}_check.csv'.format(region,event))
        
        
            #ds_sub_new.to_netcdf(path+'cyclone_intensity/emed/test.nc')
        


# %%
