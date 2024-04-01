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
df_iic = pd.read_csv(path+'fort_36_total_whem.txt')
df_iic = df_iic[['iic','lat','lon','date','idumi','slp']]
df_iic['lon'] = (df_iic['lon'] + 180) % 360 - 180

#%%
#Loop over regions
for region in locations:
    print(region)
    #Loop over eventtypes
    for event in events:

        check_rows=[]
        check_latlon = []
        print(event)
        #Open corresponding ranking for extreme event type and location
        df_ext = pd.read_csv(path+'cyclone_intensity/{}_ext_{}_CESM_modern.csv'.format(event,region))
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
            
            new_time_clim = xr.cftime_range('{:04d}-01-01'.format(cfdate.year),'{:04d}-01-01'.format(cfdate.year+1), freq='6H', calendar='noleap')[:-1]
            #cfdates = [cfdate + timedelta(hour=x) for x in np.linspace(-48,48,17)]
            file_nc = path+'extracted/TRANS.3501BP.cam.h1.{:04d}-01-01.sel_alt.nc'.format(cfdate.year)
            ds = xr.open_dataset(file_nc).squeeze()
            idx_ds = np.where(ds['time'].dt.year==cfdate.year)[0]
            ds = ds.isel(time=idx_ds)
            ds = ds.reindex(time=new_time_clim)
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby(ds.lon).squeeze()

            ds_tot = []
            
            #Open Cesm data 8 timesteps before and 8 timesteps after slp minimum
            idx_time = np.where(ds['time'].values==cfdate)[0]
            time_bnd = 8
            #If these timesteps are in multiple years open the corresponding file as well
            if idx_time<=time_bnd:
                print('-1')
                file_nc2 = path+'extracted/TRANS.3501BP.cam.h1.{:04d}-01-01.sel_alt.nc'.format(cfdate.year-1)
                if np.isin(cfdate.year-1,[1378,1583]):
                    continue
                new_time_clim = xr.cftime_range('{:04d}-01-01'.format(cfdate.year-1),'{:04d}-01-01'.format(cfdate.year), freq='6H', calendar='noleap')[:-1]
                # T850_clim_sub = T850_clim.copy(deep=True)
                # T850_clim_sub.coords['time'] = new_time_clim
                # T850_clim_sub = T850_clim_sub.reindex(time=T850_2.time)
                ds2 = xr.open_dataset(file_nc2).squeeze()
                idxs_ds = np.where(ds2['time'].dt.year==cfdate.year-1)[0]
                ds2 = ds2.isel(time=idxs_ds)
                ds2 = ds2.reindex(time=new_time_clim)
                ds2.coords['lon'] = (ds2.coords['lon'] + 180) % 360 - 180
                ds2 = ds2.sortby(ds2.lon).squeeze()

                ds = xr.concat([ds2,ds],dim='time')
                idx_time = np.where(ds['time'].values==cfdate)[0]
            elif idx_time>=(len(ds['time'])-time_bnd):
                print('+1')
                file_nc2 = path+'extracted/TRANS.3501BP.cam.h1.{:04d}-01-01.sel_alt.nc'.format(cfdate.year+1)
                if np.isin(cfdate.year+1,[1378,1583]):
                    continue

                new_time_clim = xr.cftime_range('{:04d}-01-01'.format(cfdate.year+1),'{:04d}-01-01'.format(cfdate.year+2), freq='6H', calendar='noleap')[:-1]

                ds2 = xr.open_dataset(file_nc2).squeeze()
                idxs_ds = np.where(ds2['time'].dt.year==cfdate.year+1)[0]
                ds2 = ds2.isel(time=idxs_ds)
                ds2 = ds2.reindex(time=new_time_clim)
                ds2.coords['lon'] = (ds2.coords['lon'] + 180) % 360 - 180
                ds2 = ds2.sortby(ds2.lon).squeeze()

                ds = xr.concat([ds,ds2],dim='time')


            ds_all = ds.copy(deep=True)
            #Set spatial data range    
            lon_bnd = 10
            lat_bnd = 10
            #Select all time steps from data

            ds_sub = ds_all.isel(time=np.arange(idx_time-time_bnd,idx_time+time_bnd+1))
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
                    idx_lat = np.argmin(np.abs(ds['lat'].values-lat_iic))
                    idx_lon = np.argmin(np.abs(ds['lon'].values-lon_iic))
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
                    ds_sub_lst.append(ds_sub_time)
            #Concatenate transofrmed time steps
            ds_sub_new = xr.concat(ds_sub_lst,dim='time') 
            #Save file to corresponding event and region directory
            ds_sub_new.to_netcdf(path+'cyclone_intensity/ERA5/{}/{}_rank_{:04d}_{}_all_CESM_modern.nc'.format(region,event,i,int(iic)))
        # check_latlon = pd.concat(check_latlon,axis=1).T.reset_index(drop=True)
        # check_latlon.to_csv(path+'cyclone_intensity/{}/{}_check_latlon.csv'.format(region,event))
        #check_df = pd.concat(check_rows,axis=1).T.reset_index(drop=True)
        #check_df.to_csv(path+'cyclone_intensity/{}/{}_check.csv'.format(region,event))
        
            #ds_sub_new.to_netcdf(path+'cyclone_intensity/emed/test.nc')
        


# %%
