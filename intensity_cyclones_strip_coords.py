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
events = ['compound']

#%%
df_iic = pd.read_csv(path+'fort_36_total_whem.txt')
df_iic = df_iic[['iic','lat','lon','date','idumi','slp']]
df_iic['lon'] = (df_iic['lon'] + 180) % 360 - 180

#%%
Z500_V300 = xr.open_dataset(path+'extracted/Z500_V300/Z500_V300_all.nc')
Z500_V300.coords['lon'] = (Z500_V300.coords['lon'] + 180) % 360 - 180
Z500_V300 = Z500_V300.sortby(Z500_V300.lon).squeeze()




RWP = xr.open_dataset(path+'extracted/RWP/RWP_all.nc')
RWP.coords['lon'] = (RWP.coords['lon'] + 180) % 360 - 180
RWP = RWP.sortby(RWP.lon).squeeze()

new_time = xr.cftime_range('3352-01-01','3353-01-01', freq='1MS', calendar='noleap')
Z500_V300_clim = xr.open_dataset(path+'extracted/Z500_V300/Z500_V300_clim.nc')
Z500_V300_clim.coords['lon'] = (Z500_V300_clim.coords['lon'] + 180) % 360 - 180
Z500_V300_clim.coords['time'] = new_time[:-1]
Z500_V300_clim = Z500_V300_clim.sortby(Z500_V300_clim.lon).squeeze().reindex(time=new_time)
Z500_V300_clim = Z500_V300_clim.resample(time='6H').pad().isel(time=slice(0,1460))

T850_clim = xr.open_dataset(path+'extracted/T850/T850_0005_3352_ymonmean.nc')
T850_clim.coords['lon'] = (T850_clim.coords['lon'] + 180) % 360 - 180
T850_clim = T850_clim.sortby(T850_clim.lon).squeeze().reindex(time=new_time)
T850_clim = T850_clim.resample(time='6H').pad().isel(time=slice(0,1460))
print('NetCDF files opened')

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
        df_ext = pd.read_csv(path+'cyclone_intensity/{}_ext_{}_JJA.csv'.format(event,region))
        #Transform latitudes
        df_ext['lon'] = (df_ext['lon'] + 180) % 360 - 180
        
        for i,iic in enumerate(df_ext['iic'].values[0:1000]):
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

            try:
                T850 = xr.open_dataset(path+'extracted/T850/TRANS.3501BP.cam.h1.{:04d}_T850.nc'.format(cfdate.year))
            except:
                continue
            idxs_T850 = np.where(T850['time'].dt.year==cfdate.year)[0]
            T850 = T850.isel(time=idxs_T850)
            T850.coords['lon'] = (T850.coords['lon'] + 180) % 360 - 180
            T850 = T850.sortby(T850.lon).squeeze()
            T850 = T850.reindex(time=new_time_clim)
            # T850_clim_sub = T850_clim.copy(deep=True)
            # T850_clim_sub.coords['time'] = new_time_clim
            # T850_clim_sub = T850_clim_sub.reindex(time=T850.time)
            T850 = T850 - T850_clim.T.values

            try:
                U300 = xr.open_dataset(path+'extracted/U300/TRANS.3501BP.cam.h1.{:04d}_U300.nc'.format(cfdate.year))
            except:
                continue
            idxs_U300 = np.where(U300['time'].dt.year==cfdate.year)[0]
            U300 = U300.isel(time=idxs_U300)
            U300.coords['lon'] = (U300.coords['lon'] + 180) % 360 - 180
            U300 = U300.sortby(U300.lon).squeeze()
            U300 = U300.reindex(time=new_time_clim)

            idxs_Z = np.where(Z500_V300['time'].dt.year==cfdate.year)[0]
            if len(idxs_Z)==0:
                continue
            Z500_V300_year = Z500_V300.isel(time=idxs_Z)
            Z500_V300_year = Z500_V300_year.reindex(time=new_time_clim)
            Z500_year = Z500_V300_year['Z'].to_dataset()
            # Z500_V300_clim_sub = Z500_V300_clim.copy(deep=True)
            # Z500_V300_clim_sub.coords['time'] = new_time_clim
            # Z500_V300_clim_sub = Z500_V300_clim_sub.reindex(time=Z500_V300_year.time)
            Z500_year = Z500_year - Z500_V300_clim['Z'].values
            
            V300_year = Z500_V300_year['V'].to_dataset()

            idxs_RWP = np.where(RWP['time'].dt.year==cfdate.year)[0]
            if len(idxs_RWP)==0:
                continue
            RWP_year = RWP.isel(time=idxs_RWP)
            RWP_year = RWP_year.reindex(time=new_time_clim)


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
                try:
                    T850_2 = xr.open_dataset(path+'extracted/T850/TRANS.3501BP.cam.h1.{:04d}_T850.nc'.format(cfdate.year-1))
                except:
                    continue
                idxs_T850 = np.where(T850_2['time'].dt.year==cfdate.year-1)[0]
                T850_2 = T850_2.isel(time=idxs_T850)
                T850_2.coords['lon'] = (T850_2.coords['lon'] + 180) % 360 - 180
                T850_2 = T850_2.sortby(T850_2.lon).squeeze()
                T850_2 = T850_2.reindex(time=new_time_clim)
                # T850_clim_sub = T850_clim.copy(deep=True)
                # T850_clim_sub.coords['time'] = new_time_clim
                # T850_clim_sub = T850_clim_sub.reindex(time=T850_2.time)
                T850_2 = T850_2 - T850_clim.T.values
                
                try:
                    U300_2 = xr.open_dataset(path+'extracted/U300/TRANS.3501BP.cam.h1.{:04d}_U300.nc'.format(cfdate.year-1))
                except:
                    continue
                idxs_U300 = np.where(U300_2['time'].dt.year==cfdate.year-1)[0]
                U300_2 = U300_2.isel(time=idxs_U300)
                U300_2.coords['lon'] = (U300_2.coords['lon'] + 180) % 360 - 180 
                U300_2 = U300_2.reindex(time=new_time_clim)
                U300_2 = U300_2.sortby(U300_2.lon).squeeze()

                idxs_Z = np.where(Z500_V300['time'].dt.year==cfdate.year-1)[0]
                if len(idxs_Z)==0:
                    continue
                Z500_V300_year_2 = Z500_V300.isel(time=idxs_Z)
                Z500_V300_year_2 = Z500_V300_year_2.reindex(time=new_time_clim)
                Z500_year_2 = Z500_V300_year_2['Z'].to_dataset()
                # Z500_V300_clim_sub = Z500_V300_clim.copy(deep=True)
                # Z500_V300_clim_sub.coords['time'] = new_time_clim
                # Z500_V300_clim_sub = Z500_V300_clim_sub.reindex(time=Z500_V300_year_2.time)
                Z500_year_2 = Z500_year_2 - Z500_V300_clim['Z'].values
                V300_year_2 = Z500_V300_year_2['V'].to_dataset()

                idxs_RWP = np.where(RWP['time'].dt.year==cfdate.year-1)[0]
                RWP_year_2 = RWP.isel(time=idxs_RWP)
                RWP_year_2 = RWP_year_2.reindex(time=new_time_clim)

                ds2 = xr.open_dataset(file_nc2).squeeze()
                idxs_ds = np.where(ds2['time'].dt.year==cfdate.year-1)[0]
                ds2 = ds2.isel(time=idxs_ds)
                ds2 = ds2.reindex(time=new_time_clim)
                ds2.coords['lon'] = (ds2.coords['lon'] + 180) % 360 - 180
                ds2 = ds2.sortby(ds2.lon).squeeze()

                ds = xr.concat([ds2,ds],dim='time')
                T850 = xr.concat([T850_2,T850],dim='time')
                U300 = xr.concat([U300_2,U300],dim='time')
                Z500_year = xr.concat([Z500_year_2,Z500_year],dim='time')
                V300_year = xr.concat([V300_year_2,V300_year],dim='time')
                RWP_year = xr.concat([RWP_year_2,RWP_year],dim='time')

                idx_time = np.where(ds['time'].values==cfdate)[0]
            elif idx_time>=(len(ds['time'])-time_bnd):
                print('+1')
                file_nc2 = path+'extracted/TRANS.3501BP.cam.h1.{:04d}-01-01.sel_alt.nc'.format(cfdate.year+1)
                if np.isin(cfdate.year+1,[1378,1583]):
                    continue

                new_time_clim = xr.cftime_range('{:04d}-01-01'.format(cfdate.year+1),'{:04d}-01-01'.format(cfdate.year+2), freq='6H', calendar='noleap')[:-1]

                try:
                    T850_2 = xr.open_dataset(path+'extracted/T850/TRANS.3501BP.cam.h1.{:04d}_T850.nc'.format(cfdate.year-1))
                except:
                    continue   
                idxs_T850 = np.where(T850_2['time'].dt.year==cfdate.year+1)[0]
                T850_2 = T850_2.isel(time=idxs_T850) 
                T850_2.coords['lon'] = (T850_2.coords['lon'] + 180) % 360 - 180
                T850_2 = T850_2.sortby(T850_2.lon).squeeze()
                T850_2 = T850_2.reindex(time=new_time_clim)
                # T850_clim_sub = T850_clim.copy(deep=True)
                # T850_clim_sub.coords['time'] = new_time_clim
                # T850_clim_sub = T850_clim_sub.reindex(time=T850_2.time)
                T850_2 = T850_2 - T850_clim.T.values
                
                
                try:
                    U300_2 = xr.open_dataset(path+'extracted/U300/TRANS.3501BP.cam.h1.{:04d}_U300.nc'.format(cfdate.year-1))
                except:
                    continue
                idxs_U300 = np.where(U300_2['time'].dt.year==cfdate.year+1)[0]
                U300_2 = U300_2.isel(time=idxs_U300)
                U300_2.coords['lon'] = (U300_2.coords['lon'] + 180) % 360 - 180 
                U300_2 = U300_2.sortby(U300_2.lon).squeeze()
                U300_2 = U300_2.reindex(time=new_time_clim)
                
                idxs_Z = np.where(Z500_V300['time'].dt.year==cfdate.year+1)[0]
                if len(idxs_Z)==0:
                    continue
                Z500_V300_year_2 = Z500_V300.isel(time=idxs_Z)
                Z500_V300_year_2 = Z500_V300_year_2.reindex(time=new_time_clim)
                Z500_year_2 = Z500_V300_year_2['Z'].to_dataset()
                # Z500_V300_clim_sub = Z500_V300_clim.copy(deep=True)
                # Z500_V300_clim_sub.coords['time'] = new_time_clim
                # Z500_V300_clim_sub = Z500_V300_clim_sub.reindex(time=Z500_V300_year_2.time)
                Z500_year_2 = Z500_year_2 - Z500_V300_clim['Z'].values
                V300_year_2 = Z500_V300_year_2['V'].to_dataset()
                
                idxs_RWP = np.where(RWP['time'].dt.year==cfdate.year+1)[0]
                if len(idxs_RWP)==0:
                    continue
                RWP_year_2 = RWP.isel(time=idxs_RWP)
                RWP_year_2 = RWP_year_2.reindex(time=new_time_clim)

                ds2 = xr.open_dataset(file_nc2).squeeze()
                idxs_ds = np.where(ds2['time'].dt.year==cfdate.year+1)[0]
                ds2 = ds2.isel(time=idxs_ds)
                ds2 = ds2.reindex(time=new_time_clim)
                ds2.coords['lon'] = (ds2.coords['lon'] + 180) % 360 - 180
                ds2 = ds2.sortby(ds2.lon).squeeze()

                ds = xr.concat([ds,ds2],dim='time')
                T850 = xr.concat([T850,T850_2],dim='time')
                U300 = xr.concat([U300,U300_2],dim='time')
                Z500_year = xr.concat([Z500_year,Z500_year_2],dim='time')
                V300_year = xr.concat([V300_year,V300_year_2],dim='time')
                RWP_year = xr.concat([RWP_year,RWP_year_2],dim='time')

            WS300 = np.sqrt(U300['U']**2+V300_year['V']**2).to_dataset(name='WS300')
            RWP_year = RWP_year.rename_vars({'V':'RWP'})


            ds_all = xr.merge([ds,Z500_year,WS300,T850,RWP_year],compat='override')
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
            ds_sub_new.to_netcdf(path+'cyclone_intensity/{}/{}_rank_{:04d}_{}_all.nc'.format(region,event,i,int(iic)))
        # check_latlon = pd.concat(check_latlon,axis=1).T.reset_index(drop=True)
        # check_latlon.to_csv(path+'cyclone_intensity/{}/{}_check_latlon.csv'.format(region,event))
        check_df = pd.concat(check_rows,axis=1).T.reset_index(drop=True)
        check_df.to_csv(path+'cyclone_intensity/{}/{}_check.csv'.format(region,event))
        
            #ds_sub_new.to_netcdf(path+'cyclone_intensity/emed/test.nc')
        


# %%
