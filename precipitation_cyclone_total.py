#%%
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cftime
file_total = '/storage/climatestor/PleioCEP/doensen/data/cyclone_era5_1981_2010/ERA5_prec_monsum.nc'
file_cyclone = '/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/ERA5_cyclone_spatial_stats_PRECMEAN_.nc'
#%%
ds_tot = xr.open_dataset(file_total)
ds_tot.coords['lon'] = (ds_tot.coords['lon'] + 180) % 360 - 180
ds_tot = ds_tot.sortby('lon')
ds_cyc = xr.open_dataset(file_cyclone)
ds_tot = ds_tot.sel(lat=ds_cyc.lat,lon=ds_cyc.lon,method='nearest')
newtime = [ cftime.DatetimeNoLeap(year=int(x.dt.year),month=int(x.dt.month),day=int(x.dt.day),\
            hour=int(x.dt.hour),has_year_zero=True) for x in ds_tot.time]
ds_tot['time']=newtime
#%%
days = ds_tot.time.dt.day==1
years = ds_tot.time.dt.year>=5
total = ds_tot.prec[days&years,:,:].drop_duplicates('time').sortby('time')
cyclone = ds_cyc.PRECMEAN.drop_duplicates('time').sortby('time')
common = np.intersect1d(total.time.values, cyclone.time.values)

cyclone = cyclone.sel(time=common)
total = total.sel(time=common)
#%%
freq_era = ds_cyc['PRECMEAN']/ds_tot['prec']
#freq.to_dataset(name='freq').to_netcdf('/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/cyclone_prec_explained_prec_all_era.nc')
#%%
file_total = '/storage/climatestor/PleioCEP/doensen/data/extracted/PREC/PREC_all.nc'
file_cyclone = '/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/cyclone_spatial_stats_monmean_PRECMEAN_.nc'
#%%
ds_tot = xr.open_dataset(file_total)
ds_tot.coords['lon'] = (ds_tot.coords['lon'] + 180) % 360 - 180
ds_tot = ds_tot.sortby('lon')
ds_cyc = xr.open_dataset(file_cyclone)
ds_tot = ds_tot.sel(lat=ds_cyc.lat,lon=ds_cyc.lon,method='nearest')
newtime = [ cftime.DatetimeNoLeap(year=int(x.dt.year),month=int(x.dt.month),day=int(x.dt.day),\
            hour=int(x.dt.hour),has_year_zero=True) for x in ds_tot.time]
ds_tot['time']=newtime
#%%
days = ds_tot.time.dt.day==1
years = ds_tot.time.dt.year>=5
total = ds_tot.PRECT[days&years,:,:].drop_duplicates('time').sortby('time')
cyclone = ds_cyc.PRECMEAN.drop_duplicates('time').sortby('time')
common = np.intersect1d(total.time.values, cyclone.time.values)
common_year = np.array([x.year for x in common])
common = common[(common_year>=3483)&(common_year<=3512)]
cyclone = cyclone.sel(time=common)
total = total.sel(time=common)
#%%
freq_cesm = cyclone/total

#freq.to_dataset(name='freq').to_netcdf('/storage/climatestor/PleioCEP/d
#%%
freq_era = freq_era.sel(lat=freq_cesm.lat,lon=freq_cesm.lon,method='nearest')
fig,axz=plt.subplots(2,1,figsize=(16,9),subplot_kw={'projection': ccrs.PlateCarree()})
freq_cesm_mean = freq_cesm.where((freq_cesm.time.dt.month>=12)|(freq_cesm.time.dt.month<=2)).mean(dim='time')
freq_era_mean = freq_era.where((freq_era.time.dt.month>=12)|(freq_era.time.dt.month<=2)).mean(dim='time')
ax1 = axz.flat[0]
im=ax1.pcolormesh(freq_cesm.lon,freq_cesm.lat,freq_cesm_mean*100,cmap='hot_r',
               vmin=0,vmax=100)
ax1.coastlines()
ax1.set_title('CESM')
ax2 = axz.flat[1]
ax2.pcolormesh(freq_era.lon,freq_era.lat,freq_era_mean*100,cmap='hot_r',
               vmin=0,vmax=100)
ax2.coastlines()
ax2.set_title('ERA5')
fig.subplots_adjust(left=.07,right=.95,bottom=.1,top=.97,wspace=.05,hspace=.1)
fig.suptitle('')
cbar_ax = fig.add_axes([0.235, 0.065, 0.55, 0.015])
cbar=fig.colorbar(im, cax=cbar_ax,orientation='horizontal',label='%')

fig.savefig('/storage/climatestor/PleioCEP/doensen/data/figs/cesm_vs_era_z.png',dpi=300)
