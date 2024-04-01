#%%
import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import wrf
from netCDF4 import Dataset
import cartopy.crs as ccrs
import sys
from matplotlib.animation import FuncAnimation
#%%
file_df = '/storage/climatestor/PleioCEP/doensen/data/fort_36_total_cmed.txt'
iic = 100244646
df = pd.read_csv(file_df)
# %%
idx = np.where(df['iic']==iic)
df_iic = df.iloc[idx]
# %%
file_cesm = '/storage/climatestor/PleioCEP/doensen/data/extracted/TRANS.3501BP.cam.h1.3350-01-01.sel_alt.nc'
file_wrf = '/storage/climatestor/PleioCEP/data/WRF_MedCyclones/Long/1848-01-02/wrfout_d01_1848-01-02_01:00:00'

ds = xr.open_dataset(file_cesm)
ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
ds = ds.sortby(ds['lon']).squeeze().sel(lat=slice(58,20),lon=slice(-25,50))
#%%
ncfile = Dataset(file_wrf)
slp = wrf.getvar(ncfile, "slp",wrf.ALL_TIMES)
p = wrf.getvar(ncfile, "pressure",wrf.ALL_TIMES)
z = wrf.getvar(ncfile, "z",wrf.ALL_TIMES)
ua =wrf. getvar(ncfile, "ua",wrf.ALL_TIMES)
va = wrf.getvar(ncfile, "va",wrf.ALL_TIMES)

u_850 = wrf.interplevel(ua, p, 850)
v_850 = wrf.interplevel(va, p, 850)
ws_850 = (u_850**2+v_850**2)**.5
#%%
lat_wrf,lon_wrf=wrf.latlon_coords(slp)
cart_proj = wrf.get_cartopy(slp)
savepath =  '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
#%%
def update(frame):
    lim = 196
    plt.clf()
    ax1 = plt.subplot(2,1,1,projection=ccrs.PlateCarree())
    lat_cesm = ds['lat']; lon_cesm = ds['lon']
    ax1.contour(lon_cesm,lat_cesm,ds['PSL'][((frame+1)//6+4),:,:]/100,levels=np.arange(950,1063,4),
                 colors='k',alpha=.5,transform=ccrs.PlateCarree())
    im = ax1.pcolormesh(lon_cesm,lat_cesm,ds['WS'][8,:,:],cmap='hot_r',
                   transform=ccrs.PlateCarree(),vmin=0,vmax=50)
    ax1.coastlines()
    ax1.set_title('CESM')
    ax2 = plt.subplot(2,1,2,projection=cart_proj)
    ax2.contour(lon_wrf[:lim,:],lat_wrf[:lim,:],slp[frame,:lim,:],levels=np.arange(950,1063,4),
                colors='k',alpha=.5,transform=ccrs.PlateCarree())
    ax2.pcolormesh(lon_wrf[:lim,:],lat_wrf[:lim,:],ws_850[frame,:lim,:],
                cmap='hot_r',transform=ccrs.PlateCarree(),vmin=0,vmax=50)
    ax2.coastlines()
    ax2.set_title('WRF')
    plt.subplots_adjust(left=.07,right=.95,bottom=.08,top=.92,wspace=.1,hspace=.1)
    cbar_ax = plt.axes([0.20, 0.045, 0.60, 0.015])
    cbar=plt.colorbar(im, cax=cbar_ax,orientation='horizontal',label='m/s')
    time = ds['time'][frame//6+4]
    plt.suptitle('850 hPa Wind Speed and Sea Level Pressure\n{}'.format(slp['Time'][frame].values.astype(str)[:19]))
                                                                           #int(time.dt.day),int(time.dt.month),int(time.dt.hour)))

fig = plt.figure(figsize=(12,10))
num_frames=96
animation = FuncAnimation(fig,update,frames=num_frames,repeat=False)
animation.save(savepath+'test.gif',writer='imagemagick', fps=6)
plt.show()