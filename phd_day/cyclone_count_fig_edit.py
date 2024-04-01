#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:50:32 2022

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


# %%
#CESM vs ERA5 for the Mediterranean 1981-2010

files_era = sorted(glob(path+'cyclone_count_era_????_????.nc'))
files_cesm = [path+'cyclone_count_cesm_3475_3484.nc',
             path+'cyclone_count_cesm_3485_3494.nc',
             path+'cyclone_count_cesm_3495_3504.nc',
             path+'cyclone_count_cesm_3505_3514.nc']

era = xr.open_mfdataset(files_era)
cesm = xr.open_mfdataset(files_cesm).sel(year=slice(3483,3512))
era = era.assign_coords(lon=(((era.lon + 180) % 360) - 180)).sortby('lon')
cesm = cesm.assign_coords(lon=(((cesm.lon + 180) % 360) - 180)).sortby('lon')




#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,era]
titles = ['CESM',
          'Relative Difference of CESM compared to ERA5']
fig,axz = plt.subplots(2,1,figsize=(8,9),subplot_kw={'projection':ccrs.PlateCarree()})

ax1 = axz[0]; ax2 = axz[1]
#ds_yy_mean['density'] = ds_yy_mean['density']/90
cesm = cesm.where((cesm.dayofyear>=seasons[3][0])|(cesm.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
    .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
            lon=slice(cesm.lon.min(),cesm.lon.max()))
mountains = (oro.PHIS>9810)
mountains = mountains.where(mountains==True)
levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
        0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
#nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
im1 = ax1.pcolormesh(cesm.lon,cesm.lat,cesm.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax1.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax1.coastlines()
ax1.margins(0)
ax1.set_title(titles[0])
ax1.set_global()
ax1.set_extent([-70, 50, 25, 75])
gl = ax1.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False

era = era.where((era.dayofyear>=seasons[3][0])|(era.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
era = era.where(era.density>0.01)
im2 = ax2.pcolormesh(era.lon,era.lat,(cesm.density/era.density)*100,transform=ccrs.PlateCarree(),
                cmap='bwr',vmin=0,vmax=200)
ax2.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax2.coastlines()
ax2.margins(0)
ax2.set_title(titles[1])
ax2.set_global()
ax2.set_extent([-70, 50, 25, 75])
gl = ax2.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
#fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.2)
cbar_ax1 = fig.add_axes([0.12, 0.52, 0.76, 0.03])
fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal',label='$day^-1$')

fig.subplots_adjust(left=.1,right=.9,bottom=.10,top=.915,hspace=.43)
cbar_ax2 = fig.add_axes([0.12, 0.03, 0.76, 0.03])
fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal')
cbar_ax2.set_xticklabels(['-100%','-75%','-50%','-25%','0%',
                         '25%','50%','75%','100%'])
fig.suptitle('Cyclone Density 1981-2010 DJF',fontsize=14)
fig.savefig(path+'figs/count_cesm_era5_1981_2010_relative.png',dpi=300)
#fig.tight_layout()

#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
#            dpi=300)

#%%
files_era = sorted(glob(path+'cyclone_count_era_????_????.nc'))
files_cesm = [path+'cyclone_count_cesm_3475_3484.nc',
             path+'cyclone_count_cesm_3485_3494.nc',
             path+'cyclone_count_cesm_3495_3504.nc',
             path+'cyclone_count_cesm_3505_3514.nc']

era = xr.open_mfdataset(files_era)
cesm = xr.open_mfdataset(files_cesm).sel(year=slice(3483,3512))
era = era.assign_coords(lon=(((era.lon + 180) % 360) - 180)).sortby('lon')
cesm = cesm.assign_coords(lon=(((cesm.lon + 180) % 360) - 180)).sortby('lon')




#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,era]
titles = ['CESM','ERA5']
fig,axz = plt.subplots(2,1,figsize=(8,9),subplot_kw={'projection':ccrs.PlateCarree()})

ax1 = axz[0]; ax2 = axz[1]
#ds_yy_mean['density'] = ds_yy_mean['density']/90
cesm = cesm.where((cesm.dayofyear>=seasons[3][0])|(cesm.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
    .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
            lon=slice(cesm.lon.min(),cesm.lon.max()))
mountains = (oro.PHIS>9810)
mountains = mountains.where(mountains==True)
levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
        0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
#nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
im1 = ax1.pcolormesh(cesm.lon,cesm.lat,cesm.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax1.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax1.coastlines()
ax1.margins(0)
ax1.set_title(titles[0])
ax1.set_global()
ax1.set_extent([-70, 50, 25, 75])
gl = ax1.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False

era = era.where((era.dayofyear>=seasons[3][0])|(era.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
im2=ax2.pcolormesh(era.lon,era.lat,era.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax2.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())

ax2.coastlines()
ax2.margins(0)
ax2.set_title(titles[1])
ax2.set_global()
ax2.set_extent([-70, 50, 25, 75])
gl = ax2.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
#fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.2)
# cbar_ax1 = fig.add_axes([0.12, 0.52, 0.76, 0.03])
# fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal')

fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.2)
cbar_ax2 = fig.add_axes([0.12, 0.06, 0.76, 0.03])
fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal',label='$day^-1$')

fig.suptitle('Cyclone Density 1981-2010 DJF',fontsize=14)
fig.savefig(path+'figs/count_cesm_era5_1981_2010_bothmean.png',dpi=300)

# %%
#CESM vs ERA5 for the Mediterranean 1981-2010

files_rcp85 = [path+'cyclone_count_RCP85_3572_3581.nc',
             path+'cyclone_count_RCP85_3582_3591.nc',
             path+'cyclone_count_RCP85_3592_3601.nc',]
files_cesm = [path+'cyclone_count_cesm_3475_3484.nc',
             path+'cyclone_count_cesm_3485_3494.nc',
             path+'cyclone_count_cesm_3495_3504.nc',
             path+'cyclone_count_cesm_3505_3514.nc']

rcp = xr.open_mfdataset(files_rcp85)
cesm = xr.open_mfdataset(files_cesm).sel(year=slice(3483,3512))
rcp = rcp.assign_coords(lon=(((rcp.lon + 180) % 360) - 180)).sortby('lon')
cesm = cesm.assign_coords(lon=(((cesm.lon + 180) % 360) - 180)).sortby('lon')




#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,era]
titles = ['CESM 1981-2010',
          'Relative Difference of CESM (2070-2099) compared to CESM (1981-2010)']
fig,axz = plt.subplots(2,1,figsize=(8,9),subplot_kw={'projection':ccrs.PlateCarree()})

ax1 = axz[0]; ax2 = axz[1]
#ds_yy_mean['density'] = ds_yy_mean['density']/90
cesm = cesm.where((cesm.dayofyear>=seasons[3][0])|(cesm.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
    .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
            lon=slice(cesm.lon.min(),cesm.lon.max()))
mountains = (oro.PHIS>9810)
mountains = mountains.where(mountains==True)
levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
        0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
#nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
im1 = ax1.pcolormesh(cesm.lon,cesm.lat,cesm.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax1.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax1.coastlines()
ax1.margins(0)
ax1.set_title(titles[0])
ax1.set_global()
ax1.set_extent([-70, 50, 25, 75])
gl = ax1.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False

rcp = rcp.where((rcp.dayofyear>=seasons[3][0])|(rcp.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
rcp = rcp.where(rcp.density>0.01)
im2 = ax2.pcolormesh(rcp.lon,rcp.lat,(rcp.density/cesm.density)*100,transform=ccrs.PlateCarree(),
                cmap='bwr',vmin=0,vmax=200)
ax2.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax2.coastlines()
ax2.margins(0)
ax2.set_title(titles[1])
ax2.set_global()
ax2.set_extent([-70, 50, 25, 75])
gl = ax2.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
#fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.2)
cbar_ax1 = fig.add_axes([0.12, 0.52, 0.76, 0.03])
fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal',label='$day^-1$')

fig.subplots_adjust(left=.1,right=.9,bottom=.10,top=.915,hspace=.43)
cbar_ax2 = fig.add_axes([0.12, 0.03, 0.76, 0.03])
fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal')
cbar_ax2.set_xticklabels(['-100%','-75%','-50%','-25%','0%',
                         '25%','50%','75%','100%'])
fig.suptitle('Cyclone Density DJF',fontsize=14)
fig.savefig(path+'figs/count_cesm_rcp_2070_2099_relative.png',dpi=300)
#fig.tight_layout()

#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
#            dpi=300)



#%%
files_rcp85 = [path+'cyclone_count_RCP85_3572_3581.nc',
             path+'cyclone_count_RCP85_3582_3591.nc',
             path+'cyclone_count_RCP85_3592_3601.nc',]
files_cesm = [path+'cyclone_count_cesm_3475_3484.nc',
             path+'cyclone_count_cesm_3485_3494.nc',
             path+'cyclone_count_cesm_3495_3504.nc',
             path+'cyclone_count_cesm_3505_3514.nc']

rcp = xr.open_mfdataset(files_rcp85)
cesm = xr.open_mfdataset(files_cesm).sel(year=slice(3483,3512))
rcp = rcp.assign_coords(lon=(((rcp.lon + 180) % 360) - 180)).sortby('lon')
cesm = cesm.assign_coords(lon=(((cesm.lon + 180) % 360) - 180)).sortby('lon')



#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,era]
titles = ['CESM (1981-2010)',
          'CESM (2070-2099)']
fig,axz = plt.subplots(2,1,figsize=(8,9),subplot_kw={'projection':ccrs.PlateCarree()})

ax1 = axz[0]; ax2 = axz[1]
#ds_yy_mean['density'] = ds_yy_mean['density']/90
cesm = cesm.where((cesm.dayofyear>=seasons[3][0])|(cesm.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
    .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
            lon=slice(cesm.lon.min(),cesm.lon.max()))
mountains = (oro.PHIS>9810)
mountains = mountains.where(mountains==True)
levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
        0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
#nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
im1 = ax1.pcolormesh(cesm.lon,cesm.lat,cesm.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax1.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
ax1.coastlines()
ax1.margins(0)
ax1.set_title(titles[0])
ax1.set_global()
ax1.set_extent([-70, 50, 25, 75])
gl = ax1.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False

rcp = rcp.where((rcp.dayofyear>=seasons[3][0])|(rcp.dayofyear<seasons[3][1])).\
    dropna(dim='dayofyear').mean('dayofyear').mean('year')
im2=ax2.pcolormesh(rcp.lon,era.lat,rcp.density,transform=ccrs.PlateCarree(),
                cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
ax2.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())

ax2.coastlines()
ax2.margins(0)
ax2.set_title(titles[1])
ax2.set_global()
ax2.set_extent([-70, 50, 25, 75])
gl = ax2.gridlines(draw_labels=True)


gl.top_labels = False
gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.2)
# cbar_ax1 = fig.add_axes([0.12, 0.52, 0.76, 0.03])
# fig.colorbar(im1, cax=cbar_ax1,orientation='horizontal')

#fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915,hspace=.4)
cbar_ax2 = fig.add_axes([0.12, 0.06, 0.76, 0.03])
fig.colorbar(im2, cax=cbar_ax2,orientation='horizontal',label='$day^-1$')

fig.suptitle('Cyclone Density 1981-2010 DJF',fontsize=14)
fig.savefig(path+'figs/count_cesm_rcp_1981_2010_bothmean.png',dpi=300)

# %%
