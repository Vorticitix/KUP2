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
#%%
files = sorted(glob(path+'cyclone_count_cesm_????_????.nc'))
for file in files:
    ds = xr.open_dataset(file)
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    #ds = ds.sel(dayofyear=)
    ds = ds.where((ds.dayofyear>=335)|(ds.dayofyear<60)).dropna('dayofyear')
    ds_yy_mean = ds.mean(dim='dayofyear').mean(dim='year')
    #ds_yy_mean['density'] = ds_yy_mean['density']/90
    oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
    oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
        .sel(lat=slice(ds.lat.min(),ds.lat.max()),lon=slice(ds.lon.min(),ds.lon.max()))
    mountains = (oro.PHIS>9810)
    mountains = mountains.where(mountains==True)
    
    
    levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
            0.15,0.2,0.25,0.3,0.4,0.5]
    cmap = plt.get_cmap('hot_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
    fig,ax = plt.subplots(subplot_kw={'projection':ccrs.Orthographic(central_longitude=0,central_latitude=90)},
                          figsize=(8,12))
    im = ax.pcolormesh(ds_yy_mean.lon,ds_yy_mean.lat,ds_yy_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
    fig.colorbar(im,orientation='horizontal')
    ax.coastlines()
    ax.margins(0)
    #ax.set_extent([-120, 60, 20, 85])
    ax.set_global()
    gl = ax.gridlines(draw_labels=True)
    #gl.top_labels = False
    #gl.left_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
    fig.subplots_adjust(left=.1,right=.9,bottom=.05,top=.85)
    fig.suptitle('Cyclone Density (per day) {:04d}-{:04d}'.format(int(ds.year.min())
                                                          ,int(ds.year.max()))
                 ,fontsize=20)
    fig.savefig(path+'figs/count_{:04d}_{:04d}_cesm.png'.format(int(ds.year.min()),int(ds.year.max())))
    plt.close(fig)
    print(int(ds.year.min()),int(ds.year.max()))
    #fig.tight_layout()

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


# %%

#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

fig,axz = plt.subplots(4,3,figsize=(16,9),subplot_kw={'projection':ccrs.PlateCarree()})
for i,season in enumerate(seasons):
    doy1 = season[0]
    doy2 = season[1]
    if doy2<doy1:
        era_ss = era.where((era.dayofyear>=doy2)|(era.dayofyear<doy1)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy2)|(cesm.dayofyear<doy1)).dropna('dayofyear')
    else:
        era_ss = era.where((era.dayofyear>=doy1)&(era.dayofyear<doy2)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy1)&(cesm.dayofyear<doy2)).dropna('dayofyear') 
        
    era_ss_mean = era_ss.mean(dim='dayofyear').mean(dim='year')
    cesm_ss_mean = cesm_ss.mean(dim='dayofyear').mean(dim='year')
    diff_ss_mean = (cesm_ss_mean / era_ss_mean -1)*100
    ax_cesm = axz[i,0]
    ax_era = axz[i,1]
    ax_diff = axz[i,2]
    ax_cesm.text(-0.07, 0.55, season_labels[i], va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax_cesm.transAxes)
    if i==0:
        ax_cesm.set_title('CESM 1981-2010')
        ax_era.set_title('ERA5 1981-2010')
        ax_diff.set_title('Absolute Difference')
    
    im = ax_cesm.pcolormesh(cesm_ss_mean.lon,cesm_ss_mean.lat,
                            cesm_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_era.pcolormesh(era_ss_mean.lon,era_ss_mean.lat,
                            era_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_diff.pcolormesh(diff_ss_mean.lon,diff_ss_mean.lat,
                            diff_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='bwr',vmin=-100,vmax=100)
    for ax in axz.flat:
        ax.coastlines()
        ax.margins(0)
        ax.set_extent([-12.5, 40, 27, 45])

#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
#            dpi=300)

# %%
fig,axz = plt.subplots(3,4,figsize=(16,9),subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0, central_latitude=90)})
for i,season in enumerate(seasons):
    doy1 = season[0]
    doy2 = season[1]
    if doy2<doy1:
        era_ss = era.where((era.dayofyear>=doy2)|(era.dayofyear<doy1)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy2)|(cesm.dayofyear<doy1)).dropna('dayofyear')
    else:
        era_ss = era.where((era.dayofyear>=doy1)&(era.dayofyear<doy2)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy1)&(cesm.dayofyear<doy2)).dropna('dayofyear') 
        
    era_ss_mean = era_ss.mean(dim='dayofyear').mean(dim='year')
    cesm_ss_mean = cesm_ss.mean(dim='dayofyear').mean(dim='year')
    diff_ss_mean = cesm_ss_mean - era_ss_mean
    ax_cesm = axz[0,i]
    ax_era = axz[1,i]
    ax_diff = axz[2,i]
    ax_cesm.text(0.55, -0.07, season_labels[i], va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax_cesm.transAxes)
    if i==0:
        ax_cesm.set_title('CESM 1981-2010')
        ax_era.set_title('ERA5 1981-2010')
        ax_diff.set_title('Absolute Difference')
    
    im = ax_cesm.pcolormesh(cesm_ss_mean.lon,cesm_ss_mean.lat,
                            cesm_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_era.pcolormesh(era_ss_mean.lon,era_ss_mean.lat,
                            era_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_diff.pcolormesh(diff_ss_mean.lon,diff_ss_mean.lat,
                            diff_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='bwr')
    for ax in axz.flat:
        ax.coastlines()
        ax.margins(0)
        ax.set_global()
        # ax.set_extent([-12.5, 40, 27, 45])
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm_global.png',
            dpi=300)
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



# %%
#Comparison CESM 1981-2010 vs 2070-2099
#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

fig,axz = plt.subplots(4,3,figsize=(12,5.5),subplot_kw={'projection':ccrs.PlateCarree()})
for i,season in enumerate(seasons):
    doy1 = season[0]
    doy2 = season[1]
    if doy2<doy1:
        rcp_ss = rcp.where((rcp.dayofyear>=doy2)|(rcp.dayofyear<doy1)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy2)|(cesm.dayofyear<doy1)).dropna('dayofyear')
    else:
        rcp_ss = rcp.where((rcp.dayofyear>=doy1)&(rcp.dayofyear<doy2)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy1)&(cesm.dayofyear<doy2)).dropna('dayofyear') 
        
    rcp_ss_mean = rcp_ss.mean(dim='dayofyear').mean(dim='year')
    cesm_ss_mean = cesm_ss.mean(dim='dayofyear').mean(dim='year')
    diff_ss_mean =  rcp_ss_mean - cesm_ss_mean
    ax_cesm = axz[i,0]
    ax_rcp = axz[i,1]
    ax_diff = axz[i,2]
    ax_cesm.text(-0.07, 0.55, season_labels[i], va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax_cesm.transAxes)
    if i==0:
        ax_cesm.set_title('CESM 1981-2010')
        ax_rcp.set_title('RCP85 2070-2099')
        ax_diff.set_title('Absolute Difference')
    
    im = ax_cesm.pcolormesh(cesm_ss_mean.lon,cesm_ss_mean.lat,
                            cesm_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_rcp.pcolormesh(rcp_ss_mean.lon,rcp_ss_mean.lat,
                            rcp_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_diff.pcolormesh(diff_ss_mean.lon,diff_ss_mean.lat,
                            diff_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='bwr',vmin=-0.1,vmax=0.1)
    for ax in axz.flat:
        ax.coastlines()
        ax.margins(0)
        ax.set_extent([-12.5, 40, 27, 45])
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/rcp_vs_cesm.png',
            dpi=300)
# %%
fig,axz = plt.subplots(4,3,figsize=(10,12),subplot_kw={'projection':ccrs.Orthographic(central_longitude=0, central_latitude=90)})
for i,season in enumerate(seasons):
    doy1 = season[0]
    doy2 = season[1]
    if doy2<doy1:
        rcp_ss = rcp.where((rcp.dayofyear>=doy2)|(rcp.dayofyear<doy1)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy2)|(cesm.dayofyear<doy1)).dropna('dayofyear')
    else:
        rcp_ss = rcp.where((rcp.dayofyear>=doy1)&(rcp.dayofyear<doy2)).dropna('dayofyear')
        cesm_ss = cesm.where((cesm.dayofyear>=doy1)&(cesm.dayofyear<doy2)).dropna('dayofyear') 
        
    rcp_ss_mean = rcp_ss.mean(dim='dayofyear').mean(dim='year')
    cesm_ss_mean = cesm_ss.mean(dim='dayofyear').mean(dim='year')
    diff_ss_mean =  rcp_ss_mean - cesm_ss_mean
    ax_cesm = axz[i,0]
    ax_rcp = axz[i,1]
    ax_diff = axz[i,2]
    ax_cesm.text(-0.07, 0.55, season_labels[i], va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax_cesm.transAxes)
    if i==0:
        ax_cesm.set_title('CESM 1981-2010')
        ax_rcp.set_title('CESM 2070-2099')
        ax_diff.set_title('Absolute Difference')
    
    im = ax_cesm.pcolormesh(cesm_ss_mean.lon,cesm_ss_mean.lat,
                            cesm_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_rcp.pcolormesh(rcp_ss_mean.lon,rcp_ss_mean.lat,
                            rcp_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    ax_diff.pcolormesh(diff_ss_mean.lon,diff_ss_mean.lat,
                            diff_ss_mean.density,transform=ccrs.PlateCarree(),
                  cmap='bwr',vmin=-0.1,vmax=0.1)
    for ax in axz.flat:
        ax.coastlines()
        ax.margins(0)
        ax.set_global()
        # ax.set_extent([-12.5, 40, 27, 45])
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/rcp_vs_cesm_global.png',
            dpi=300)
# %%
