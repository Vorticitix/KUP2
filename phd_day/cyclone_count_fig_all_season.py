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
import matplotlib.gridspec as gridspec
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


# %%

#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,era]
titles = ['CESM','ERA5']
fig = plt.figure(figsize=(16,12))
grids = gridspec.GridSpec(2, 2, figure=fig)
k=0
for i in [0,1]:
    for j in [0,1]:
        
        subgrid1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=grids[i, j])
        ax1 = plt.subplot(subgrid1[0],projection=ccrs.PlateCarree())
        ax2 = plt.subplot(subgrid1[1],projection=ccrs.PlateCarree())
        if seasons[k][1]<seasons[k][0]:
            era_season = era.where((era.dayofyear>=seasons[k][0])|(era.dayofyear<seasons[k][1])).\
                dropna(dim='dayofyear').mean('dayofyear').mean('year')
            cesm_season = cesm.where((cesm.dayofyear>=seasons[k][0])|(cesm.dayofyear<seasons[k][1])).\
                dropna(dim='dayofyear').mean('dayofyear').mean('year')
        else:
            era_season = era.where((era.dayofyear>=seasons[k][0])&(era.dayofyear<seasons[k][1])).\
                dropna(dim='dayofyear').mean('dayofyear').mean('year')
            cesm_season = cesm.where((cesm.dayofyear>=seasons[k][0])&(cesm.dayofyear<seasons[k][1])).\
                dropna(dim='dayofyear').mean('dayofyear').mean('year')     
        oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
        oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
            .sel(lat=slice(era_season.lat.min(),era_season.lat.max()),
            lon=slice(era_season.lon.min(),era_season.lon.max()))   
        mountains = (oro.PHIS>9810)
        mountains = mountains.where(mountains==True)   
        levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
            0.15,0.2,0.25,0.3,0.4,0.5]
        cmap = plt.get_cmap('hot_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
        #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
        im = ax1.pcolormesh(cesm_season.lon,cesm_season.lat,cesm_season.density,transform=ccrs.PlateCarree(),
                        cmap='hot_r',norm=norm)
        ax1.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
        ax1.coastlines()
        ax1.margins(0)
        ax1.set_title(titles[0])
        ax1.set_global()
        ax1.set_extent([-70, 50, 25, 75])
        gl = ax1.gridlines(draw_labels=True)
        
        im = ax2.pcolormesh(era_season.lon,era_season.lat,era_season.density,transform=ccrs.PlateCarree(),
                        cmap='hot_r',norm=norm)
        ax2.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
        ax2.coastlines()
        ax2.margins(0)
        ax2.set_title(titles[1])
        ax2.set_global()
        ax2.set_extent([-70, 50, 25, 75])
        gl = ax2.gridlines(draw_labels=True)
    
    
        gl.top_labels = False
        gl.right_labels = False
        k+=1
#%%
        






        k+=1
#ds_yy_mean['density'] = ds_yy_mean['density']/90
    data = datas[i]
    data = data.where((data.dayofyear>=seasons[3][0])|(data.dayofyear<seasons[3][1])).\
        dropna(dim='dayofyear').mean('dayofyear').mean('year')
    oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
    oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
        .sel(lat=slice(data.lat.min(),data.lat.max()),
             lon=slice(data.lon.min(),data.lon.max()))
    mountains = (oro.PHIS>9810)
    mountains = mountains.where(mountains==True)
    levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
            0.15,0.2,0.25,0.3,0.4,0.5]
    cmap = plt.get_cmap('hot_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
    im = ax.pcolormesh(data.lon,data.lat,data.density,transform=ccrs.PlateCarree(),
                    cmap='hot_r',norm=norm)
    ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.margins(0)
    ax.set_title(titles[i])
    ax.set_global()
    ax.set_extent([-70, 50, 25, 75])
    gl = ax.gridlines(draw_labels=True)
    
    
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
fig.subplots_adjust(left=.07,right=.95,bottom=.11,top=.9,hspace=0.21)
cbar_ax = fig.add_axes([0.12, 0.05, 0.76, 0.02])
fig.colorbar(im, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
fig.suptitle('Cyclone Frequency (per day) 1981-2010 DJF',fontsize=14)
#fig.savefig(path+'figs/count_cesm_EGU.png',dpi=300)
#fig.tight_layout()

#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
#            dpi=300)


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
#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
datas = [cesm,rcp]
titles = ['CESM 1981-2010','CESM 2070-2099 (RCP8.5 Scenario)']
fig,axz = plt.subplots(2,1,figsize=(9,8),subplot_kw={'projection':ccrs.PlateCarree()})
for i,ax in enumerate(axz.flat):
    
#ds_yy_mean['density'] = ds_yy_mean['density']/90
    data = datas[i]
    data = data.where((data.dayofyear>=seasons[3][0])|(data.dayofyear<seasons[3][1])).\
        dropna(dim='dayofyear').mean('dayofyear').mean('year')
    if (datas[i].density.values==cesm.density.values).all():
        cesm_store = data.copy(deep=False)
    oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
    oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
        .sel(lat=slice(data.lat.min(),data.lat.max()),
             lon=slice(data.lon.min(),data.lon.max()))
    mountains = (oro.PHIS>9810)
    mountains = mountains.where(mountains==True)
    levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
            0.15,0.2,0.25,0.3,0.4,0.5]
    cmap = plt.get_cmap('hot_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
    if (datas[i].density.values==cesm.density.values).all():
        im = ax.pcolormesh(data.lon,data.lat,data.density,transform=ccrs.PlateCarree(),
                        cmap='hot_r',norm=norm,vmin=0.01,vmax=.5)
    elif (datas[i].density.values==rcp.density.values).all():
                im = ax.pcolormesh(data.lon,data.lat,data.density-cesm_store.density.values,transform=ccrs.PlateCarree(),
                        cmap='bwr',vmin=-0.1,vmax=0.1)
    ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
    ax.coastlines()
    #ax.margins(0)
    ax.set_title(titles[i])
    ax.set_global()
    ax.set_extent([-70, 50, 25, 75])
    gl = ax.gridlines(draw_labels=True)
    
    
    gl.top_labels = False
    gl.right_labels = False
    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
fig.subplots_adjust(left=.1,right=.9,bottom=.13,top=.915)
cbar_ax = fig.add_axes([0.12, 0.06, 0.76, 0.03])
fig.colorbar(im, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
fig.suptitle('Cyclone Density (per day) DJF',fontsize=14)
fig.savefig(path+'figs/count_cesm_rcp85.png',dpi=300)
#fig.tight_layout()

#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
#            dpi=300)



# %%
