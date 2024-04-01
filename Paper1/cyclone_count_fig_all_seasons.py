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
#CESM vs rcp5 for the Mediterranean 1981-2010

files_era = sorted(glob(path+'cyclone_count_era_????_????_monmean.nc'))
files_cesm = [path+'cyclone_count_cesm_3475_3484_monmean.nc',
             path+'cyclone_count_cesm_3485_3494_monmean.nc',
             path+'cyclone_count_cesm_3495_3504_monmean.nc',
             path+'cyclone_count_cesm_3505_3514_monmean.nc']
files_rcp = [
             path+'cyclone_count_RCP85_3572_3581.nc',
             path+'cyclone_count_RCP85_3582_3591.nc',
             path+'cyclone_count_RCP85_3592_3601.nc']

era = xr.open_mfdataset(files_era)
cesm = xr.open_mfdataset(files_cesm)
cesm = cesm.where((cesm.time.dt.year>=3483)&(cesm.time.dt.year<=3512),drop=True)
rcp = xr.open_mfdataset(files_rcp)
era = era.assign_coords(lon=(((era.lon + 180) % 360) - 180)).sortby('lon')
cesm = cesm.assign_coords(lon=(((cesm.lon + 180) % 360) - 180)).sortby('lon')
rcp = rcp.assign_coords(lon=(((rcp.lon + 180) % 360) - 180)).sortby('lon')

# %%

diffs = ['absolute','relative']
#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
seasons_ = [(3,5),(6,8),(9,11),(12,2)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

titles = ['CESM','ERA5']
for diff_label in diffs:
    fig,axz = plt.subplots(4,3,figsize=(12,9),subplot_kw={'projection':ccrs.PlateCarree()})
    for i,label in enumerate(season_labels):
        ax_cesm = axz[i,0]
        ax_era = axz[i,1]
        ax_diff = axz[i,2]
    #ds_yy_mean['density'] = ds_yy_mean['density']/90
        if seasons[i][1]<seasons[i][0]:
            cesm_season = cesm.where((cesm.time.dt.month>=seasons_[i][0])|(cesm.time.dt.month<seasons_[i][1])).\
                mean('time')
            era_season = era.where((era.time.dt.month>=seasons_[i][0])|(era.time.dt.month<seasons_[i][1])).\
                mean('time')
        else:
            cesm_season = cesm.where((cesm.time.dt.month>=seasons_[i][0])&(cesm.time.dt.month<seasons_[i][1])).\
                mean('time')
            era_season = era.where((era.time.dt.month>=seasons_[i][0])&(era.time.dt.month<seasons_[i][1])).\
                mean('time')
        if diff_label == 'absolute':
            diff = cesm_season-era_season
        elif diff_label == 'relative':
            diff = cesm_season/era_season
        levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
                0.15,0.2,0.25,0.3,0.4,0.5]
        cmap = plt.get_cmap('hot_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
        #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
        im = ax_cesm.pcolormesh(cesm_season.lon,cesm_season.lat,cesm_season.density,transform=ccrs.PlateCarree(),
                        cmap='hot_r',norm=norm)
        ax_era.pcolormesh(era_season.lon,era_season.lat,era_season.density,transform=ccrs.PlateCarree(),
                        cmap='hot_r',norm=norm)
        if diff_label=='absolute':
            im2 = ax_diff.pcolormesh(diff.lon,diff.lat,diff.density,transform=ccrs.PlateCarree(),
                            cmap='seismic',vmin=-.2,vmax=.2)
        elif diff_label=='relative':
            diff = diff.where(cesm_season.density>0.005)
            diff = diff.where((diff.density>0.5)&(diff.density<2))
            im2 = ax_diff.pcolormesh(diff.lon,diff.lat,diff.density,transform=ccrs.PlateCarree(),
                            cmap='seismic',vmin=-0,vmax=2)
            
        gl = ax_cesm.gridlines(draw_labels=True)
        
        
        gl.top_labels = False
        gl.right_labels = False
        if i!=3:
            gl.bottom_labels = False
        gl = ax_era.gridlines(draw_labels=True)
        
        
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False
        if i!=3:
            gl.bottom_labels = False
        gl = ax_diff.gridlines(draw_labels=True)
        
        
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False
        if i!=3:
            gl.bottom_labels = False
    
        ax_cesm.set_title('CESM {}'.format(label))
        ax_era.set_title('ERA5 {}'.format(label))
        if diff_label == 'absolute':
            ax_diff.set_title('Absolute difference {}'.format(label))
        elif diff_label == 'relative':
            ax_diff.set_title('Relative difference {}'.format(label))
            
        
            
            
    
    oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
    oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
        .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
             lon=slice(cesm.lon.min(),cesm.lon.max()))
    mountains = (oro.PHIS>9810)
    mountains = mountains.where(mountains==True)
    for ax in axz.flat:
        ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
        ax.coastlines()
        ax.margins(0)
        ax.set_global()
        ax.set_extent([-70, 50, 25, 75])
    
    
        #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
    #fig.colorbar(im,orientation='horizontal')
    fig.subplots_adjust(left=.07,right=.95,bottom=.1,top=.93,wspace=.05,hspace=.05)
    cbar_ax = fig.add_axes([0.07, 0.065, 0.56, 0.015])
    cbar=fig.colorbar(im, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
    cbar.ax.tick_params(labelsize=12)
    cbar_ax = fig.add_axes([0.67, 0.065, 0.29, 0.015])
    cbar=fig.colorbar(im2, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
    formatter = mticker.FuncFormatter(lambda value, pos: f'{value:.2f}')
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.xaxis.set_major_formatter(formatter)
    fig.suptitle('Cyclone Frequency ($\mathregular{day^{-1}}$) 1981-2010 DJF',fontsize=16)
    if diff_label=='absolute':
        fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/count_cesm_all_seasons_absolute.png',dpi=400)
    elif diff_label=='relative':
        fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/count_cesm_all_seasons_relative.png',dpi=400)
    #fig.tight_layout()
    
    #fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/era_vs_cesm.png',
    #            dpi=300)



# %%
#Define first and last dayofyear for each season
seasons = [(60,152),(152,244),(244,335),(335,60)]
season_labels = ['MAM','JJA','SON','DJF']

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)

titles = ['CESM (1981-2010)','CESm RCP8.5 (2070-2099)']
fig,axz = plt.subplots(4,3,figsize=(12,9),subplot_kw={'projection':ccrs.PlateCarree()})
for i,label in enumerate(season_labels):
    ax_cesm = axz[i,0]
    ax_rcp = axz[i,1]
    ax_diff = axz[i,2]
#ds_yy_mean['density'] = ds_yy_mean['density']/90
    if seasons[i][1]<seasons[i][0]:
        cesm_season = cesm.where((cesm.time.dt.month>=seasons[i][0])|(cesm.time.dt.month<seasons[i][1])).\
            dropna(dim='dayofyear').mean('dayofyear').mean('year')
        rcp_season = rcp.where((rcp.dayofyear>=seasons[i][0])|(rcp.dayofyear<seasons[i][1])).\
            dropna(dim='dayofyear').mean('dayofyear').mean('year')
    else:
        cesm_season = cesm.where((cesm.time.dt.month>=seasons[i][0])&(cesm.time.dt.month<seasons[i][1])).\
            dropna(dim='dayofyear').mean('dayofyear').mean('year')
        rcp_season = rcp.where((rcp.dayofyear>=seasons[i][0])&(rcp.dayofyear<seasons[i][1])).\
            dropna(dim='dayofyear').mean('dayofyear').mean('year')
    diff = cesm_season-rcp_season
    levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,
            0.15,0.2,0.25,0.3,0.4,0.5]
    cmap = plt.get_cmap('hot_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    #nl_cmap = nlcmap(plt.get_cmap('hot_r'),levels)
    im = ax_cesm.pcolormesh(cesm_season.lon,cesm_season.lat,cesm_season.density,transform=ccrs.PlateCarree(),
                    cmap='hot_r',norm=norm)
    ax_rcp.pcolormesh(rcp_season.lon,rcp_season.lat,rcp_season.density,transform=ccrs.PlateCarree(),
                    cmap='hot_r',norm=norm)
    im2 = ax_diff.pcolormesh(diff.lon,diff.lat,diff.density,transform=ccrs.PlateCarree(),
                    cmap='bwr',vmin=-.2,vmax=.2)
    gl = ax_cesm.gridlines(draw_labels=True)
    
    
    gl.top_labels = False
    gl.right_labels = False
    if i!=3:
        gl.bottom_labels = False
    gl = ax_rcp.gridlines(draw_labels=True)
    
    
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    if i!=3:
        gl.bottom_labels = False
    gl = ax_diff.gridlines(draw_labels=True)
    
    
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    if i!=3:
        gl.bottom_labels = False

    ax_cesm.set_title('{}'.format(label))
    ax_rcp.set_title('{}'.format(label))
    ax_diff.set_title('Difference {}'.format(label))
    
        
        

oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')\
    .sel(lat=slice(cesm.lat.min(),cesm.lat.max()),
         lon=slice(cesm.lon.min(),cesm.lon.max()))
mountains = (oro.PHIS>9810)
mountains = mountains.where(mountains==True)
for ax in axz.flat:
    ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=1.5,transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.margins(0)
    ax.set_global()
    ax.set_extent([-70, 50, 25, 75])

    #gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
#fig.colorbar(im,orientation='horizontal')
fig.subplots_adjust(left=.07,right=.95,bottom=.1,top=.93,wspace=.05,hspace=.05)
cbar_ax = fig.add_axes([0.06, 0.065, 0.55, 0.015])
cbar=fig.colorbar(im, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
cbar.ax.tick_params(labelsize=12)
cbar_ax = fig.add_axes([0.67, 0.065, 0.3, 0.015])
cbar=fig.colorbar(im2, cax=cbar_ax,orientation='horizontal',label='$day^{-1}$')
cbar.ax.tick_params(labelsize=12)
fig.suptitle('Cyclone Frequency (per day) 1981-2010 DJF',fontsize=16)
#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/count_cyclone/figs/count_cesm_all_seasons_rcp.png',dpi=300)
#fig.tight_layout()



# %%
