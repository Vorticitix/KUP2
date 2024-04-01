#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:05:07 2022

@author: doensen
"""

import xarray  as xr
from matplotlib import pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import sys
from scipy.ndimage.filters import gaussian_filter
import scipy
from matplotlib.colors import LinearSegmentedColormap
import string

p = [-.5 ,-.2, .2, .5]
f = lambda x: np.interp(x, p, [0, 0.5, 0.5, 1])

cmap = LinearSegmentedColormap.from_list('map_white', 
              list(zip(np.linspace(0,1), plt.cm.seismic(f(np.linspace(min(p), max(p)))))))

path = '/storage/climatestor/PleioCEP/doensen/data/extracted/'
path_cyclone = '/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/'

def average_2D_matrix(old_array):
    array=old_array.copy()
    #get shape of 2D
    y,x = array.shape
    for i in range(1,x-1):
        for j in range(1,y-1):
            if np.isnan(array[j,i])==False:
                array[j,i]=np.nanmean(array[j-1:j+2,i-1:i+2])
            else:
                array[j,i]=np.nan
    array[0,:]=np.nan;array[y-1,:]=np.nan;array[:,0]=np.nan;array[:,x-1]=np.nan
    return array

oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/oro.nc')['PHIS']/9.81
oro.coords['lon']  =  (oro.coords['lon'] + 180) % 360 - 180
oro = oro.sortby(oro.lon).sel(lat=slice(25,85),lon=slice(-100,60))
oro = (oro>1000)
oro = oro.where(oro==True)

#%%
path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'
file = 'corr_1D_{}_3D_{}_{}_whem_no_rolling.nc'

var_dic_gen = {
                'PC1':['T850','PREC','WS','Z500'],
               'PC2':['T850','PREC','WS','Z500'],
               'PC3':['T850','PREC','WS','Z500'],
               'PC4':['T850','PREC','WS','Z500'],
               }

# Create a figure with 2x2 subplots

for season in ['DJF','JJA']:
    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(16, 12), subplot_kw={'projection': ccrs.PlateCarree()})
    #fig.suptitle('Correlation of {} against Mediterranean mean of other variables {}'.format(key,season))
    # Flatten the subplots array for easier iteration        
    # Loop over the subplots and plot the Plate Carree projection in each subplot


    for i, key in enumerate(var_dic_gen):
        for j, var in enumerate(var_dic_gen[key]):
            ax = axs[i,j]
            if j==0:
                ax.text(-0.03, 0.55, key, va='bottom', ha='center',
                rotation='vertical', rotation_mode='anchor',
                transform=ax.transAxes)
            # Set the extent of the map to show the whole world
            ds = xr.open_dataset(path_corr+file.format(key,var_dic_gen[key][j],season))
            var=list(ds.keys())[0]
            corr = ds[var]
            lons_3D = corr.lon.values; lats_3D = corr.lat.values
            #ax.set_global()
    
            # Draw coastlines and gridlines on the map
            ax.coastlines(zorder=1)
            # ax.gridlines(draw_labels=True)
            im=ax.pcolormesh(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                            cmap=cmap,vmin=-.75,vmax=.75,zorder=2)
            ax.pcolormesh(oro.lon,oro.lat,oro,transform=ccrs.PlateCarree(),
                           cmap='Greys',vmin=0,vmax=2,zorder=3)
            # Add a dummy title and axis labels
            ax.set_title(var_dic_gen[key][i])
            ax.label_outer()
                
        
        # Adjust the spacing between subplots
        #plt.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.subplots_adjust(left=.02,right=.98,bottom=0.1,top=.95,hspace=.02,wspace=.02)
    cax = fig.add_axes([0.02, 0.05, 0.96, 0.025])
    fig.colorbar(im, cax=cax,orientation='horizontal')
    sys.exit()
    # Show the plot
    fig.savefig(path_corr+'test_figs/all_PC_non_cyclone_{}.png'.format(season))
    plt.close()

#%%
var_dic_gen = {
               # 'SST' : ['PRECMEAN','DENSITY','GZ','WSMEAN']
                'PC1':['DENSITY','PRECMEAN','WSMEAN'],
                'PC2':['DENSITY','PRECMEAN','WSMEAN'],
                'PC3':['DENSITY','PRECMEAN','WSMEAN'],
                'PC4':['DENSITY','PRECMEAN','WSMEAN'],
               }

abet = list(string.ascii_lowercase)
ylabels = ['PC1 (NAO-like)','PC2 (EA-like)','PC3 (EAWR-like)','PC4 (SCAN-like)']
titles = ['Cyclone Frequency','Cyclone-related Precipitation','Cyclone-related wind speed']
path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'
file = 'corr_1D_{}_3D_{}_{}_whem_no_rolling.nc'
for season in ['DJF','JJA']:
    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(12, 8.5), subplot_kw={'projection': ccrs.PlateCarree()})
    #fig.suptitle('Correlation of {} against Mediterranean mean of other variables {}'.format(key,season))
    # Flatten the subplots array for easier iteration        
    # Loop over the subplots and plot the Plate Carree projection in each subplot

    c=0
    for i, key in enumerate(var_dic_gen):
        for j, var in enumerate(var_dic_gen[key]):
            ax = axs[i,j]
            ax.text(0.05, 1.15, abet[c]+')', transform=ax.transAxes,
              fontsize=14, fontweight='bold', va='top', ha='right')
            ax.set_extent([-100,60,25,85])
            c+=1
            if j==0:
                ax.text(-0.0225, 0.55, ylabels[i], va='bottom', ha='center',
                rotation='vertical', rotation_mode='anchor',
                transform=ax.transAxes)
            # Set the extent of the map to show the whole world
            ds = xr.open_dataset(path_corr+file.format(key,var_dic_gen[key][j],season))
            var=list(ds.keys())[0]
            corr = ds[var]
            lons_3D = corr.lon.values; lats_3D = corr.lat.values
            #ax.set_global()
    
            # Draw coastlines and gridlines on the map
            ax.coastlines(zorder=4)
            # ax.gridlines(draw_labels=True)
            im=ax.pcolormesh(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                            cmap=cmap,vmin=-.5,vmax=.5,zorder=2)
            ax.pcolormesh(oro.lon,oro.lat,oro,transform=ccrs.PlateCarree(),
                           cmap='Greys',vmin=0,vmax=2,zorder=3)
            # Add a dummy title and axis labels
            ax.set_title(titles[j])
            ax.label_outer()
                
        
        # Adjust the spacing between subplots
        #plt.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.subplots_adjust(left=.02,right=.98,bottom=0.08,top=.93,hspace=.02,wspace=.02)
    cax = fig.add_axes([0.02, 0.05, 0.96, 0.025])
    fig.colorbar(im, cax=cax,orientation='horizontal')
    fig.suptitle('Correlation Teleconnections and cyclone-related features {}'.format(season))
    # Show the plot
    fig.savefig(path_corr+'test_figs/all_PC_cyclone_no_rolling_{}.png'.format(season))
    plt.close()
#%%
var_dic_gen = {
               # 'SST' : ['PRECMEAN','DENSITY','GZ','WSMEAN']
                 'PC1':['PRECMEAN','DENSITY','GZ','WSMEAN'],
                'PC2':['PRECMEAN','DENSITY','GZ','WSMEAN'],
                'PC3':['PRECMEAN','DENSITY','GZ','WSMEAN'],
                'PC4':['PRECMEAN','DENSITY','GZ','WSMEAN'],
               }

path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'
file = 'corr_1D_{}_3D_{}_{}_whem_no_rolling.nc'
# Create a figure with 2x2 subplots
for region in ['whem']:
    for season in ['DJF','JJA']:
        for key in var_dic_gen:
            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree()})
            fig.suptitle('Correlation of {} against Mediterranean mean of other variables {}'.format(key,season))
            # Flatten the subplots array for easier iteration
            axs = axs.flatten()
            
            # Loop over the subplots and plot the Plate Carree projection in each subplot
            for i, ax in enumerate(axs):
        
                # Set the extent of the map to show the whole world
                ds = xr.open_dataset(path_corr+file.format(key,var_dic_gen[key][i],season,region))
                var=list(ds.keys())[0]
                corr = ds[var]
                lons_3D = corr.lon.values; lats_3D = corr.lat.values
                #ax.set_global()
        
                # Draw coastlines and gridlines on the map
                ax.coastlines()
                ax.gridlines()
                im=ax.pcolormesh(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                                cmap=cmap,vmin=-.5,vmax=.5)
                # Add a dummy title and axis labels
                ax.set_title(var_dic_gen[key][i])
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                
            
            # Adjust the spacing between subplots
            #plt.subplots_adjust(wspace=0.05, hspace=0.05)
            plt.subplots_adjust(left=.05,right=.95,bottom=0.1,top=.9,wspace=.05)
            cax = fig.add_axes([0.05, 0.05, 0.9, 0.025])
            fig.colorbar(im, cax=cax,orientation='horizontal')
            #fig.tight_layout()
            # Show the plot
            fig.savefig(path_corr+'test_figs/{}_cyclone_{}_{}_no_rolling.png'.format(key,season,region))
            plt.close()
#%%
var_dic_gen = {
                'PRECMEAN':['PREC','WS','RWP','Z500'],
               'WSMEAN':['PREC','WS','RWP','Z500'],
               'ZDEP':['PREC','WS','RWP','Z500'],
               'GZ':['PREC','WS','RWP','Z500'],
               }
file = 'corr_1D_{}_3D_{}_{}_{}_99p_XL.nc'
path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'

# Create a figure with 2x2 subplots
for region in ['med']:
    for season in ['DJF','JJA']:
        for key in var_dic_gen:
            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(16, 10), subplot_kw={'projection': ccrs.PlateCarree()})
            fig.suptitle('Correlation of {} against Mediterranean mean of other variables {}'.format(key,season))
            # Flatten the subplots array for easier iteration
            axs = axs.flatten()
            
            # Loop over the subplots and plot the Plate Carree projection in each subplot
            for i, ax in enumerate(axs):
        
                # Set the extent of the map to show the whole world
                ds = xr.open_dataset(path_corr+file.format(key,var_dic_gen[key][i],season,region))
                var=list(ds.keys())[0]
                corr = ds[var]
                lons_3D = corr.lon.values; lats_3D = corr.lat.values
                #ax.set_global()
        
                # Draw coastlines and gridlines on the map
                ax.coastlines()
                ax.gridlines()
                im=ax.pcolormesh(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                                cmap=cmap,vmin=-.5,vmax=.5)
                # Add a dummy title and axis labels
                ax.set_title(var_dic_gen[key][i])
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                
            
            # Adjust the spacing between subplots
            #plt.subplots_adjust(wspace=0.05, hspace=0.05)
            plt.subplots_adjust(left=.05,right=.95,bottom=0,top=.9,wspace=.05)
            fig.colorbar(im, ax=axs.ravel().tolist(),orientation='horizontal')
            #fig.tight_layout()
            # Show the plot
            fig.savefig(path_corr+'test_figs/{}_cyclone_{}_{}.png'.format(key,season,region))
            plt.close()
            
#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
cmap = colors.LinearSegmentedColormap.from_list('blue_white_red', [(0, 'blue'), (0.3, 'white'), (0.5, 'white'), (0.7, 'white'), (1, 'red')])

# create a discrete set of values ranging from -1 to 1
values = np.linspace(-.5, .5, num=11)

# create a colorbar with the discrete values mapped to colors
norm = colors.BoundaryNorm(values, cmap.N)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
file = 'corr_1D_{}_3D_{}_{}_{}_XL.nc'
path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(11, 4), 
                        subplot_kw={'projection': ccrs.PlateCarree()})
corr1=xr.open_dataset(path_corr+'corr_1D_WSMEAN_3D_Z500_DJF_med_XL.nc')
corr2=xr.open_dataset(path_corr+'corr_1D_WSMEAN_3D_RWP_DJF_med_XL.nc')
lons_3D = corr1.lon.values; lats_3D = corr1.lat.values
var=list(corr1.keys())[0]

ax1 = axs.flat[0]
im=ax1.pcolormesh(lons_3D,lats_3D,corr2[var],transform=ccrs.PlateCarree(),
    cmap=cmap,norm=norm)
ax1.set_extent([-50,40,25,75])
ax1.coastlines()
ax1.set_title(' with Rossby wave amplitude')
gl = ax1.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

ax2 = axs.flat[1]
im=ax2.pcolormesh(lons_3D,lats_3D,corr1[var],transform=ccrs.PlateCarree(),
    cmap=cmap,norm=norm)
ax2.set_extent([-50,40,25,75])
ax2.coastlines()
ax2.set_title('with 500 hPa geopotential height')
gl = ax2.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.025])

# add a horizontal colorbar to the new axis
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
fig.subplots_adjust(left=.1,right=.9,bottom=.1,top=.95)
cbar.ax.set_xlabel('correlation coefficient')
fig.suptitle()
# %%
