#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:05:07 2022

@author: doensen
"""

import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import sys
from scipy.ndimage.filters import gaussian_filter
import scipy
from matplotlib.colors import LinearSegmentedColormap

p = [-.6, -.2, .2, .6]
f = lambda x: np.interp(x, p, [0, 0.5, 0.5, 1])

cmap = LinearSegmentedColormap.from_list('map_white', 
              list(zip(np.linspace(0,1), plt.cm.bwr(f(np.linspace(min(p), max(p)))))))

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
#%%
seasons = {'MAM':[3,4,5],
              'JJA':[6,7,8],
              'SON':[9,10,11],
              'DJF':[12,1,2]}

dic_3D = {
    'PRECMEAN':path_cyclone+'cyclone_stats_monmean_PRECMEAN.nc',
        'Z500':path+'Z500_V300/Z500_anom.nc',
          'PREC':path+'PREC/PREC_anom.nc',
            'WS':path+'WS/WS_anom.nc',
            'T850':path+'T850/T850_0005_3514_anom.nc',
            'RWP':path+'RWP/RWP_all_anom.nc',
            'GZ':path_cyclone+'cyclone_stats_monmean_GZ.nc',
            'CZ':path_cyclone+'cyclone_stats_monmean_CZ.nc',
            'ZRAD':path_cyclone+'cyclone_stats_monmean_ZRAD.nc',
            'ZDEP':path_cyclone+'cyclone_stats_monmean_ZDEP.nc',
    
            #'PRECMAX':path_cyclone+'cyclone_stats_all_monmean_XL_99p_PRECMAX.nc',
            'WSMEAN':path_cyclone+'cyclone_stats_monmean_WSMEAN.nc',
            'PSL':path_cyclone+'cyclone_stats_monmean_PSL.nc',
            #'WSMAX':path_cyclone+'cyclone_stats_all_monmean_XL_99p_WSMAX.nc',
           'DENSITY':'/storage/climatestor/PleioCEP/doensen/data/count_cyclone/density_monmean.nc'
          }
#infile = path_cyclone + 'cyclone_spatial_stats_monmean.nc'
dic_1D = {
        'PC1':path+'EOF_season/PC_mode_01.nc',
        'PC2':path+'EOF_season/PC_mode_02.nc',
        'PC3':path+'EOF_season/PC_mode_03.nc',
        'PC4':path+'EOF_season/PC_mode_04.nc',
        # 'PC5':path+'EOF_season/PC_mode_05.nc',
        # 'PC6':path+'EOF_season/PC_mode_06.nc',
        # 'PC7':path+'EOF_season/PC_mode_07.nc',
        # 'PC8':path+'EOF_season/PC_mode_08.nc',
        'SST':path+'SST/SST_monmean_anom_med.nc'
        }

for key_3D in dic_3D:
    print(key_3D)
    #The spatial variable that is chosen to be correlated against a 1D timeseries
    var_3D_all = xr.open_dataset(dic_3D[key_3D]).squeeze()
    var_3D_all = var_3D_all.assign_coords(lon=(((var_3D_all.lon + 180) % 360) - 180)).sortby('lon')
    var_3D_all.attrs['history']=''
    var_3D_all = var_3D_all.sel(lat=slice(85,25),lon=slice(-100,60))
    lats_3D = var_3D_all.lat.values; lons_3D = var_3D_all.lon.values
    #var_3D_all = var_3D_all.where(var_3D_all.time.dt.year<=3352)

    #var_3D = var_3D.drop_vars('time_bnds')
    for key_1D in dic_1D:
    #The 1D variable that is compared against the spatial 3D field
        print(key_1D)
        var_1D_all = xr.open_dataset(dic_1D[key_1D]).squeeze()
        var_1D_all.attrs['history']=''

        #var_1D_all = var_1D_all.assign_coords(lon=(((vfrom matplotlib.colors import LinearSegmentedColormap
        #var_1D = var_1D_all.sel(lat=lats_3D,lon=lons_3D, method='nearest')

        for key in ['DJF','JJA']:
            #Select relevant season
            if key=='DJF':
                var_3D = var_3D_all.where((var_3D_all.time.dt.month>=seasons[key][0])|(var_3D_all.time.dt.month<=seasons[key][-1]),drop=True)
                var_1D = var_1D_all.where((var_1D_all.time.dt.month>=seasons[key][0])|(var_1D_all.time.dt.month<=seasons[key][-1]),drop=True)
            else:
                var_3D = var_3D_all.where((var_3D_all.time.dt.month>=seasons[key][0])&(var_3D_all.time.dt.month<=seasons[key][-1]),drop=True)
                var_1D = var_1D_all.where((var_1D_all.time.dt.month>=seasons[key][0])&(var_1D_all.time.dt.month<=seasons[key][-1]),drop=True)        
            
            #Convert from cftime to float as date
            var_3D_time = np.array([x.dt.year.values + x.dt.month.values/12 - 1/12 for x in var_3D.time])
            var_1D_time = np.array([x.dt.year.values + x.dt.month.values/12 - 1/12 for x in var_1D.time])
            #Remove duplictae values in time dimension and assign float to time dimension
            u, c = np.unique(var_3D_time, return_counts=True)
            sim  = u[c == 1]
            idx = np.where(np.isin(var_3D_time,sim)==True)[0]

            var_3D = var_3D.assign_coords(time=var_3D_time)
            var_3D = var_3D.isel(time=idx)
            
            u, c = np.unique(var_1D_time, return_counts=True)
            sim  = u[c == 1]
            idx = np.where(np.isin(var_1D_time,sim)==True)[0]
            
            var_1D = var_1D.assign_coords(time=var_1D_time)
            var_1D = var_1D.isel(time=idx)
            
            #Only take pre industrial values (pre 1850)
            var_3D = var_3D.where(var_3D.time<=3352,drop=True)
            var_1D = var_1D.where(var_1D.time<=3352,drop=True)
            #Make time dimension length equal
            if len(var_3D.time)>len(var_1D.time):
                var_1D = var_1D.reindex(time=var_3D.time)
            else:
                var_3D = var_3D.reindex(time=var_1D.time)
            
            #Define key to select data variable
            varkey_3D = list(var_3D.keys())[-1]; varkey_1D = list(var_1D.keys())[-1]
            
            mask = var_3D[varkey_3D].notnull().sum(dim='time')/len(var_3D.time)>.15
            var_3D_dat = var_3D[varkey_3D]#.rolling(time=90,min_periods=1).mean(dim='time',skipna=True)
            var_1D_dat = var_1D[varkey_1D]#.rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
                
            corr = xr.corr(var_3D_dat,var_1D_dat,dim='time').where(mask)
            fig,ax = plt.subplots(figsize=(16,8),subplot_kw={'projection': ccrs.PlateCarree()})
            fig.suptitle('Correlation 1D: {} , 3D: {}, Season: {}'.format(key_1D,key_3D,key),fontsize=16)
            
            im=ax.contourf(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                        levels=np.linspace(-.5,.5,21),cmap=cmap,extend='both')
            # ax.contour(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
            #             cmap='bwr',vmin=-.5,vmax=.5,levels=np.linspace(-1,1,21))
            fig.colorbar(im,orientation='horizontal')
            ax.coastlines()
            corr.to_netcdf(path_cyclone+'corr/corr_1D_{}_3D_{}_{}_whem_no_rolling.nc'.format(key_1D,key_3D,key))
            fig.tight_layout()
            fig.savefig(path_cyclone+'corr/corr_1D_{}_3D_{}_{}_whem_no_rolling.png'.format(key_1D,key_3D,key))
            plt.close(fig)

#%%

path_corr='/storage/climatestor/PleioCEP/doensen/data/cyclone_spatial_stats/corr/'
file = 'corr_1D_{}_3D_{}_{}_whem_alt.nc'
var_dic_gen = {
                'PC1':['PREC','WS','RWP','Z500'],
               'PC2':['T850','PREC','WS','Z500'],
               'PC3':['T850','PREC','WS','Z500'],
               'PC4':['T850','PREC','WS','RWP'],
               }

# Create a figure with 2x2 subplots
for region in ['med']:
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
            fig.subplots_adjust(left=.05,right=.95,bottom=0.1,top=.9)
            cax = fig.add_axes([0.05, 0.05, 0.9, 0.025])
            fig.colorbar(im, cax=cax,orientation='horizontal')

            # Show the plot
            fig.savefig(path_corr+'test_figs/{}_non_cyclone_{}_{}.png'.format(key,season,region))
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
file = 'corr_1D_{}_3D_{}_{}_whem_alt.nc'
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
            sys.exit()
            #fig.tight_layout()
            # Show the plot
            fig.savefig(path_corr+'test_figs/{}_cyclone_{}_{}_alt.png'.format(key,season,region))
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
