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
        'Z500':path+'Z500_V300/detrend/Z500_RCP85_anom_detrend.nc',
          'PREC':path+'PREC/detrend/PREC_RCP85_anom_detrend.nc',
            'WS':path+'WS/detrend/WS_RCP85_anom_detrend.nc',
            'T850':path+'T850/detrend/T850_RCP85_3507_3601_anom_detrend.nc',
            'RWP':path+'RWP/detrend/RWP_RCP85_all_anom_detrend.nc',
           'PSL':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_PSL_detrend.nc',
           'GZ':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_GZ_detrend.nc',
           'CZ':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_CZ_detrend.nc',
           'ZRAD':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_ZRAD_detrend.nc',
           'ZDEP':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_ZDEP_detrend.nc',
           'PRECMEAN':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_PRECMEAN_detrend.nc',
           'PRECMAX':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_PRECMAX_detrend.nc',
           'WSMEAN':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_WSMEAN_detrend.nc',
           'WSMAX':path_cyclone+'detrend/cyclone_stats_RCP85_all_monmean_XL_99p_WSMAX_detrend.nc',
          }

dic_1D = {
        'Z500':path+'Z500_V300/detrend/Z500_RCP85_anom_med_detrend.nc',
        'PREC':path+'PREC/detrend/PREC_RCP85_anom_med_detrend.nc',
          'WS':path+'WS/detrend/WS_RCP85_anom_med_detrend.nc',
          'T850':path+'T850/detrend/T850_RCP85_3507_3601_anom_med_detrend.nc',
          'RWP':path+'RWP/detrend/RWP_RCP85_all_anom_med_detrend.nc',
          'NAOI':path+'try_PSL/NAOI_RCP85.nc',
        'PSL':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_PSL_detrend.nc',
          'GZ':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_GZ_detrend.nc',
          'CZ':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_CZ_detrend.nc',
          'ZRAD':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_ZRAD_detrend.nc',
          'ZDEP':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_ZDEP_detrend.nc',
          'PRECMEAN':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_PRECMEAN_detrend.nc',
          'PRECMAX':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_PRECMAX_detrend.nc',
          'WSMEAN':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_WSMEAN_detrend.nc',
          'WSMAX':path_cyclone+'detrend/cyclone_stats_RCP85_all_med_XL_99p_WSMAX_detrend.nc'}

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
            #var_3D = var_3D.where(var_3D.time<=3352,drop=True)
            #var_1D = var_1D.where(var_1D.time<=3352,drop=True)
            #Make time dimension length equal
            if len(var_3D.time)>len(var_1D.time):
                var_1D = var_1D.reindex(time=var_3D.time)
            else:
                var_3D = var_3D.reindex(time=var_1D.time)
            
            #Define key to select data variable
            varkey_3D = list(var_3D.keys())[-1]; varkey_1D = list(var_1D.keys())[-1]
            
            mask = var_3D[varkey_3D].notnull().sum(dim='time')/len(var_3D.time)>.15
            var_3D_dat = var_3D[varkey_3D].rolling(time=90,min_periods=1).mean(dim='time',skipna=True)
            var_1D_dat = var_1D[varkey_1D].rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
                
            corr = xr.corr(var_3D_dat,var_1D_dat,dim='time').where(mask)
            fig,ax = plt.subplots(figsize=(16,8),subplot_kw={'projection': ccrs.PlateCarree()})
            fig.suptitle('RCP85 Correlation 1D: {} , 3D: {}, Season: {}'.format(key_1D,key_3D,key),fontsize=16)
            
            im=ax.contourf(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
                        levels=np.linspace(-.5,.5,21),cmap=cmap,extend='both')
            # ax.contour(lons_3D,lats_3D,corr,transform=ccrs.PlateCarree(),
            #             cmap='bwr',vmin=-.5,vmax=.5,levels=np.linspace(-1,1,21))
            fig.colorbar(im,orientation='horizontal')
            ax.coastlines()
            corr.to_netcdf(path_cyclone+'corr/detrend/corr_RCP85_1D_{}_3D_{}_{}_med_99p_XL_detrend.nc'.format(key_1D,key_3D,key))
            fig.tight_layout()
            fig.savefig(path_cyclone+'corr/detrend/corr_RCP85_1D_{}_3D_{}_{}_med_99p_XL_detrend.png'.format(key_1D,key_3D,key))
            plt.close(fig)
            






# %%
