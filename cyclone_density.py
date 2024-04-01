#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
from my_tools import *
from matplotlib import cm,colors
import cartopy.crs as ccrs
import xarray as xr
from geopy.distance import great_circle
import matplotlib.ticker as mticker

# %%

#set path and filename to analyze 
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone_3005_3514/"
filez = ["fort_36_3505_3514"]
#output filename
filez_to_write = ['cyclone_count_cesm_test.nc']
#
start_year = 1
dy = 10

for j,file in enumerate(filez):
    header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
             'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmean',
		 'precmax','preccmean','preclmean','wsmean','wsmax']

    
    df = pd.read_csv(path+file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
    print(df)    
    #df = df.where(df['date']>1e7).dropna()
    df['lat']=ilat_to_lat(df['ilat'])
    df['lon']=ilon_to_lon(df['ilon'])
    df = df.where((df['lon']>=240)|(df['lon']<=60)).dropna().reset_index()
    df = df.where(df['date']>50000).dropna().reset_index()
    df['month'] = [int(x[4:6]) for x in df['date'].astype(str).str.zfill(10)]
    df['year'] = [int(x[:4]) for x in df['date'].astype(str).str.zfill(10)]
    print(df)
    #df = df.where((df['year']>=start_year)&(df['year']<=start_year + dy -1)).dropna()
    df['dayofyear'] = [pd.to_datetime('2022-{}-{}'.format(x[4:6],x[6:8])).dayofyear\
                       for x in df['date'].astype(str).str.zfill(10)]
    print(df)
    lats_uq = df['lat'].sort_values(ascending=False).dropna().unique()
    lons_uq_sorted = df['lon'].sort_values().dropna().unique()
    lons_uq = np.append(lons_uq_sorted[lons_uq_sorted>60]
                        ,lons_uq_sorted[lons_uq_sorted<=60])
    lats_idx = np.arange(len(lats_uq))
    lons_idx = np.arange(len(lons_uq))
    range_years =np.unique(df['year'].values)
    zeros_1 = np.zeros((len(lats_uq),len(lons_uq),len(range_years),365))
    zeros_frac = np.zeros((len(lats_uq),len(lons_uq),len(range_years),365))
    #zeros_2 = np.zeros((len(lats_uq),len(lons_uq),12))\
    count=0
    for j,year in enumerate(range_years):
        print(year)
        df_year = df.where(df['year']==year).dropna()
        for dayofyear in np.arange(1,366):
            df_dayofyear = df_year.where(df['dayofyear']==dayofyear).dropna()
            for i,latlon in df_dayofyear[['lat','lon','zrad']].iterrows():
                count+=1
                lat_fix,lon_fix,zrad = latlon.values
                lat_idx = np.where(lat_fix==lats_uq)[0][0]
                lon_idx = np.where(lon_fix==lons_uq)[0][0]
                lon_idx_min = lon_idx - 14 ; lon_idx_max = lon_idx + 16
                lat_idx_min = lat_idx - 6 ; lat_idx_max = lat_idx + 7
                if lon_idx_min<0:
                    lon_idx_min=0
                elif lon_idx_max>len(lons_uq):
                    lon_idx_max=len(lons_uq)
                if lat_idx_min<0:
                    lat_idx_min = 0
                elif lat_idx_max>len(lats_uq):
                    lat_idx_max=len(lats_uq)
                
                dist_1 = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
                    for lon in lons_uq[lon_idx_min:lon_idx_max]] for lat in lats_uq[lat_idx_min:lat_idx_max]])
                
                #dist_2 = np.array([[great_circle((lat,lon),(lat_fix,lon_fix)).km\
                #         for lon in lons_uq] for lat in lats_uq])
                dist_1_bool = (dist_1<=zrad*1000).astype(int)
                #dist_2_bool = (dist_2<=zrad*1000).astype(int)
                #zeros_2[:,:,month-1] = np.add(dist_2_bool,zeros_2[:,:,month-1])
                zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,dayofyear-1] = np.add(dist_1_bool,zeros_1[lat_idx_min:lat_idx_max,lon_idx_min:lon_idx_max,j,dayofyear-1])
                # if (zeros_1[:,:,0] == zeros_2[:,:,0]).all()==False:
                #     print('Warning: Selected grid box not wide enough to cover entire cyclone')
                # if count%100==0:
                #     print(count)
        zeros_frac[:,:,:,dayofyear-1] = zeros_1[:,:,:,dayofyear-1]/4
    ds_count = xr.Dataset(data_vars=dict(density=(['lat','lon','year','dayofyear'],zeros_frac)),
                          coords=dict(lon=(['lon'],lons_uq),
                                      lat=(['lat'],lats_uq),
                                      dayofyear=(['dayofyear'],np.arange(1,366)),
                                      year=(['year'],range_years)))
    sys.exit()
    ds_count.to_netcdf(path+filez_to_write[j])
# %%

dic = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
nc_files = ['cyclone_count_cesm_PD.nc','RCP85/cyclone_count_cesm_RCP85.nc','cyclone_count_diff_frac_RCP85.nc']
ax_titles = ['Present Day (1980-2009)','RCP8.5 Scenario (2070-2099)','Relative Difference']
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')
oro_bool = (oro.PHIS<9810).sel(lat=slice(85,19),lon=slice(-120,60))
ticks_1 = np.array([0,0.01,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.4,0.5])
ticks_2 = np.linspace(0,200,9)
ticks_2_labels= ['0%','25%','50%','75%','100%','125%','150%','175%','200%']
cmap_1 = cm.get_cmap('hot_r')
cmap_2 = cm.get_cmap('bwr')
norm_1 = colors.BoundaryNorm(ticks_1, cmap_1.N,clip=True)
norm_2 = colors.BoundaryNorm(ticks_2, cmap_2.N,clip=True)
fig,axz = plt.subplots(4,3,figsize=(10,5.5),subplot_kw={'projection':ccrs.PlateCarree()})
for i,file in enumerate(nc_files):
    ds = xr.open_dataset(path+file)
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    ds['season'] = (['dayofyear'],pd.date_range('2022-01-01','2022-12-31'))
    ds.season.values = ds.season.dt.season
    #ds['season'] = pd.date_range('2022-01-01','2022-12-31')
    #ds = ds.assign_coords(season=ds['season'].dt.season.values)#.swap_dims({'dayofyear':'season'})    
    ds = ds.groupby('season').mean()
    

    for j,season in enumerate(['DJF','MAM','JJA','SON']):
        ax = axz[j,i]
        if j==0:
            ax.set_title(ax_titles[i],fontsize=14)
        if i==0:
            ax.text(-0.05, 0.5,season,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = ax.transAxes,
                 rotation=90, fontsize=12)
        ds_season = ds.sel(season=season)
        lon = ds_season.lon.values ; lat = ds_season.lat.values  
        #fraction = (ds.density /(dic[month]*4*30)).where(oro_bool.values)
        if i<2:
            im1 = ax.pcolormesh(lon,lat,ds_season.density.values,transform=ccrs.PlateCarree(),
                               cmap=cmap_1,norm=norm_1,)

        elif i==2:
            ds_season = ds_season.where(ds_bool.sel(season=season).density>0.005)
            ds_season.density.values = ds_season.density.values*100
            im2 = ax.pcolormesh(lon,lat,ds_season.density.values,transform=ccrs.PlateCarree(),
                                           cmap=cmap_2,norm=norm_2)

        ax.coastlines()
        ax.margins(0)
        ax.set_extent([-12.5, 40, 27, 45])
        mountains = (oro.PHIS>9810).where(oro.PHIS>9810,np.nan)
        ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=2)
    if i==0:
        cax1 = fig.add_axes([0.043, 0.2, 0.015, 0.6])
        fig.colorbar(im1,cax=cax1)
        cax1.yaxis.set_ticks_position('left')
        cax1.set_yticks(ticks_1)
        cax1.set_ylabel('Relative Amount of Cyclone Presence Time')
        cax1.yaxis.set_label_position("right")
        ds_bool = ds.copy()
    if i==2:
        cax2 = fig.add_axes([0.94, 0.2, 0.015, 0.6])
        cax2.set_yticks(ticks_2)
        #cax2.set_xticklabels(ticks_2_labels)
        cax2.yaxis.set_ticks_position('right')
        cax2.set_ylabel('Relative Difference PD vs RCP8.5')
        cax2.yaxis.set_label_position("right")
        cax2.text(-0.5, 0.5,'Relative Difference PD vs RCP8.5 [%]',
             horizontalalignment='center',
             verticalalignment='center',
             transform = cax2.transAxes,
             rotation=90, fontsize=10)
        fig.colorbar(im2,cax=cax2)
fig.subplots_adjust(left=0.10,bottom=0.015,top=0.9,right=0.915,
                    wspace=0.02,hspace=0.025)
fig.suptitle('Cyclone Density CESM',fontsize=16)
fig.savefig('/storage/mirrored/homes/doensen/poster_figs/cyclone_density.png',dpi=300)
plt.close(fig)

    

    #fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/"+save_files[i])
# %%
path = "/storage/climatestor/PleioCEP/doensen/data/"
nc_files = ['cyclone/cyclone_count_cesm_PD.nc','era5/cyclone_count_cesm_era5.nc']
ax_titles = ['CESM','ERA5']
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')
oro_bool = (oro.PHIS<9810).sel(lat=slice(85,19),lon=slice(-120,60))
ticks_1 = np.array([0,0.01,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.4,0.5])
ticks_2 = np.linspace(0,2,9)
cmap_1 = cm.get_cmap('hot_r')
cmap_2 = cm.get_cmap('bwr')
norm_1 = colors.BoundaryNorm(ticks_1, cmap_1.N,clip=True)
norm_2 = colors.BoundaryNorm(ticks_2, cmap_2.N,clip=True)
fig,axz = plt.subplots(2,1,figsize=(6.5,5),subplot_kw={'projection':ccrs.PlateCarree()})
for i,file in enumerate(nc_files):
    ax = axz[i]
    ax.margins(0)
    ds = xr.open_dataset(path+file)
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
    ds['season'] = (['dayofyear'],pd.date_range('2022-01-01','2022-12-31'))
    ds.season.values = ds.season.dt.season
    #ds['season'] = pd.date_range('2022-01-01','2022-12-31')
    #ds = ds.assign_coords(season=ds['season'].dt.season.values)#.swap_dims({'dayofyear':'season'})    
    ds_mean = ds.mean(dim='dayofyear')
    lon = ds.lon.values ; lat = ds.lat.values 
    if i<2:
        im1 = ax.pcolormesh(lon,lat,ds_mean.density.values,transform=ccrs.PlateCarree(),
                               cmap=cmap_1,norm=norm_1)
    # elif i ==2:
    #     ds_season = ds_season.where(ds_bool.sel(season=season).density>0.005)
    #     im2 = ax.pcolormesh(lon,lat,ds_mean.density.values,transform=ccrs.PlateCarree(),
    #                                  cmap=cmap_2,norm=norm_2)  
    if i==0:
        cax1 = fig.add_axes([0.075, 0.2, 0.025, 0.6])
        fig.colorbar(im1,cax=cax1)
        cax1.yaxis.set_ticks_position('left')
        cax1.set_yticks(ticks_1)
        cax1.set_ylabel('Relative Amount of Cyclone Presence Time')
        cax1.yaxis.set_label_position("right")
        ds_bool = ds.copy()
    # if i==2:
    #     cax2 = fig.add_axes([0.925, 0.2, 0.025, 0.6])
    #     cax2.set_yticks(ticks_2)
    #     cax2.yaxis.set_ticks_position('right')
    #     cax2.set_ylabel('Relative Difference PD vs RCP8.5')
    #     cax2.yaxis.set_label_position("right")
    #     cax2.text(-0.25, 0.5,'Relative Difference PD vs RCP8.5',
    #          horizontalalignment='center',
    #          verticalalignment='center',
    #          transform = cax2.transAxes,
    #          rotation=90, fontsize=10)
    #     fig.colorbar(im2,cax=cax2)
    fig.subplots_adjust(left=0.15,bottom=0.05,top=0.9,right=0.925,
                        wspace=0.02,hspace=0.07)
    ax.coastlines()
    ax.set_title(ax_titles[i])
    mountains = (oro.PHIS>9810).where(oro.PHIS>9810,np.nan)
    ax.pcolormesh(oro.lon,oro.lat,mountains,cmap='Greys',vmin=0,vmax=2)
    ax.set_extent([-120, 60, 20, 85])    
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.left_labels = False
    gl.xlocator = mticker.FixedLocator([-120,-90,-60,-30,0,30 ,60])
    if i ==0:
        gl.bottom_labels = False
    ax.label_outer()
fig.suptitle('          Cyclone Density Comparison CESM vs ERA5',fontsize=16)
fig.savefig('/storage/mirrored/homes/doensen/poster_figs/cyclone_density_era5.png',dpi=300)
    
    
    


# %%
dic = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
nc_files = ['cyclone_count_frac.nc']
save_files = ['cyclone_count_frac.png']
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro2.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')
oro_bool = (oro.PHIS<9810).sel(lat=slice(85,19),lon=slice(-120,60))
for i,file in enumerate(nc_files):
    fig,axz = plt.subplots(4,3,figsize=(16,9),subplot_kw={'projection':ccrs.PlateCarree()})
    for month in np.arange(1,13):
        ax = axz.flat[month-1]
        ds = xr.open_dataset(path+file).sel(month=month)
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
        lon = ds.lon.values ; lat = ds.lat.values  
        fraction = (ds.density /(dic[month]*4*30)).where(oro_bool.values)
        im = ax.pcolormesh(lon,lat,fraction.values,transform=ccrs.PlateCarree(),
                           vmin=-0.4,vmax=0.4,cmap='bwr')
        ax.coastlines()
        ax.margins(0)
    fig.subplots_adjust(left=0.02,bottom=0.015,top=0.965,right=0.935,
                        wspace=0.02,hspace=0.075)
    cax = fig.add_axes([0.95, 0.2, 0.025, 0.6])
    fig.colorbar(im,cax=cax)
    sys.exit()
    fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/"+save_files[i])
    
#%%


dic = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
nc_files = ['cyclone_count_cesm_PD.nc','RCP85/cyclone_count_cesm_RCP85.nc']
save_files = ['cyclone_count_cesm_PD.png','RCP85/cyclone_count_cesm_RCP85.png']
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')
oro_bool = (oro.PHIS<9810).sel(lat=slice(85,19),lon=slice(-120,60))
ticks= np.array([0,0.01,0.02,0.04,0.06,0.1,0.15,0.2,0.3,0.4,0.5])
cmap = cm.get_cmap('hot_r')
norm = colors.BoundaryNorm(ticks, cmap.N,clip=True)
for i,file in enumerate(nc_files):
    fig,axz = plt.subplots(4,3,figsize=(16,9),subplot_kw={'projection':ccrs.PlateCarree()})
    for month in np.arange(1,13):
        ax = axz.flat[month-1]
        ds = xr.open_dataset(path+file).sel(month=month)
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
        lon = ds.lon.values ; lat = ds.lat.values  
        fraction = (ds.density /(dic[month]*4*30)).where(oro_bool.values)
        im = ax.pcolormesh(lon,lat,fraction.values,transform=ccrs.PlateCarree(),
                           cmap=cmap,norm=norm)
        ax.coastlines()
        ax.margins(0)
        ax.set_title('Month = {}'.format(month))
        ax.set_extent([-15, 45, 25, 50])
    fig.subplots_adjust(left=0.02,bottom=0.015,top=0.965,right=0.935,
                        wspace=0.02,hspace=0.15)
    cax = fig.add_axes([0.94, 0.2, 0.025, 0.6])
    sys.exit()
    fig.colorbar(im,cax=cax)
    fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/"+save_files[i])

# %%

dic = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
nc_files = ['cyclone_count_diff.nc']
save_files = ['cyclone_count_diff.png']
oro = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro2.nc')
oro = oro.assign_coords(lon=(((oro.lon + 180) % 360) - 180)).sortby('lon')
oro_bool = (oro.PHIS<9810).sel(lat=slice(85,19),lon=slice(-120,60))
for i,file in enumerate(nc_files):
    fig,axz = plt.subplots(4,3,figsize=(16,9),subplot_kw={'projection':ccrs.PlateCarree()})
    for month in np.arange(1,13):
        ax = axz.flat[month-1]
        ds = xr.open_dataset(path+file).sel(month=month)
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')
        lon = ds.lon.values ; lat = ds.lat.values  
        fraction = (ds.density /(dic[month]*4*30)).where(oro_bool.values)
        im = ax.pcolormesh(lon,lat,fraction.values,transform=ccrs.PlateCarree(),
                           vmin=-0.3,vmax=0.3,cmap='bwr')
        ax.coastlines()
        ax.margins(0)
        ax.set_extent([-15, 45, 25, 50])
        ax.set_title('Month = {}'.format(month))
    fig.subplots_adjust(left=0.02,bottom=0.015,top=0.965,right=0.935,
                        wspace=0.02,hspace=0.115)
    cax = fig.add_axes([0.94, 0.2, 0.025, 0.6])
    fig.colorbar(im,cax=cax)
    sys.exit()
    fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/"+save_files[i])
