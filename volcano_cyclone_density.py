#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import sys
from scipy.stats import ttest_ind
import xarray as xr
from datetime import datetime
import cartopy.crs as ccrs
from matplotlib.colors import Normalize, BoundaryNorm
from scipy.stats import ttest_rel


path = '/storage/climatestor/PleioCEP/doensen/data/'
path2 = '/storage/climatestor/PleioCEP/doensen/data/extracted/'

def day_of_year_to_month(day_of_year):
    date = datetime.strptime(str(day_of_year), '%j')
    return date.month

#%%
#Load in list of strongest eruptions
erup = pd.read_csv(path+'cyclone_volcano/eruptions.txt')
#Sort eruption chronologically and test wether there is an eruption 10 years before an eruption and 5 years after an eruption
erup_sort = erup.sort_values('time')
bool_arr = np.full(erup.time.values.shape,True)
for i in np.arange(len(erup.time)):
    if i < len(erup_sort.time.values)-1:
        diff_lower = np.abs(erup_sort.time.values[i-1]-erup_sort.time.values[i])
        diff_higher = np.abs(erup_sort.time.values[i+1]-erup_sort.time.values[i])
        if (diff_lower<10)|(diff_higher<5):
            bool_arr[i]=False
    else:
        bool_arr[i]=False

#Only keep eruptions that do not overlap other eruptions
erup_sep_sort = erup_sort[bool_arr==True]
erup_sep = erup_sep_sort.sort_values('colmass',ascending=False)
#erup_sep = erup_sep[(erup_sep['lat']<20)&(erup_sep['lat']>-20)]
#Select 30 strongest eruptions
strong = erup_sep['time'][:20]
# %%
pre_tot = []
post_tot = []
post1_tot = []
post2_tot = []
post3_tot = []
post4_tot = []
post5_tot = []

for year in strong.values:
    year = int(year)
    lower =  int(round(year,-1)-15)
    higher = int(round(year,-1)+14)
    print(higher)
    print(lower)
    ds1 = xr.open_dataset(path+'count_cyclone/cyclone_count_cesm_{:04d}_{:04d}.nc'.format(lower,higher-20))
    ds2 = xr.open_dataset(path+'count_cyclone/cyclone_count_cesm_{:04d}_{:04d}.nc'.format(lower+10,higher-10))
    ds3 = xr.open_dataset(path+'count_cyclone/cyclone_count_cesm_{:04d}_{:04d}.nc'.format(lower+20,higher))
    ds = xr.concat([ds1,ds2,ds3],dim='date')
    # months = [day_of_year_to_month(int(x-.25)) for x in ds.dayofyear]
    # ds['dayofyear'] = months
    # ds = ds.rename({'dayofyear':'month'})
    ds = ds.where((ds.time.dt.month<=2)|(ds.time.dt.month>=12)).dropna(dim='date')
    sys.exit()
    ds['year'] = ds.time.dt.year
    ds = ds.mean(dim='date')
    pre_erup = ds.sel(year = slice(year-5,year-1))
    post_erup = ds.sel(year=slice(year+1,year+5))
    post_erup1 = ds.sel(year=year+1)
    post_erup2 = ds.sel(year=year+2)
    post_erup3 = ds.sel(year=year+3)
    post_erup4 = ds.sel(year=year+4)
    post_erup5 = ds.sel(year=year+5)
    pre_tot.append(pre_erup)
    post_tot.append(post_erup)
    post1_tot.append(post_erup1)
    post2_tot.append(post_erup2)
    post3_tot.append(post_erup3)
    post4_tot.append(post_erup4)
    post5_tot.append(post_erup5)

pre_xr = xr.concat(pre_tot,dim='year').mean(dim='year').sortby('lon')
post1_xr = xr.concat(post1_tot,dim='year').mean(dim='year')
post2_xr = xr.concat(post2_tot,dim='year').mean(dim='year')
post3_xr = xr.concat(post3_tot,dim='year').mean(dim='year')
post4_xr = xr.concat(post4_tot,dim='year').mean(dim='year')
post5_xr = xr.concat(post5_tot,dim='year').mean(dim='year')
post_xr = xr.concat(post_tot,dim='year').mean(dim='year').sortby('lon')
#%%
lat = pre_xr.lat; lon= pre_xr.lon.sortby('lon')
pre_xr_a = xr.concat(pre_tot,dim='year').transpose('year','lat','lon').sortby('lon')
post_xr_a = xr.concat(post_tot,dim='year').transpose('year','lat','lon').sortby('lon')

statres, pval = ttest_rel(pre_xr_a.density.values, post_xr_a.density.values)
latitude, longitude = np.meshgrid(lat, lon)

# pval = xr.Dataset(data_vars=dict(pval=(['lat','lon'],pval)), 
#                     coords=[('lat', lat.values), ('lon', lon.values)])
# pval = pval.sortby('lon')
fig,ax = plt.subplots(figsize=(16,9),subplot_kw={'projection': ccrs.PlateCarree()})
im=ax.pcolormesh(lon.values,lat.values,pval<0.05,transform=ccrs.PlateCarree())
ax.coastlines()
fig.colorbar(im,orientation='horizontal')

#%%

levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
fig,axz = plt.subplots(2,1,figsize=(12,9),subplot_kw={'projection': ccrs.PlateCarree()})
#total = [pre_xr,post1_xr,post2_xr,post3_xr,post4_xr,post5_xr]
total = [pre_xr,post_xr]
total = [x.sortby('lon') for x in total]
total = [x.where(x.density>0.01) for x in total]
for i,ax in enumerate(axz.flat):
    if i == 0:
        ax.set_extent([-90,50,20,80])
        im = ax.pcolormesh(lon.values,lat.values,total[i].density,transform=ccrs.PlateCarree(),
                      cmap='hot_r',norm=norm)
        ax.coastlines()
        fig.colorbar(im, ax=ax,orientation='vertical')
        
    else:
        ax.set_extent([-90,50,20,80])
        #total[i] = total[i].where(pval.values<0.05)
        im = ax.pcolormesh(lon.values,lat.values,total[i].density/total[0].density.values -1,transform=ccrs.PlateCarree(),
                      cmap='bwr',vmin=-.5,vmax=.5)
        
        ax.contourf(
            lon.values, lat.values, pval<0.05,
            transform=ccrs.PlateCarree(),
            colors='none',
            levels=[.5,1.5],
            hatches=['..']
        )
        ax.coastlines()
        fig.colorbar(im, ax=ax,orientation='vertical')

fig.tight_layout()
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_volcano/figs/cyclone_density.png',dpi=300)

# %%
lat = pre_xr.lat; lon= pre_xr.lon
levels=[0,0.01,0.02,0.03,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5]
cmap = plt.get_cmap('hot_r')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
fig,axz = plt.subplots(2,3,figsize=(16,7),subplot_kw={'projection': ccrs.PlateCarree()})
total = [pre_xr,post1_xr,post2_xr,post3_xr,post4_xr,post5_xr]
#total = [pre_xr,post_xr]
total = [x.where(x.density>0.01) for x in total]
for i,ax in enumerate(axz.flat):
    if i == 0:
        ax.set_extent([-90,50,20,80])
        im = ax.pcolormesh(lon.values,lat.values,total[i].density,transform=ccrs.PlateCarree(),
                      cmap='hot_r',norm=norm)
        ax.coastlines()
        fig.colorbar(im, ax=ax,orientation='horizontal')
        
    else:
        ax.set_extent([-90,50,20,80])
        im = ax.pcolormesh(lon.values,lat.values,total[i].density/total[0].density.values -1,transform=ccrs.PlateCarree(),
                      cmap='bwr',vmin=-.5,vmax=.5)
        ax.coastlines()
        fig.colorbar(im, ax=ax,orientation='horizontal')

fig.tight_layout()
#fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_volcano/figs/cyclone_density.png',dpi=300)
# %%
