#%%
import xarray as xr
from glob import glob
import sys
import numpy as np
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/ERA5/'
JJA = False
#%%
# regions = ['emed','cmed','wmed','med']
# events = ['ws','prec','compound']
if JJA == True:
    regions = ['cmed_JJA']
    jjastring = '_JJA'
else:
    regions = ['emed','cmed']
    jjastring = ''
events = ['wsmax','precmax','compoundmax']
for region in regions:
    for event in events:
        tot_nc = []
        files = sorted(glob(path+'{}/{}_rank_*_all.nc'.format(region,event)))
        for i,file in enumerate(files):
            print(i)
            #rank = int(file[76:80])
            ds = xr.open_dataset(file)
            ds = ds.expand_dims(dim={"rank": 1})
            tot_nc.append(ds)
            if i+1==10:
                ds_tot = xr.concat(tot_nc,dim='rank')
                ds_mean = ds_tot.mean(dim='rank')
                ds_tot.to_netcdf(path+'{}_{}_top_{:04d}_all_ERA5.nc'.format(region,event,i+1))
                ds_mean.to_netcdf(path+'{}_{}_top_{:04d}_ERA5.nc'.format(region,event,i+1))
            if i+1==len(files):
                ds_tot = xr.concat(tot_nc,dim='rank')
                ds_mean = ds_tot.mean(dim='rank')
                ds_tot.to_netcdf(path+'{}_{}_top_{:04d}_all_ERA5.nc'.format(region,event,20))
                ds_mean.to_netcdf(path+'{}_{}_top_{:04d}_ERA5.nc'.format(region,event,20))
            # if i+1==len(files):
            #     ds_tot = xr.concat(tot_nc,dim='rank')
            #     ds_mean = ds_tot.mean(dim='rank')
            #     ds_tot.to_netcdf(path+'{}_{}_top_{:04d}_all.nc'.format(region,event,1000))
            #     ds_mean.to_netcdf(path+'{}_{}_top_{:04d}.nc'.format(region,event,1000))
            # if i+1>1000:
            #     break
        
# %%
