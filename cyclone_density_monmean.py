#%%
import xarray as xr
import numpy as np
from glob import glob
from datetime import datetime
import sys
from cftime import DatetimeNoLeap
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/count_cyclone/'
name = 'cyclone_count_cesm_????_????.nc'
files = sorted(glob(path+name))
def day_of_year_to_month(day_of_year):
    date = datetime.strptime(str(day_of_year), '%j')
    return date.month
# %%
total = []
for file in files:
    print(file)
    ds = xr.open_dataset(file)
    months = [day_of_year_to_month(int(x-.25)) for x in ds.dayofyear]
    ds['dayofyear'] = months
    ds = ds.rename({'dayofyear':'month'})
    ds_mean = ds.groupby(ds.month).mean()
    years = ds_mean.year
    months = ds_mean.month
    date_list = []
    count = 0
    for year in years:
        for month in months:
            date = DatetimeNoLeap(year, month, 1)
            date_list.append(date)
    ds_flat = ds_mean.stack(time=['year','month'])
    ds_flat = ds_flat.assign_coords(time=date_list)
    ds_flat = ds_flat.transpose('time','lat','lon').sortby('lon')
    sys.exit()
    total.append(ds_flat)

ds_total = xr.concat(total,dim='time')
ds_total.to_netcdf(path+'density_monmean.nc')

# %%
