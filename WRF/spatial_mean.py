#%%
import xarray as xr
import numpy as np
from glob import glob
from matplotlib import pyplot as ptl
import sys
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/Long_1853_1887'
dirs = glob(path)
var = "T2"
# %%
for dir in dirs:
    means = []
    files=sorted(glob(dir+'/{}/*.nc'.format(var)))
    for file in files:
        print(file[78:-3])
        da = xr.open_dataset(file)[var]
        try:
            da = da.where((da['XLAT']<=47)&(da['XLAT']>=28)&
                        (da['XLONG']>=-10)&(da['XLONG']<=40))
        except:
            continue
        weights = np.cos(np.deg2rad(da.XLAT)) 
        a = da.weighted(weights).mean(('south_north','west_east'))
        means.append(a)
    mean = xr.concat(means,dim='Time')
    mean.to_netcdf(dir+'/{}/{}_mean_med.nc'.format(var,var))
# %%
