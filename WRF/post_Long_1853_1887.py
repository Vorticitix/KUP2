import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import wrf 
from netCDF4 import Dataset
from glob import glob
import sys
from pathlib import Path
path = '/storage/climatestor/PleioCEP/data/WRF_MedCyclones/Long_1853_1887/'
path_to_save = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/'
files=sorted(glob(path+'*/wrfout*'))

var = 'T2'
#%%
for file in files:
    if int(file[66:70])<1885:
        continue
.    if '.gz' in file:
        continue
    try:
        wrfout = Dataset(file)
    except:
        continue
    Path(path_to_save+'Long_1853_1887/{}/'.format(var)).mkdir(parents=True,exist_ok=True)
    da = wrf.getvar(wrfout,var,wrf.ALL_TIMES)
    da.attrs['projection'] = str(da.attrs['projection'])
    init_time = da.Time[0].dt.strftime("%Y-%m-%d")
    print(str(init_time.values))
    ds = da.to_dataset(name=var)
    ds.to_netcdf(path_to_save+'Long_1853_1887/{}/{}_{}.nc'.format(var,var,str(init_time.values)))
