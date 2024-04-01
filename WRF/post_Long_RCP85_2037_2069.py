import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
import wrf 
from netCDF4 import Dataset
from glob import glob
import sys
from pathlib import Path
path = '/storage/climatestor/PleioCEP/data/WRF_MedCyclones/Long_RCP85_2037_2069/'
path_to_save = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/'
files=sorted(glob(path+'*/wrfout*'))

var = 'T2'
#%%
for file in files:
    wrfout = Dataset(file)
    Path(path_to_save+'Long_RCP85_2037_2069/{}/'.format(var)).mkdir(parents=True,exist_ok=True)
    da = wrf.getvar(wrfout,var,wrf.ALL_TIMES)
    da.attrs['projection'] = str(da.attrs['projection'])
    init_time = da.Time[0].dt.strftime("%Y-%m-%d")
    print(str(init_time.values))
    ds = da.to_dataset(name=var)
    ds.to_netcdf(path_to_save+'Long_RCP85_2037_2069/{}/{}_{}.nc'.format(var,var,str(init_time.values)))
