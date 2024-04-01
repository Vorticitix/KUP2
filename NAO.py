#%%
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from eofs.standard import Eof
from eofs.examples import example_data_path
from glob import glob

#%%
path = "/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/atm/hist/"
files_1 = sorted(glob(path+'TRANS.3501BP.cam.h0.347?-0[12].nc'))
files_2 = sorted(glob(path+'TRANS.3501BP.cam.h0.347?-12.nc'))
files = sorted(files_1+files_2)
# %%
ds = xr.open_mfdataset(files)
slp = ds.PSL
# %%
