#%%
# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""

import xarray as xr
import numpy as np
from scipy.signal import argrelextrema
import pandas as pd
#%%
path = '/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/forcing_user'
file='/Evolk_EVA_distribution_1502BP-2015_extended.nc'
ds = xr.open_dataset(path+file)
colmass = ds.colmass.max('lat')
colmass[colmass<0.00002]=np.nan

idx = argrelextrema(colmass.values,np.greater)
quan = colmass[idx]
loc = ds.MMRVOLC[idx[0],:,:].sum('lev').idxmax('lat').to_pandas()
time = colmass.time[idx]

df = quan.sortby(quan, ascending=False).to_pandas()
df = pd.concat([df,loc],axis=1)
df.columns = ['colmass','lat']
df.to_csv('/storage/climatestor/PleioCEP/doensen/data/cyclone_volcano/eruptions.txt')
#ds['time'] =ds['time']-1502
max_per_month = ds.max('lat')
max_sorted = max_per_month[['time','colmass']].sortby(max_per_month.colmass,ascending=False)

# %%
