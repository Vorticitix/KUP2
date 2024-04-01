#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 09:15:49 2022

@author: doensen
"""
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

#%%
path = '/storage/climatestor/PleioCEP/doensen/data/extracted/'

file_past = 'total_mean_past_scenario.nc'
file_future = 'total_mean_future_scenario.nc'


past = xr.open_dataset(path+file_past).squeeze()
future = xr.open_dataset(path+file_future).squeeze()

past  = past.assign_coords(time=pd.date_range('1850-01-01','2012-12-01',freq='1MS'))
future = future.assign_coords(time=pd.date_range('2005-01-01','2099-12-01',freq='1MS'))

past_T = past.TBOT.rolling(time=12).mean()
past_co2 = past.co2vmr

future_T = future.TBOT.rolling(time=12).mean()
future_co2 = future.co2vmr

#%%
#time = pd.date_range('1850-01-01','2099-01-01',freq='MS')

fig,ax = plt.subplots()
ax.plot(past_T.time,past_T,label='Past')
ax.plot(future_T.time,future_T,label='RCP8.5')
ax.grid()
ax.set_ylabel('Global Mean T [K]')
ax.set_xlabel('year')
ax.legend()
fig.savefig('/storage/mirrored/homes/doensen/tmp_figs/comparison_temp.png',dpi=300)

fig,ax = plt.subplots()
ax.plot(past_co2.time,past_co2*1000000,label='Past')
ax.plot(future_co2.time,future_co2*1000000,label='RCP8.5')
ax.grid()
ax.set_ylabel('CO2 volume mixing ratio [ppm]')
ax.set_xlabel('year')
ax.legend()
fig.savefig('/storage/mirrored/homes/doensen/tmp_figs/comparison_co2.png',dpi=300)
