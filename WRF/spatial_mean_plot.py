#%%
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
from glob import glob 
import sys

#%%
var = 'T2'
path = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/Long*/{}/*mean*'.format(var)
files = sorted(glob(path))
year_range=[1853,1885,1917,1949,1981]
fig,axz = plt.subplots(2,3,figsize=(16,9))
for i,file in enumerate(files):
    if i == len(files)-1:
        break
    ax = axz.flat[i]
    file_begin = files[i]
    file_end = files[i+1]
    ds_begin = xr.open_dataset(file_begin)
    ds_end = xr.open_dataset(file_end)
    da_begin = ds_begin.T2.resample(Time='1MS').mean()
    da_end = ds_end.T2.resample(Time='1MS').mean()
    da_begin = da_begin[(da_begin.Time.dt.year>=year_range[i])&
                        (da_begin.Time.dt.year<=year_range[i]+1)]
    da_end = da_end[(da_end.Time.dt.year>=year_range[i])&
                    (da_end.Time.dt.year<=year_range[i]+1)]
    ax.plot(da_begin.Time,da_begin.values,label='Ending Run')
    ax.plot(da_end.Time,da_end.values,label='Beginning Run')
    #ax.set_ylim([285,295])
    ax.grid()
    ax.set_ylim([280,300])
    if i == 0:
        box = ax.get_position()

        
        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    for label in axz.flat[i].get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')

fig.suptitle('2m Temperature Mediterranean')
fig.tight_layout()
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/figs/T2_med_mean.png',dpi=300)
plt.close(fig)
    
#%%
var = 'TSLB'
path = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/Long*/{}/*mean*'.format(var)
files = sorted(glob(path))
year_range=[1853,1885,1917,1949,1981]
fig,axz = plt.subplots(2,3,figsize=(16,9))
for i,file in enumerate(files):
    if i == len(files)-1:
        break
    ax = axz.flat[i]
    file_begin = files[i]
    file_end = files[i+1]
    ds_begin = xr.open_dataset(file_begin)
    ds_end = xr.open_dataset(file_end)
    da_begin = ds_begin.TSLB.resample(Time='1MS').mean()
    da_end = ds_end.TSLB.resample(Time='1MS').mean()
    da_begin = da_begin[(da_begin.Time.dt.year>=year_range[i])&
                        (da_begin.Time.dt.year<=year_range[i]+1)]
    da_end = da_end[(da_end.Time.dt.year>=year_range[i])&
                    (da_end.Time.dt.year<=year_range[i]+1)]
    for j in range(4):
        ax.plot(da_begin.Time,da_begin.values[:,j],label='Ending Run Layer {}'.format(j+1))
        ax.plot(da_end.Time,da_end.values[:,j],label='Beginning Run Layer {}'.format(j+1))
    #ax.set_ylim([285,295])
    ax.grid()
    ax.set_ylim([275,305])
    if i == 0:        
        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=3,prop={'size':6})
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_horizontalalignment('right')

fig.suptitle('Soil Temperature Mediterranean')
fig.tight_layout()
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/figs/TSLB_med_mean.png',dpi=300)
plt.close(fig)

#%%
path_T2 = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/Long/T2/*mean*'
files_T2 = sorted(glob(path_T2))
path_TSLB = '/storage/climatestor/PleioCEP/doensen/data/extracted/WRF/Long/TSLB/*mean*'
files_TSLB = sorted(glob(path_TSLB))

T2 = xr.open_dataset(files_T2[0]).resample(Time='1MS').mean()
TSLB = xr.open_dataset(files_TSLB[0]).resample(Time='1MS').mean()
T2_da = T2.T2; TSLB_da = TSLB.TSLB
T2_da = T2_da[T2.Time.dt.year>1870]
TSLB_da = TSLB_da[TSLB_da.Time.dt.year>1870]
fig,axz = plt.subplots(2,1,figsize=(16,9))
ax1 = axz.flat[0]
ax1.plot(T2_da.Time,T2_da.values)
ax2 = axz.flat[1]
for i in range(4):
    ax2.plot(TSLB_da.Time,TSLB_da.values[:,i])

