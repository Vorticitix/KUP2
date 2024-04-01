import pandas as pd
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
from my_tools import *
import cartopy.crs as ccrs
import xarray as xr
# %%

path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
file = "fort_30_total_era5.txt"

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep']

df = pd.read_csv(path+file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
df['lat']=ilat_to_lat(df['ilat'])
df['lon']=ilon_to_lon(df['ilon'])
df['month']= [int(float(x[4:6])) for x in df['date'].astype(str)]


ds = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro2.nc')

# %%
colorz = [(0,'#006600'),(0.125,'#33cc33'),
          (0.25,'#ffff00'),(0.375,'#ff6600'),(0.5,'#ff0000'),(0.625,'#ff6600'),(0.75,'#ffff00'),
          (0.875,'#33cc33'),(1,'#006600')]   
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom', colorz, N=1024)

def cov_plot(df):
    columns = df.columns
    fig,axz = plt.subplots(len(columns),len(columns),figsize=(18,13))
    for y in range(len(columns)):
        for x in range(len(columns)):
            ax=axz[x,y]
            if y==x:
                df[columns[x]].plot.hist(ax=ax,bins=50)
                ax.set_ylabel([])
            if y<x:
                ax.scatter(df[columns[y]],df[columns[x]],
                           marker='.',s=0.5)
    
            if x<y:
                corr = df[columns[[y,x]]].corr().values[1,0]
                ax.text(0.5 , 0.5 , '{:.3f}'.format(corr),
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                fontsize=16)
                ax.set_facecolor(cmap((corr+1)/2))
            ax.set_xlabel(columns[y])
            ax.set_ylabel(columns[x])
            ax.label_outer()
    fig = fig.tight_layout()
    return fig
        
            
# %%

df_atl = df.where(((df['lat']>20)&(df['lat']<75))&((df['lon']>300)|(df['lon']<60)))\
    .dropna(how='any')
df_iic_mean = df_atl[['iic','cz','gz','zrad','zdep','slp','precmax']]\
    .groupby('iic').mean()

fig = cov_plot(df_iic_mean)

# %%
df_atl = df.where(((df['lat']>20)&(df['lat']<75))&((df['lon']>300)|(df['lon']<60)))\
    .dropna(how='any')
df_atl_jja = df_atl.where(df['month'].isin([6,7,8])).dropna()
df_iic_min_jja = df_atl_jja[['iic','cz']].groupby('iic').min()
df_iic_max_jja = df_atl_jja[['iic','gz','zrad','zdep']].groupby('iic').max()
df_iic_ext_jja = df_iic_min_jja.join(df_iic_max_jja)\
    .reindex(columns=['cz','gz','zrad','zdep'])

fig = cov_plot(df_iic_ext_jja)
df_atl_djf = df_atl.where(df['month'].isin([12,1,2]))
df_iic_min_djf = df_atl_djf[['iic','cz']].groupby('iic').min()
df_iic_max_djf = df_atl_djf[['iic','gz','zrad','zdep']].groupby('iic').max()
df_iic_ext_djf = df_iic_min_djf.join(df_iic_max_djf)\
    .reindex(columns=['cz','gz','zrad','zdep'])
fig = cov_plot(df_iic_ext_djf)

# fig,axz = plt.subplots(2,2,figsize=(16,9))
# columns = df_iic_ext_jja.columns
# for i,ax in enumerate(axz.flat):
#     jja = df_iic_ext_jja[columns[i]]
#     djf = df_iic_ext_djf[columns[i]]
#     ax.hist(jja,bins=50,alpha=.5,label='JJA')
#     ax.hist(djf,bins=50,alpha=.5,label='DJF')
#     ax.legend(loc='upper right')
#     ax.set_title(columns[i])

# fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/hist_atl_era5.png")
# plt.close(fig)

# %%

df_med = df.where(((df['lat']>28)&(df['lat']<47))&((df['lon']>350)|(df['lon']<45)))\
    .dropna(how='any')
df_iic_mean = df_med[['iic','cz','gz','zrad','zdep','slp','precmax']]\
    .groupby('iic').mean()

fig = cov_plot(df_iic_mean)

# %%
df_med = df.where(((df['lat']>28)&(df['lat']<47))&((df['lon']>345)|(df['lon']<45)))\
    .dropna(how='any')
df_med_jja = df_med.where(df['month'].isin([6,7,8])).dropna()
df_iic_min_jja = df_med_jja[['iic','cz']].groupby('iic').min()
df_iic_max_jja = df_med_jja[['iic','gz','zrad','zdep']].groupby('iic').max()
df_iic_ext_jja = df_iic_min_jja.join(df_iic_max_jja)\
    .reindex(columns=['cz','gz','zrad','zdep'])

fig = cov_plot(df_iic_ext_jja)
df_med_djf = df_med.where(df['month'].isin([12,1,2]))
df_iic_min_djf = df_med_djf[['iic','cz']].groupby('iic').min()
df_iic_max_djf = df_med_djf[['iic','gz','zrad','zdep']].groupby('iic').max()
df_iic_ext_djf = df_iic_min_djf.join(df_iic_max_djf)\
    .reindex(columns=['cz','gz','zrad','zdep'])
fig = cov_plot(df_iic_ext_djf)

# fig,axz = plt.subplots(2,2,figsize=(16,9))
# columns = df_iic_ext_jja.columns
# for i,ax in enumerate(axz.flat):
#     jja = df_iic_ext_jja[columns[i]]
#     djf = df_iic_ext_djf[columns[i]]
#     ax.hist(jja,bins=50,alpha=.5,label='JJA')
#     ax.hist(djf,bins=50,alpha=.5,label='DJF')
#     ax.legend(loc='upper right')
#     ax.set_title(columns[i])

# fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/hist_med_era5.png")
# plt.close(fig)

# %%
path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
file = "fort_30_total.txt"

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep']

df_cesm = pd.read_csv(path+file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
df_cesm['lat']=ilat_to_lat(df_cesm['ilat'])
df_cesm['lon']=ilon_to_lon(df_cesm['ilon'])
df_cesm['month']= [int(float(x[4:6])) for x in df_cesm['date'].astype(str)]

#%%
df_cesm_atl = df_cesm.where(((df_cesm['lat']>20)&(df_cesm['lat']<75))&((df_cesm['lon']>300)|(df_cesm['lon']<60)))\
    .dropna(how='any')
df_cesm_atl_jja = df_cesm_atl.where(df_cesm['month'].isin([6,7,8])).dropna()
df_cesm_iic_min_jja = df_cesm_atl_jja[['iic','cz']].groupby('iic').min()
df_cesm_iic_max_jja = df_cesm_atl_jja[['iic','gz','zrad','zdep']].groupby('iic').max()
df_cesm_iic_ext_jja = df_cesm_iic_min_jja.join(df_cesm_iic_max_jja)\
    .reindex(columns=['cz','gz','zrad','zdep'])

#fig = cov_plot(df_cesm_iic_ext_jja)
df_cesm_atl_djf = df_cesm_atl.where(df_cesm['month'].isin([12,1,2]))
df_cesm_iic_min_djf = df_cesm_atl_djf[['iic','cz']].groupby('iic').min()
df_cesm_iic_max_djf = df_cesm_atl_djf[['iic','gz','zrad','zdep']].groupby('iic').max()
df_cesm_iic_ext_djf = df_cesm_iic_min_djf.join(df_cesm_iic_max_djf)\
    .reindex(columns=['cz','gz','zrad','zdep'])
#fig = cov_plot(df_cesm_iic_ext_djf)

fig,axz = plt.subplots(2,2,figsize=(16,9))
columns = df_cesm_iic_ext_jja.columns
for i,ax in enumerate(axz.flat):
    era5 = df_iic_ext_jja[columns[i]]
    cesm = df_cesm_iic_ext_jja[columns[i]]
    ax.hist(era5,bins=50,alpha=.5,label='ERA5')
    ax.hist(cesm,bins=50,alpha=.5,label='CESM')
    ax.legend(loc='upper right')
    ax.set_title(columns[i])



