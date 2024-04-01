import pandas as pd
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
from my_tools import *
import cartopy.crs as ccrs
import xarray as xr
# %%

path = "/storage/climatestor/PleioCEP/doensen/data/cyclone_0005_0604/"
file = "test.txt"

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmean',
         'precmax','preccmean','preclmean','wsmean','wsmax']

df = pd.read_csv(path+file,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
df['lat']=ilat_to_lat(df['ilat'])
df['lon']=ilon_to_lon(df['ilon'])
df['month']= [int(float(x[4:6])) for x in df['date'].astype(str).str.zfill(8)]

ds = xr.open_dataset('/storage/climatestor/PleioCEP/doensen/data/extracted/oro.nc')

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
df_iic_min_jja = df_atl_jja[['iic','cz','slp']].groupby('iic').min()
df_iic_max_jja = df_atl_jja[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext_jja = df_iic_min_jja.join(df_iic_max_jja)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax'])

fig = cov_plot(df_iic_ext_jja)

df_atl_djf = df_atl.where(df['month'].isin([12,1,2]))
df_iic_min_djf = df_atl_djf[['iic','cz','slp']].groupby('iic').min()
df_iic_max_djf = df_atl_djf[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext_djf = df_iic_min_djf.join(df_iic_max_djf)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax'])
fig = cov_plot(df_iic_ext_djf)

# fig,axz = plt.subplots(3,2,figsize=(16,9))
# columns = df_iic_ext_jja.columns
# for i,ax in enumerate(axz.flat):
#     jja = df_iic_ext_jja[columns[i]]
#     djf = df_iic_ext_djf[columns[i]]
#     ax.hist(jja,bins=50,alpha=.5,label='JJA')
#     ax.hist(djf,bins=50,alpha=.5,label='DJF')
#     ax.legend(loc='upper right')
#     ax.set_title(columns[i])
# fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/hist_atl_cesm.png")
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
df_iic_min_jja = df_med_jja[['iic','cz','slp']].groupby('iic').min()
df_iic_max_jja = df_med_jja[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext_jja = df_iic_min_jja.join(df_iic_max_jja)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax'])

fig = cov_plot(df_iic_ext_jja)
df_med_djf = df_med.where(df['month'].isin([12,1,2]))
df_iic_min_djf = df_med_djf[['iic','cz','slp']].groupby('iic').min()
df_iic_max_djf = df_med_djf[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext_djf = df_iic_min_djf.join(df_iic_max_djf)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax'])
fig = cov_plot(df_iic_ext_djf)

# fig,axz = plt.subplots(3,2,figsize=(16,9))
# columns = df_iic_ext_jja.columns
# for i,ax in enumerate(axz.flat):
#     jja = df_iic_ext_jja[columns[i]]
#     djf = df_iic_ext_djf[columns[i]]
#     ax.hist(jja,bins=50,alpha=.5,label='JJA')
#     ax.hist(djf,bins=50,alpha=.5,label='DJF')
#     ax.legend(loc='upper right')
#     ax.set_title(columns[i])
# fig.savefig("/storage/climatestor/PleioCEP/doensen/figs/hist_med_cesm.png")
# plt.close(fig)

# %%
# df_sub=df.where((df['lat']>35)&(df['lat']<52)&\
#                   (df['lon']>294)&(df['lon']<308)).dropna()
df_sub = df.dropna()
df_count = df_sub.groupby(['lat','lon']).count()['iic'].unstack(fill_value=0)
# df_count = df_count.where(df_count<100,100)
# df_count[0] = df_count[0].where(df_count[0]<1000)
lats_uq = df_count.index
lons_uq = df_count.columns
# ds_sub = ds.sel(lat=lats_uq,method='nearest')['PHIS']<9810
# df_count = df_count.where(ds_sub.values,np.nan)
# df_count_stacked = df_count.stack(dropna=False).reset_index()[::-1]
# nan_bool = np.isnan(df_count_stacked[0])
# idx_nan = df_count_stacked[nan_bool][['lat','lon']]


# idx = [not np.isin(df.iloc[i]['lat'],idx_nan['lat'])\
#         &np.isin(df.iloc[i]['lon'],idx_nan['lon']) for i in df.index.values]

# df_non_oro = df[idx]

# for i in df_non_oro.index.values:
#     boolean = np.isin(df_non_oro.iloc[i]['lat'],idx_nan['lat'])&\
#     np.isin(df_non_oro.iloc[i]['lon'],idx_nan['lon'])
#     if boolean:
#         df_non_oro.iloc[i] = np.nan
        
# [df_non_oro.iloc[i]=np.nan for i in df.index.values if (np.isin(df_non_oro.iloc[i]['lat'],idx_nan['lat'])&np.isin(df_non_oro.iloc[i]['lon'],idx_nan['lon']))]
# [df_non_oro.iloc[i]=np.nan for i in df.index.values]# if (np.isin(df_non_oro.iloc[i]['lat'],idx_nan['lat'])&np.isin(df_non_oro.iloc[i]['lon'],idx_nan['lon']))]

# df_non_oro = df_non_oro.dropna(how='all')
        
                
            
# lats = np.unique(lats)[::-1]
# lons = np.unique(lons)
# counts = counts.reshape(int(len(counts)/144),144)[::-1,:]

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.Orthographic(0, 90)})
ax.gridlines()
ax.pcolormesh(lons_uq,lats_uq,df_count.values,transform=ccrs.PlateCarree(),
              cmap='hot_r',vmin=0,vmax=50)
#ax.contour(ds.lon.values,ds.lat.values,(ds['PHIS']>9810).values,transform=ccrs.PlateCarree(),
#               cmap='Greys')
ax.coastlines()
# %%
df_iic_min = df_sub[['iic','cz','slp']].groupby('iic').min()
df_iic_max = df_sub[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext = df_iic_min.join(df_iic_max)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax','precmean'])
df_sub_cz = df_iic_ext.where((df_iic_ext['cz']>60)&(df_iic_ext['cz']<70)).dropna()
fig = cov_plot(df_sub_cz)


# %%

season_dic = {"MAM":[3,4,5],
              "JJA":[6,7,8],
              "SON":[9,10,11],
              "DJF":[12,1,2]}
df_sub_month = df_med[:]
df_sub_month['month']=[int(str(int(date))[4:6]) for date in df_sub_month['date']]

fig,axz = plt.subplots(2,2,figsize=(16,9),
                       subplot_kw={'projection':ccrs.Orthographic(20, 35)})
for k,key in enumerate(list(season_dic.keys())):
    idx = np.isin(df_sub_month['month'],season_dic[key])
    df_sub_season = df_sub_month[idx]
    # df_sub_season=df_sub_season.where(((df_sub_season['lat']>28)&(df_sub_season['lat']<47))&\
    #                  ((df_sub_season['lon']>350)|(df_sub_season['lon']<40))).dropna()
    df_count_med = df_sub_season.groupby(['lat','lon']).count()['iic'].unstack(fill_value=0)
    # df_count[0] = df_count[0].where(df_count[0]<100,100)
    # df_count[0] = df_count[0].where(df_count[0]<1000)
    lats_uq = df_count_med.index
    lons_uq = df_count_med.columns
    
    ax = axz.flat[k]
    ax.gridlines()
    ax.pcolormesh(lons_uq,lats_uq,df_count_med.values,transform=ccrs.PlateCarree(),
                  cmap='hot_r')
    ax.coastlines()
    ax.set_title(key)

    
# %%
lat = 31.263157894736842
lon = 15.000000

df_date = df.where((df['lat']==lat)&(df['lon']==lon)).dropna()

# df_date = df_date[['iic','date','agetot','cz','gz','zrad','zdep','slp','precmean']]
# df_date.to_csv("../data/cyclone/test.csv")

nc_file = "/storage/climatestor/PleioCEP/kim/cesm1.2_data/archive/TRANS.3501BP_3/atm/hist/TRANS.3501BP.cam.h1.3481-01-01-00000.nc"
ds = xr.open_dataset(nc_file)
z1000_file = "/storage/climatestor/PleioCEP/doensen/data/interpolated/TRANS.3501BP.cam.h1.3481-01-01.z1000_alt.nc"
Z1000 = xr.open_dataset(z1000_file).sel(time=ds.time[1146:1154]).sel(lat=slice(28,47),lon=lons_uq)
ds_sub = ds.sel(time=ds.time[161:170])
# PSL = ds_sub.PSL.sel(lat=slice(30,18),lon=slice(270,300))
fig,axz = plt.subplots(2,4,figsize=(16,9),
                           subplot_kw={'projection':ccrs.Orthographic(20, 35)})
for i,time in enumerate(Z1000.time):
    if i==10:
        break
    ax = axz.flat[i]
    var = Z1000.sel(time=time).squeeze().Z3
    ax.pcolormesh(ds.lon.values,ds["PHIS"]>1000,cmap='greys')
    im = ax.pcolormesh(var.lon.values,var.lat.values,var.values,transform=ccrs.PlateCarree())
    ax.coastlines()


fig.colorbar(im)
fig.tight_layout()




