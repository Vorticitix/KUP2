#%%
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy
import seaborn
import cartopy.crs as ccrs
from scipy.stats import boxcox

#%%
#Load cyclone script data
model = 'CESM'
region = 'emed'
path = '/storage/climatestor/PleioCEP/doensen/data/'
if model == 'ERA5':
    file = 'fort_36_total_{}_ERA5.txt'.format(region)
elif model == 'CESM':
    file = 'fort_36_total_{}.txt'.format(region)

#Read data using pandas. Only use data from DJF and remove invalid data
df = pd.read_csv(path+file).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
if model == 'ERA5':
    df = df.where((df['year']>=1981)&(df['year']<=2010)).dropna(how='all')
elif model == 'CESM':
    df = df.where((df['year']>=3483)&(df['year']<=3512)).dropna(how='all')
df['lon'] = (df['lon'] + 180) % 360 - 180
if region == 'emed':    
    df = df.where((df['lat']<42)|(df['lon']<23))
#df = df.where((df['lat']<40))
#df = df.where(df['day']==1).dropna(how='all')
df = df.where((df['month']>=6) & (df['month']<=8)).dropna(how='all')
df['yearmonth']=np.around(df['year'].values + df['month'].values/12 -1/12,3)
df['datetime']=df['date'].astype(int).astype(str).str.zfill(8)+(df['idumi']//10000).astype(int).astype(str).str.zfill(2)
df.index = np.arange(len(df['iic']))
df_ = df[:]

#%%

#Set index where lowest slp is observed per cyclone 
idx = df_.groupby('iic')['slp'].idxmin().values[:-1]
#Set time window over which precipitation has to be summed together



df = df_[['iic','year','wsmean','wsmax','precmean','precmax','lat','lon','datetime','slp']].iloc[idx]
df['wsnorm']=(df['wsmean']-df['wsmean'].mean())/df['wsmean'].std()
df['wsmaxnorm']=(df['wsmax']-df['wsmax'].mean())/df['wsmax'].std()
df['precmean']=df['precmean'].where(df['precmean']>1)
df['precmax']=df['precmax'].where(df['precmax']>1)
df['prectrans'] = df['precmean']**(1/3)
df['prectransmax'] = df['precmax']**(1/2)
df['precnorm']=(df['prectrans']-df['prectrans'].mean())/df['prectrans'].std()
df['precmaxnorm']=(df['prectransmax']-df['prectransmax'].mean())/df['prectransmax'].std()
#%%
#prec_ext = df.where(df['preccum']>=df['preccum'].quantile(q=0.1)).dropna().sort_values('preccum',ascending=False)
#ws_ext = df.where(df['wsmean']>=df['wsmean'].quantile(q=0.1)).dropna().sort_values('wsmean',ascending=False)

prec_ext = df.where(df['precmean']>=df['precmean'].quantile(q=0.5)).dropna().sort_values('precmean',ascending=False)
ws_ext = df.where(df['wsmean']>=df['wsmean'].quantile(q=0.5)).dropna().sort_values('wsmean',ascending=False)
cmpnd_ext = np.intersect1d(prec_ext['iic'],ws_ext['iic'])
cmpnd_df = df.where(df['iic'].isin(cmpnd_ext)).dropna()
cmpnd_df['euc_dist'] = (cmpnd_df['wsnorm']**2+cmpnd_df['precnorm']**2)**.5
cmpnd_df = cmpnd_df.sort_values('euc_dist',ascending=False)
names = ['precmean','wsmean']

precmax_ext = df.where(df['precmax']>=df['precmax'].quantile(q=0.6)).dropna().sort_values('precmax',ascending=False)
wsmax_ext = df.where(df['wsmax']>=df['wsmax'].quantile(q=0.6)).dropna().sort_values('wsmax',ascending=False)
cmpndmax_ext = np.intersect1d(precmax_ext['iic'],wsmax_ext['iic'])
cmpndmax_df = df.where(df['iic'].isin(cmpndmax_ext)).dropna()
cmpndmax_df['euc_dist'] = (cmpndmax_df['wsmaxnorm']**2+cmpndmax_df['precmaxnorm']**2)**.5
cmpndmax_df = cmpndmax_df.sort_values('euc_dist',ascending=False)
names = ['precmean','wsmean']

if model == 'ERA5':
    prec_ext.to_csv(path+'cyclone_intensity/prec_ext_{}_ERA5_JJA.csv'.format(region))
    ws_ext.to_csv(path+'cyclone_intensity/ws_ext_{}_ERA5_JJA.csv'.format(region))
    cmpnd_df.to_csv(path+'cyclone_intensity/compound_ext_{}_ERA5_JJA.csv'.format(region))

    precmax_ext.to_csv(path+'cyclone_intensity/precmax_ext_{}_ERA5_JJA.csv'.format(region))
    wsmax_ext.to_csv(path+'cyclone_intensity/wsmax_ext_{}_ERA5_JJA.csv'.format(region))
    cmpndmax_df.to_csv(path+'cyclone_intensity/compoundmax_ext_{}_ERA5_JJA.csv'.format(region))
elif model == 'CESM':
    prec_ext.to_csv(path+'cyclone_intensity/prec_ext_{}_CESM_modern_JJA.csv'.format(region))
    ws_ext.to_csv(path+'cyclone_intensity/ws_ext_{}_CESM_modern_JJA.csv'.format(region))
    cmpnd_df.to_csv(path+'cyclone_intensity/compound_ext_{}_CESM_modern_JJA.csv'.format(region))

    precmax_ext.to_csv(path+'cyclone_intensity/precmax_ext_{}_CESM_modern_JJA.csv'.format(region))
    wsmax_ext.to_csv(path+'cyclone_intensity/wsmax_ext_{}_CESM_modern_JJA.csv'.format(region))
    cmpndmax_df.to_csv(path+'cyclone_intensity/compoundmax_ext_{}_CESM_modern_JJA.csv'.format(region))

da_tot = []
for j,var_ in enumerate([prec_ext,ws_ext]):
    da = xr.DataArray(dims=['lat','lon'],
                  coords={'lat':sorted(df_['lat'].unique(),reverse=True),
                          'lon':sorted(df_['lon'].unique())})
    var = var_.groupby(['lat','lon'])[names[j]].max()
    midx = var.index
    for lat,lon in midx:
        da.loc[lat,lon]=var.loc[lat].loc[lon]
    da.coords['lon'] = (da.coords['lon'] + 180) % 360 - 180
    da = da.sortby(da.lon)
    da_tot.append(da)

#%%
fig,axz = plt.subplots(2,1,figsize=(12,10),subplot_kw={'projection':ccrs.PlateCarree()})
for i,ax in enumerate(axz.flat):
    da=da_tot[j]
    #if i ==0:
    im=ax.pcolormesh(da.lon,da.lat,da)
    #else:
    #    ax.pcolormesh(da.lon,da.lat,da)
    ax.coastlines()
    fig.colorbar(im,ax=ax)
    ax.set_title(names[j])
#cbar_ax = fig.add_axes([0.05, 0.065, 0.9, 0.015])
#fig.colorbar(im,cax=cbar_ax,orientation='horizontal')

# %%