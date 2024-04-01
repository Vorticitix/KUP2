#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:59:39 2022

@author: doensen
"""

import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import cftime
from scipy import stats
import sys

#%%


path = "/storage/climatestor/PleioCEP/doensen/data/"
# file_past = 'fort_36_total_med.txt'
file_med = 'fort_36_total_med.txt'
file_rcp = 'fort_36_total_med_RCP85.txt'


med = pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
rcp = pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
rcp = rcp.where(rcp['year']>=3515).dropna()
df = pd.concat([med,rcp])
df = df.where(df['year']>=100).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
#df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
df['yearmonth']=df['year'] + df['month']/12 -1/12 
df['cen'] = (df['year']-1502)//100


#%%
#df_groupby= df.groupby('yearseason')
df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
df_gb_min = df[['slp','iic']].groupby('iic').min()
df_gb_yas = df[['yearmonth','iic','cen']].groupby('iic').first()
df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                            &(df_gb_iic['zrad']!=0)
                            &(df_gb_iic['zdep']!=0)).dropna()


#%%
radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
#radius_cen = df_gb_iic[['zrad','cen']].groupby('cen').median()
# fig,ax = plt.subplots()
# ax.plot(radius.index-1502,radius.values)
#ax.plot(radius_cen.index*100,radius_cen.values)
# ax.grid()
#ax.set_xticks(radius.index)
#ax.set_ylim([0.15,0.25])

#%%
slp = df_gb_iic[['slp','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
# fig,ax = plt.subplots()
# ax.plot(slp.index-1502,slp.values)
# ax.grid()

#%%
depth = df_gb_iic[['zdep','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
# fig,ax = plt.subplots()
# ax.plot(depth.index-1502,depth.values)
#%%
prec = df_gb_iic[['precmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
#prec = prec[prec.index<200].dropna()
# fig,ax = plt.subplots()
# ax.plot(prec.index,prec.values)
# ax.grid()
#%%
precc = df_gb_iic[['preccmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
# fig,ax = plt.subplots()
# ax.plot(precc.index,precc.values)
# ax.grid()
#%%
precl = df_gb_iic[['preclmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
# fig,ax = plt.subplots()
# ax.plot(precl.index,precl.values)
# ax.grid()
#%%
ws = df_gb_iic[['wsmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
#ws = ws[ws.index<200].dropna()
# fig,ax = plt.subplots()
# ax.plot(ws.index,ws.values)



#%%
full_time = radius.index[(radius.index>=5)&(radius.index<=3352)]
T850 = xr.open_dataset(path+'extracted/T850/T850_0005_3352_anom_med_mean.nc').squeeze()
# T850 = T850.assign_coords(lon=(((T850.lon.values + 180) % 360) - 180)).sortby('lon')
time_arr = T850.time.values
time = np.array([x.year + x.month/12-1/12 for x in time_arr])
T850 = T850.assign_coords(time=time)
T850 = T850.reindex(time=full_time)

# T850_med = T850.sel(lat=slice(47,28),lon=slice(-10,40))
# weights_med = np.cos(np.deg2rad(T850_med.lat))
# T850_med_weighted = T850_med.weighted(weights_med)
# T850_med_mean = T850_med_weighted.mean(('lat','lon'))


#%%
seasons_dic = {1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',6:'JJA',
               7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON',12:'DJF'}
names = ['Mediterranean [47 °N - 28 °N, 10 °W - 40 °E]',
         'Atlantic [70 °N - 30 °N, 65 °W - 10 °W]']
areas  = ['med','atl']
quants = [.5,.9,.99]
for i,area in enumerate(areas):
    path = "/storage/climatestor/PleioCEP/doensen/data/"
    # file_past = 'fort_36_total_med.txt'
    file = 'fort_36_total_{}.txt'.format(area)

    #df_past = pd.read_csv(path+file_past,index_col='index').drop(['Unnamed: 0','level_0'],axis=1).reset_index(drop=True)
    #df_fut = pd.read_csv(path+file_fut,index_col='index').drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    #df_fut['iic'] =df_fut['iic']+100000
    #df = pd.concat([df_past,df_fut]).reset_index()
    df = pd.read_csv(path+file).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    df = df.where(df['year']>=100).dropna()
    df = df.where((df['month']>=12) | (df['month']<=2)).dropna()
    df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
    df['yearmonth']=df['year'] + df['month']/12 -1/12 
    df['cen'] = (df['year']-1502)//100
    if area == 'atl':
        df = df.where(df['lon']<=350).dropna()
        
    
    df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
    df_gb_min = df[['slp','iic']].groupby('iic').min()
    df_gb_yas = df[['yearmonth','iic','cen']].groupby('iic').first()
    df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
    df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                                &(df_gb_iic['zrad']!=0)
                                &(df_gb_iic['zdep']!=0)).dropna()
    
    for j,quant in enumerate(quants):
        radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
        slp = df_gb_iic[['slp','yearmonth']].groupby('yearmonth').quantile(1-quant).rolling(90).mean()
        depth = df_gb_iic[['zdep','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
        prec = df_gb_iic[['precmean','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
        ws = df_gb_iic[['wsmean','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
        fig,axz=plt.subplots(6,1,figsize=(16,12))
        index = radius.index - 1502
        full_time = radius.index
        T850 = xr.open_dataset(path+'extracted/T850/T850_0005_3352_anom_med_mean.nc').squeeze()
        # T850 = T850.assign_coords(lon=(((T850.lon.values + 180) % 360) - 180)).sortby('lon')
        time_arr = T850.time.values
        time = np.array([x.year + x.month/12-1/12 for x in time_arr])
        T850 = T850.assign_coords(time=time)
        T850 = T850.reindex(time=full_time).T.rolling(time=90).mean()
        axz[0].plot(index,T850.values)
        axz[0].set_ylabel('T850 Anomaly [K]')        
        axz[1].plot(index,radius.values)
        axz[1].set_ylabel('Radius [1000 km]')
        axz[2].plot(index,depth.values)
        axz[2].set_ylabel('Depth [gpm]')
        axz[3].plot(index,slp.values)
        axz[3].set_ylabel('SLP [hPa]')
        axz[4].plot(index,prec.values)
        axz[4].set_ylabel('Precipitation [mm]')
        axz[5].plot(index,ws.values)
        axz[5].set_ylabel('Wind Speed [m/s]')
        
        for ax in axz.flat:
            ax.grid()
            ax.label_outer()
            ax.set_xlim([-1500,2100])
        
        sys.exit()
        fig.subplots_adjust(left=.05,right=.99,top=.95,bottom=.035,hspace=.075)
        # fig.suptitle('Cyclones {}: {}th Percentile Cyclone Stastistics DJF - 30 year running mean'\
                    #  .format(names[i],int(quant*100)))
        # fig.savefig(path+'figs/{}_{}th_3500.png'.format(area,int(quant*100)),dpi=300)
    

#%%
"""
Create scatter plot between temperature anomalies and cyclone related variables
"""
radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').median()
full_time = radius.index[(radius.index>=5)&(radius.index<=3352)]
T850 = xr.open_dataset(path+'extracted/T850_fldmean_med_anom.nc').squeeze()
time_arr = T850.time.values
time = np.array([x.year + x.month/12-1/12 for x in time_arr])
T850 = T850.assign_coords(time=time)
T850 = T850.reindex(time=full_time)


def plot_corr(x,y):
    x = x.T.values; y = np.squeeze(y.values)
    fig,ax=plt.subplots(figsize=(12,8))
    ax.scatter(x, y)
    idx = np.isfinite(x) & np.isfinite(y)
    x = x[idx]; y = y[idx]
    z = np.polyfit(x,y,1)
    d = np.poly1d(z)
    r,p = stats.pearsonr(x,y)
    r = round(r,3); p = round(p,3)
    ax.text(0.99,1.03,d, va='top',ha="right",
            color="k",transform = ax.transAxes,size=15)
    ax.text(0.99,0.95,'r = {}'.format(r), va='top',ha="right",
            color="k",transform = ax.transAxes,size=15)
    ax.text(0.99,0.91,'p = {}'.format(p), va='top',ha="right",
            color="k",transform = ax.transAxes,size=15)
    ax.plot(x,d(x),color='r')
    return fig,ax
    

#%%
slp = df_gb_iic[['slp','yearmonth']].groupby('yearmonth').quantile(.01)
slp = slp[(slp.index>=2000)&(slp.index<=3352)]
fig,ax = plot_corr(T850,slp)
ax.set_ylabel('Sea Level Pressure [hPa]')
ax.set_xlabel('850 hPa Temperature Anomaly [K]')
fig.suptitle('Correlation 1st Percentile Monthly Pressure vs 850 hPa Temperature montly Anomaly')
fig.savefig(path+'figs/corr_slp.png',dpi=300)
#%%
prec = df_gb_iic[['precmean','yearmonth']].groupby('yearmonth').quantile(.99)
prec = prec[(prec.index>=2000)&(prec.index<=3352)]
fig,ax = plot_corr(T850,prec)
ax.set_ylabel('Cyclone Precipitation [mm]')
ax.set_xlabel('850 hPa Temperature Anomaly [K]')
fig.suptitle('Correlation 99st Percentile Monthly Cyclone Precip. vs 850 hPa Temperature montly Anomaly')
fig.savefig(path+'figs/corr_prec.png',dpi=300)


#%%
precc = df_gb_iic[['preccmean','yearmonth']].groupby('yearmonth').quantile(.99)
precc = precc[(precc.index>=2000)&(precc.index<=3352)]
plot_corr(T850,precc)


#%%
precl = df_gb_iic[['preclmean','yearmonth']].groupby('yearmonth').quantile(.99)
precl = precl[(precl.index>=2000)&(precl.index<=3352)]
plot_corr(T850,precl)


#%%
depth = df_gb_iic[['zdep','yearmonth']].groupby('yearmonth').quantile(.99)
depth = depth[(depth.index>=2000)&(depth.index<=3352)]
plot_corr(T850, depth)
#%%
radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(.99)
radius = radius[(radius.index>=2000)&(radius.index<=3352)]
plot_corr(T850, radius)

#%%
radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(.99)
radius = radius[(radius.index>=2000)&(radius.index<=3352)]
fig,ax=plt.subplots()
ax.scatter(T850.T.values, radius.values.squeeze())
#%%
ws = df_gb_iic[['wsmean','yearmonth']].groupby('yearmonth').quantile(.99)
ws = ws[(ws.index>=2000)&(ws.index<=3352)]
plot_corr(T850, ws)

