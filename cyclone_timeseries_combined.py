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

locs = ['med','wmed','cmed','emed']
period = 'past'


for track_loc in locs:
    print(track_loc)
    path = "/storage/climatestor/PleioCEP/doensen/data/"
    # file_past = 'fort_36_total_med.txt'
    file_med = 'fort_36_total_{}_full.txt'.format(track_loc)


    df=pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    #df_rcp=pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    # med = pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    # rcp = pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    # rcp = rcp.where(rcp['year']>=3515).dropna()
    # df = pd.concat([med,rcp])
    df = df.where(df['year']>=5).dropna(how='all')
    df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
    #df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
    df['yearmonth']=np.around(df['year'].values + df['month'].values/12 -1/12,3)
    #df['cen'] = (df['year']-1502)//100



    #df_groupby= df.groupby('yearseason')
    df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
    df_gb_min = df[['slp','iic','agetot']].groupby('iic').min()
    df_gb_yas = df[['yearmonth','iic']].groupby('iic').first()
    df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
    df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                                &(df_gb_iic['zrad']!=0)
                                &(df_gb_iic['zdep']!=0)).dropna(how='all')

    full_time = np.around(np.sort(np.concatenate([np.arange(5,3602),
                                            np.arange(5+1/12,3602+1/12),
                                           np.arange(6-1/12,3603-1/12)])),3)
    agetot = df_gb_iic[['agetot','yearmonth']]
    agetot = agetot.groupby('yearmonth').sum().reindex(full_time).rolling(90,min_periods=30).mean()
    agetot_mean = agetot[agetot.index<=3352].dropna().mean().values[0]
    agetot_std = agetot[agetot.index<=3352].dropna().std().values[0]
    agetot_plus_1std = agetot_mean + agetot_std; agetot_min_1std = agetot_mean - agetot_std
    agetot_plus_2std = agetot_mean + 2*agetot_std; agetot_min_2std = agetot_mean - 2*agetot_std
    agetot_plus_3std = agetot_mean + 3*agetot_std; agetot_min_3std = agetot_mean - 3*agetot_std


    radius = df_gb_iic[['zrad','yearmonth']]
    radius=radius.groupby('yearmonth').quantile(.99).reindex(full_time).rolling(90,min_periods=30).mean()
    radius_mean = radius[radius.index<=3352].dropna().mean().values[0]
    radius_std = radius[radius.index<=3352].dropna().std().values[0]
    radius_plus_1std = radius_mean + radius_std; radius_min_1std = radius_mean - radius_std
    radius_plus_2std = radius_mean + 2*radius_std; radius_min_2std = radius_mean - 2*radius_std
    radius_plus_3std = radius_mean + 3*radius_std; radius_min_3std = radius_mean - 3*radius_std

    #radius_cen = df_gb_iic[['zrad','cen']].groupby('cen').median()
    # fig,ax = plt.subplots()
    # ax.plot(radius.index-1502,radius.values)
    #ax.plot(radius_cen.index*100,radius_cen.values)
    # ax.grid()
    #ax.set_xticks(radius.index)
    #ax.set_ylim([0.15,0.25])

    slp = df_gb_iic[['slp','yearmonth']]
    slp = slp.groupby('yearmonth').quantile(1-.99).reindex(full_time).rolling(90,min_periods=30).mean()
    slp_mean = slp[slp.index<=3352].dropna().mean().values[0]
    slp_std = slp[slp.index<=3352].dropna().std().values[0]
    slp_plus_1std = slp_mean + slp_std; slp_min_1std = slp_mean - slp_std
    slp_plus_2std = slp_mean + 2*slp_std; slp_min_2std = slp_mean - 2*slp_std
    slp_plus_3std = slp_mean + 3*slp_std; slp_min_3std = slp_mean - 3*slp_std

    # fig,ax = plt.subplots()
    # ax.plot(slp.index-1502,slp.values)
    # ax.grid()

    depth = df_gb_iic[['zdep','yearmonth']]
    depth = depth.groupby('yearmonth').quantile(.99).reindex(full_time).rolling(90,min_periods=30).mean()
    depth_mean = depth[depth.index<=3352].dropna().mean().values[0]
    depth_std = depth[depth.index<=3352].dropna().std().values[0]
    depth_plus_1std = depth_mean + depth_std; depth_min_1std = depth_mean - depth_std
    depth_plus_2std = depth_mean + 2*depth_std; depth_min_2std = depth_mean - 2*depth_std
    depth_plus_3std = depth_mean + 3*depth_std; depth_min_3std = depth_mean - 3*depth_std

    # fig,ax = plt.subplots()
    # ax.plot(depth.index-1502,depth.values)

    prec = df_gb_iic[['precmean','yearmonth']]
    prec = prec.groupby('yearmonth').quantile(.99).reindex(full_time).rolling(90,min_periods=30).mean()
    prec_mean = prec[prec.index<=3352].dropna().mean().values[0]
    prec_std = prec[prec.index<=3352].dropna().std().values[0]
    prec_plus_1std = prec_mean + prec_std; prec_min_1std = prec_mean - prec_std
    prec_plus_2std = prec_mean + 2*prec_std; prec_min_2std = prec_mean - 2*prec_std
    prec_plus_3std = prec_mean + 3*prec_std; prec_min_3std = prec_mean - 3*prec_std

    #prec = prec[prec.index<200].dropna()
    # fig,ax = plt.subplots()
    # ax.plot(prec.index,prec.values)
    # ax.grid()
    precc = df_gb_iic[['preccmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
    # fig,ax = plt.subplots()
    # ax.plot(precc.index,precc.values)
    # ax.grid()
    precl = df_gb_iic[['preclmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
    # fig,ax = plt.subplots()
    # ax.plot(precl.index,precl.values)
    # ax.grid()
    ws = df_gb_iic[['wsmean','yearmonth']]
    ws = ws.groupby('yearmonth').quantile(.99).reindex(full_time).rolling(90,min_periods=30).mean()
    ws_mean = ws[ws.index<=3352].dropna().mean().values[0]
    ws_std = ws[ws.index<=3352].dropna().std().values[0]
    ws_plus_1std = ws_mean + ws_std; ws_min_1std = ws_mean - ws_std
    ws_plus_2std = ws_mean + 2*ws_std; ws_min_2std = ws_mean - 2*ws_std
    ws_plus_3std = ws_mean + 3*ws_std; ws_min_3std = ws_mean - 3*ws_std




    
    T850 = xr.open_dataset(path+'extracted/T850/T850_0005_3601_anom_full_{}.nc'.format(track_loc)).squeeze()
    #T850_RCP85 = xr.open_dataset(path+'extracted/T850/T850_RCP85_3507_3601_anom_{}.nc'.format(temp_loc)).squeeze()
    
    time_arr = T850.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    T850 = T850.assign_coords(time=time)
    T850 = T850.reindex(time=full_time)
    T850 = T850.T.rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
    T850_std = float(T850.sel(time=slice(5,3352)).std())
    T850_plus_1std = T850_std; T850_min_1std = -T850_std
    T850_plus_2std = 2*T850_std; T850_min_2std = -2*T850_std
    T850_plus_3std = 3*T850_std; T850_min_3std = -3*T850_std
    
    RWP = xr.open_dataset(path+'extracted/RWP/RWP_all_monmean_{}.nc'.format(track_loc)).squeeze()    
    time_arr = RWP.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    u, c = np.unique(time, return_counts=True)
    sim  = u[c == 1]
    idx = np.where(np.isin(time,sim)==True)[0]
    RWP = RWP.assign_coords(time=time)
    RWP = RWP.isel(time=idx)
    RWP = RWP.reindex(time=full_time)
    RWP = RWP.V.rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
    RWP_mean = float(RWP.sel(time=slice(0,3352)).mean())
    RWP_std = float(RWP.sel(time=slice(0,3352)).std())
    RWP_plus_1std = RWP_mean + RWP_std; RWP_min_1std = RWP_mean - RWP_std
    RWP_plus_2std = RWP_mean + 2*RWP_std; RWP_min_2std = RWP_mean - 2*RWP_std
    RWP_plus_3std = RWP_mean + 3*RWP_std; RWP_min_3std = RWP_mean - 3*slp_std
    
    
    NAOI = xr.open_dataset(path+'extracted/try_PSL/NAOI_full.nc').squeeze()  
    NAOI = NAOI.rename({'PSL':"NAOI"})
    time_arr = NAOI.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    u, c = np.unique(time, return_counts=True)
    sim  = u[c == 1]
    idx = np.where(np.isin(time,sim)==True)[0]
    NAOI = NAOI.assign_coords(time=time)
    NAOI = NAOI.isel(time=idx)
    NAOI = NAOI.reindex(time=full_time)
    NAOI = NAOI.NAOI.rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
    NAOI_mean = float(NAOI.sel(time=slice(0,3352)).mean())
    NAOI_std = float(NAOI.sel(time=slice(0,3352)).std())
    NAOI_plus_1std = NAOI_std; NAOI_min_1std = -NAOI_std
    NAOI_plus_2std = 2*NAOI_std; NAOI_min_2std = -2*NAOI_std
    NAOI_plus_3std = 3*NAOI_std; NAOI_min_3std = -3*NAOI_std
    
    print(full_time)
    
    seasons_dic = {1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',6:'JJA',
                7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON',12:'DJF'}
    names = ['Mediterranean [47 °N - 28 °N, 10 °W - 40 °E]',
            'Atlantic [30 °N - 30 °N, 65 °W - 10 °W]']
    path = "/storage/climatestor/PleioCEP/doensen/data/"
    # file_past = 'fort_36_total_med.txt'
    
    #df_past = pd.read_csv(path+file_past,index_col='index').drop(['Unnamed: 0','level_0'],axis=1).reset_index(drop=True)
    #df_fut = pd.read_csv(path+file_fut,index_col='index').drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    #df_fut['iic'] =df_fut['iic']+100000
    #df = pd.concat([df_past,df_fut]).reset_index()
    # df = pd.read_csv(path+file).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
    # df = df.where(df['year']>=100).dropna(how='all')
    # df = df.where((df['month']>=12) | (df['month']<=2)).dropna()
    # df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
    # df['yearmonth']=df['year'] + df['month']/12 -1/12 
    # df['cen'] = (df['year']-1502)//100
    
    # df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
    # df_gb_min = df[['slp','iic']].groupby('iic').min()
    # df_gb_yas = df[['yearmonth','iic','cen']].groupby('iic').first()
    # df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
    # df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
    #                             &(df_gb_iic['zrad']!=0)
    #                             &(df_gb_iic['zdep']!=0)).dropna()
    
    
    #radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
    #slp = df_gb_iic[['slp','yearmonth']].groupby('yearmonth').quantile(1-quant).rolling(90).mean()
    #depth = df_gb_iic[['zdep','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
    #prec = df_gb_iic[['precmean','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
    #ws = df_gb_iic[['wsmean','yearmonth']].groupby('yearmonth').quantile(quant).rolling(90).mean()
    fig,axz=plt.subplots(8,1,figsize=(17,10))
    plt.rcParams['axes.labelsize'] = 10
    index = radius.index - 1502
    #full_time = radius.index
    #T850 = xr.open_dataset(path+'extracted/T850/T850_0005_3352_anom_med_mean.nc').squeeze()
    # T850 = T850.assign_coords(lon=(((T850.lon.values + 180) % 360) - 180)).sortby('lon')
    #time_arr = T850.time.values
    #time = np.array([x.year + x.month/12-1/12 for x in time_arr])
    #axz[0].plot(index,T850.values)
    axz[0].set_ylabel('T850 [K]')
    axz[0].set_ylim([-.5,.5])
    axz[0].fill_between(index, 0, T850.values, where=T850.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[0].fill_between(index, 0, T850.values, where=T850.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    #axz[0].fill_between(index, 0, T850_plus_1std, color='darkornage', alpha=0.7)
    axz[0].fill_between(index, 0, T850_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[0].fill_between(index, 0, T850_plus_3std, color='lightgrey', alpha=0.1)
    #axz[0].fill_between(index, 0, T850_min_1std, color='lightgrey', alpha=0.7)
    axz[0].fill_between(index, 0, T850_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[0].fill_between(index, 0, T850_min_3std, color='lightgrey', alpha=0.1)
    axz[0].axhline(y=0,color='k')
    axz[0].axvline(x=2005,color='g')        
    axz[1].set_ylabel('RWP Envelope [m/s]')
    axz[1].fill_between(index, RWP_mean, RWP.values, where=RWP.values >= RWP_mean, facecolor='red', interpolate=True,zorder=2)
    axz[1].fill_between(index, RWP_mean, RWP.values, where=RWP.values <= RWP_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[1].fill_between(index, RWP_mean, RWP_plus_1std, color='lightgrey', alpha=0.7)
    axz[1].fill_between(index, RWP_mean, RWP_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[1].fill_between(index, RWP_mean, RWP_plus_3std, color='lightgrey', alpha=0.1)
    #axz[1].fill_between(index, RWP_mean, RWP_min_1std, color='lightgrey', alpha=0.7)
    axz[1].fill_between(index, RWP_mean, RWP_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[1].fill_between(index, RWP_mean, RWP_min_3std, color='lightgrey', alpha=0.1)
    axz[1].axhline(y=RWP_mean,color='k')
    axz[1].axvline(x=2005,color='g')      
    axz[2].set_ylabel('NAOI [-]')
    axz[2].fill_between(index, 0, NAOI.values, where=NAOI.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[2].fill_between(index, 0, NAOI.values, where=NAOI.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    #axz[2].fill_between(index, 0, NAOI_plus_1std, color='lightgrey', alpha=0.7)
    axz[2].fill_between(index, 0, NAOI_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[2].fill_between(index, 0, NAOI_plus_3std, color='lightgrey', alpha=0.1)
    #axz[2].fill_between(index, 0, NAOI_min_1std, color='lightgrey', alpha=0.7)
    axz[2].fill_between(index, 0, NAOI_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[2].fill_between(index, 0, NAOI_min_3std, color='lightgrey', alpha=0.1)
    
    axz[2].axhline(y=0,color='k')
    axz[2].axvline(x=2005,color='g')    
    #axz[1].plot(index,radius.values)
    axz[3].set_ylabel('Total Cyclone \nTime Steps [$\mathregular{month^{-1}}$]')
    axz[3].axhline(y=agetot_mean,color='k')
    axz[3].axvline(x=2005,color='g')   
    #axz[3].set_ylim([70,95])
    axz[3].fill_between(index, agetot_mean, agetot.values.squeeze(), where=agetot.values.squeeze() >= agetot_mean, facecolor='red', interpolate=True,zorder=2)
    axz[3].fill_between(index, agetot_mean, agetot.values.squeeze(), where=agetot.values.squeeze() <= agetot_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[3].fill_between(index, agetot_mean, agetot_plus_1std, color='lightgrey', alpha=0.7)
    axz[3].fill_between(index, agetot_mean, agetot_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[3].fill_between(index, agetot_mean, agetot_plus_3std, color='lightgrey', alpha=0.1)
    #axz[3].fill_between(index, agetot_mean, agetot_min_1std, color='lightgrey', alpha=0.7)
    axz[3].fill_between(index, agetot_mean, agetot_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[3].fill_between(index, agetot_mean, agetot_min_3std, color='lightgrey', alpha=0.1)
    #axz[2].plot(index,depth.values)
    axz[4].set_ylabel('Depth [dam]')
    axz[4].axhline(y=depth_mean,color='k')
    axz[4].axvline(x=2005,color='g')   
    axz[4].fill_between(index, depth_mean, depth.values.squeeze(), where=depth.values.squeeze() >= depth_mean, facecolor='red', interpolate=True,zorder=2)
    axz[4].fill_between(index, depth_mean, depth.values.squeeze(), where=depth.values.squeeze() <= depth_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[4].fill_between(index, depth_mean, depth_plus_1std, color='lightgrey', alpha=0.7)
    axz[4].fill_between(index, depth_mean, depth_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[4].fill_between(index, depth_mean, depth_plus_3std, color='lightgrey', alpha=0.1)
    #axz[4].fill_between(index, depth_mean, depth_min_1std, color='lightgrey', alpha=0.7)
    axz[4].fill_between(index, depth_mean, depth_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[4].fill_between(index, depth_mean, depth_min_3std, color='lightgrey', alpha=0.1)
    #axz[3].plot(index,slp.values)
    axz[5].set_ylabel('SLP [hPa]')
    axz[5].axhline(y=slp_mean,color='k')
    axz[5].axvline(x=2005,color='g')   
    axz[5].fill_between(index, slp_mean, slp.values.squeeze(), where=slp.values.squeeze() >= slp_mean, facecolor='red', interpolate=True,zorder=2)
    axz[5].fill_between(index, slp_mean, slp.values.squeeze(), where=slp.values.squeeze() <= slp_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[5].fill_between(index, slp_mean, slp_plus_1std, color='lightgrey', alpha=0.7)
    axz[5].fill_between(index, slp_mean, slp_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[5].fill_between(index, slp_mean, slp_plus_3std, color='lightgrey', alpha=0.1)
    #axz[5].fill_between(index, slp_mean, slp_min_1std, color='lightgrey', alpha=0.7)
    axz[5].fill_between(index, slp_mean, slp_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[5].fill_between(index, slp_mean, slp_min_3std, color='lightgrey', alpha=0.1)
    #axz[4].plot(index,prec.values)
    axz[6].set_ylabel('Precipitation [mm]')
    axz[6].axhline(y=prec_mean,color='k')
    axz[6].axvline(x=2005,color='g')   
    #axz[4].set_ylim([18,21.5])
    axz[6].fill_between(index, prec_mean, prec.values.squeeze(), where=prec.values.squeeze() >= prec_mean, facecolor='red', interpolate=True,zorder=2)
    axz[6].fill_between(index, prec_mean, prec.values.squeeze(), where=prec.values.squeeze() <= prec_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[6].fill_between(index, prec_mean, prec_plus_1std, color='lightgrey', alpha=0.7)
    axz[6].fill_between(index, prec_mean, prec_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[6].fill_between(index, prec_mean, prec_plus_3std, color='lightgrey', alpha=0.1)
    #axz[6].fill_between(index, prec_mean, prec_min_1std, color='lightgrey', alpha=0.7)
    axz[6].fill_between(index, prec_mean, prec_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[6].fill_between(index, prec_mean, prec_min_3std, color='lightgrey', alpha=0.1)
    #axz[5].plot(index,ws.values)
    axz[7].set_ylabel('Wind Speed [m/s]')
    axz[7].axhline(y=ws_mean,color='k')
    axz[7].axvline(x=2005,color='g')   
    axz[7].fill_between(index, ws_mean, ws.values.squeeze(), where=ws.values.squeeze() >= ws_mean, facecolor='red', interpolate=True,zorder=2)
    axz[7].fill_between(index, ws_mean, ws.values.squeeze(), where=ws.values.squeeze() <= ws_mean, facecolor='blue', interpolate=True,zorder=2)
    #axz[7].fill_between(index, ws_mean, ws_plus_1std, color='lightgrey', alpha=0.7)
    axz[7].fill_between(index, ws_mean, ws_plus_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[7].fill_between(index, ws_mean, ws_plus_3std, color='lightgrey', alpha=0.1)
    #axz[7].fill_between(index, ws_mean, ws_min_1std, color='lightgrey', alpha=0.7)
    axz[7].fill_between(index, ws_mean, ws_min_2std, color='grey', alpha=0.7,edgecolor='none')
    #axz[7].fill_between(index, ws_mean, ws_min_3std, color='lightgrey', alpha=0.1)
    
    
    for ax in axz.flat:
        ax.grid()
        ax.set_xlabel('Year')
        ax.label_outer()
        if period=='past':
            ax.set_xlim([-1400,1800])
        elif period=='fut':
            ax.set_xlim([1700,2100])
    
    fig.subplots_adjust(left=.05,right=.99,top=.95,bottom=.045,hspace=.075)
    fig.suptitle(track_loc)
    #fig.suptitle('Temperature Anomalies and 1% most extreme cyclone statistics over Europe in DJF (30 Year Average) (1400 BC - 1800 AD)')
    fig.savefig(path+'figs/Timeseries_future_99p_track_{}_extended_{}.png'.format(track_loc,period),dpi=300)
    #plt.close(fig)

    



    import matplotlib
    df_cov = pd.DataFrame({
                        "T850 Anomaly [K]":T850.values.squeeze(),
                        "RWP ENvelope [m/s]": RWP.values.squeeze(),
                        "NAOI [-]": NAOI.values.squeeze(),
                        "Total Cyclone \nTime Steps [$\mathregular{month^{-1}}$]" : agetot.values.squeeze(),
                        "Depth [dam]" : depth.values.squeeze(),
                        "SLP [hPa]" : slp.values.squeeze(),
                        "Precipitation [mm]" : prec.squeeze(),
                        "Wind Speed [m/s]" : ws.squeeze()})
    
    df_cov_past = df_cov[df_cov.index<=3352].dropna(how='all')
    
    colorz = [(0,'#006600'),(0.125,'#33cc33'),
            (0.25,'#ffff00'),(0.375,'#ff6600'),(0.5,'#ff0000'),(0.625,'#ff6600'),(0.75,'#ffff00'),
            (0.875,'#33cc33'),(1,'#006600')]   
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom', colorz, N=1024)
    
    def cov_plot(df):
        columns = df.columns
        fig,axz = plt.subplots(len(columns),len(columns),figsize=(16,9))
        for y in range(len(columns)):
            for x in range(len(columns)):
                ax=axz[x,y]
                if y==x:
                    df[columns[x]].plot.hist(ax=ax,bins=30)
                    ax.set_ylabel([])
                    # ax.text(0.5 , 0.5 , '1',
                    # horizontalalignment='center',
                    # verticalalignment='center',
                    # transform=ax.transAxes,
                    # fontsize=16)
                    # ax.set_facecolor("lightgrey")
                    # ax.set_yticklabels([])
                if y<x:
                    ax.scatter(df[columns[y]],df[columns[x]],
                            marker='.',s=0.2)
        
                if x<y:
                    corr = df[columns[[y,x]]].corr().values[1,0]
                    ax.text(0.5 , 0.5 , '{:.3f}'.format(corr),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes,
                    fontsize=16)
                    ax.set_facecolor(cmap((corr+1)/2))
                ax.set_xlabel(columns[y],fontsize=8)
                ax.set_ylabel(columns[x],fontsize=8)
                ax.label_outer()
        fig.tight_layout()
        fig.savefig(path+'figs/cov_plot_past_track_{}_extended_{}.png'.format(track_loc,period),dpi=300)
        return fig
            
    fig = cov_plot(df_cov_past)
    del(NAOI)
    #plt.close(fig)


# %%
