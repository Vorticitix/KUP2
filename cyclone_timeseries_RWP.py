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
track_loc = 'natl'


path = "/storage/climatestor/PleioCEP/doensen/data/"
# file_past = 'fort_36_total_med.txt'
file_med = 'fort_36_total_{}_full.txt'.format(track_loc)


df=pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# med = pd.read_csv(path+file_med).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# rcp = pd.read_csv(path+file_rcp).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
# rcp = rcp.where(rcp['year']>=3515).dropna()
# df = pd.concat([med,rcp])
df = df.where(df['year']>=100).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
#df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
df['yearmonth']=df['year'] + df['month']/12 -1/12 
df['cen'] = (df['year']-1502)//100



#df_groupby= df.groupby('yearseason')
df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
df_gb_min = df[['slp','iic','agetot']].groupby('iic').min()
df_gb_yas = df[['yearmonth','iic','cen']].groupby('iic').first()
df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                            &(df_gb_iic['zrad']!=0)
                            &(df_gb_iic['zdep']!=0)).dropna(how='all')

#%%
period='fut'
agetot = df_gb_iic[['agetot','yearmonth']].groupby('yearmonth').sum().rolling(90).mean()
agetot_mean = agetot[agetot.index<=3352].dropna().mean().values[0]

radius = df_gb_iic[['zrad','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
radius_mean = radius[radius.index<=3352].dropna().mean().values[0]
#radius_cen = df_gb_iic[['zrad','cen']].groupby('cen').median()
# fig,ax = plt.subplots()
# ax.plot(radius.index-1502,radius.values)
#ax.plot(radius_cen.index*100,radius_cen.values)
# ax.grid()
#ax.set_xticks(radius.index)
#ax.set_ylim([0.15,0.25])

slp = df_gb_iic[['slp','yearmonth']].groupby('yearmonth').quantile(1-.99).rolling(90).mean()
slp_mean = slp[slp.index<=3352].dropna().mean().values[0]
# fig,ax = plt.subplots()
# ax.plot(slp.index-1502,slp.values)
# ax.grid()

depth = df_gb_iic[['zdep','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
depth_mean = depth[depth.index<=3352].dropna().mean().values[0]
# fig,ax = plt.subplots()
# ax.plot(depth.index-1502,depth.values)

prec = df_gb_iic[['precmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
prec_mean = prec[prec.index<=3352].dropna().mean().values[0]
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
ws = df_gb_iic[['wsmean','yearmonth']].groupby('yearmonth').quantile(.99).rolling(90).mean()
ws_mean = ws[ws.index<=3352].dropna().mean().values[0]
#ws = ws[ws.index<200].dropna()
# fig,ax = plt.subplots()
# ax.plot(ws.index,ws.values)

regions = ['natl','satl','med','eur']
for rwp_loc in regions:
    full_time = radius.index
    RWP = xr.open_dataset(path+'extracted/RWP/RWP_all_monmean_{}.nc'.format(rwp_loc)).squeeze()
    #RWP = RWP.rename({'__xarray_dataarray_variable__':"RWP"})
    # T850 = T850.assign_coords(lon=(((T850.lon.values + 180) % 360) - 180)).sortby('lon')
    time_arr = RWP.time.values
    time = np.array([x.year + x.month/12-1/12 for x in time_arr])
    u, c = np.unique(time, return_counts=True)
    sim  = u[c == 1]
    idx = np.where(np.isin(time,sim)==True)[0]
    RWP = RWP.assign_coords(time=time)
    RWP = RWP.isel(time=idx)
    RWP = RWP.reindex(time=full_time)
    RWP = RWP.V.rolling(time=90,min_periods=85).mean(dim='time',skipna=True)
    RWP_mean = float(RWP.sel(time=slice(0,3352)).mean())
    
    seasons_dic = {1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',6:'JJA',
                   7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON',12:'DJF'}
    names = ['Mediterranean [47 °N - 28 °N, 10 °W - 40 °E]',
             'Atlantic [70 °N - 30 °N, 65 °W - 10 °W]']
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
    fig,axz=plt.subplots(6,1,figsize=(17,10))
    index = radius.index - 1502
    full_time = radius.index
    #T850 = xr.open_dataset(path+'extracted/T850/T850_0005_3352_anom_med_mean.nc').squeeze()
    # T850 = T850.assign_coords(lon=(((T850.lon.values + 180) % 360) - 180)).sortby('lon')
    #time_arr = T850.time.values
    #time = np.array([x.year + x.month/12-1/12 for x in time_arr])
    #axz[0].plot(index,T850.values)
    axz[0].set_ylabel('RWP Envelope [m/s]')
    #axz[0].set_ylim([-.5,.5])
    axz[0].fill_between(index, RWP_mean, RWP.values, where=RWP.values >= RWP_mean, facecolor='red', interpolate=True)
    axz[0].fill_between(index, RWP_mean, RWP.values, where=RWP.values <= RWP_mean, facecolor='blue', interpolate=True)
    axz[0].axhline(y=RWP_mean,color='k')
    axz[0].axvline(x=2005,color='g')        
    #axz[1].plot(index,radius.values)
    axz[1].set_ylabel('Total Cyclone \nTime Steps [$\mathregular{month^{-1}}$]')
    axz[1].axhline(y=agetot_mean,color='k')
    axz[1].axvline(x=2005,color='g')   
    #axz[1].set_ylim([70,95])
    axz[1].fill_between(index, agetot_mean, agetot.values.squeeze(), where=agetot.values.squeeze() >= agetot_mean, facecolor='red', interpolate=True)
    axz[1].fill_between(index, agetot_mean, agetot.values.squeeze(), where=agetot.values.squeeze() <= agetot_mean, facecolor='blue', interpolate=True)
    #axz[2].plot(index,depth.values)
    axz[2].set_ylabel('Depth [dam]')
    axz[2].axhline(y=depth_mean,color='k')
    axz[2].axvline(x=2005,color='g')   
    axz[2].fill_between(index, depth_mean, depth.values.squeeze(), where=depth.values.squeeze() >= depth_mean, facecolor='red', interpolate=True)
    axz[2].fill_between(index, depth_mean, depth.values.squeeze(), where=depth.values.squeeze() <= depth_mean, facecolor='blue', interpolate=True)
    #axz[3].plot(index,slp.values)
    axz[3].set_ylabel('SLP [hPa]')
    axz[3].axhline(y=slp_mean,color='k')
    axz[3].axvline(x=2005,color='g')   
    axz[3].fill_between(index, slp_mean, slp.values.squeeze(), where=slp.values.squeeze() >= slp_mean, facecolor='red', interpolate=True)
    axz[3].fill_between(index, slp_mean, slp.values.squeeze(), where=slp.values.squeeze() <= slp_mean, facecolor='blue', interpolate=True)
    #axz[4].plot(index,prec.values)
    axz[4].set_ylabel('Precipitation [mm]')
    axz[4].axhline(y=prec_mean,color='k')
    axz[4].axvline(x=2005,color='g')   
    #axz[4].set_ylim([18,21.5])
    axz[4].fill_between(index, prec_mean, prec.values.squeeze(), where=prec.values.squeeze() >= prec_mean, facecolor='red', interpolate=True)
    axz[4].fill_between(index, prec_mean, prec.values.squeeze(), where=prec.values.squeeze() <= prec_mean, facecolor='blue', interpolate=True)
    #axz[5].plot(index,ws.values)
    axz[5].set_ylabel('Wind Speed [m/s]')
    axz[5].axhline(y=ws_mean,color='k')
    axz[5].axvline(x=2005,color='g')   
    axz[5].fill_between(index, ws_mean, ws.values.squeeze(), where=ws.values.squeeze() >= ws_mean, facecolor='red', interpolate=True)
    axz[5].fill_between(index, ws_mean, ws.values.squeeze(), where=ws.values.squeeze() <= ws_mean, facecolor='blue', interpolate=True)
    
    
    for ax in axz.flat:
        ax.grid()
        ax.set_xlabel('Year')
        ax.label_outer()
        if period=='past':
            ax.set_xlim([-1400,1800])
        elif period=='fut':
            ax.set_xlim([1700,2100])
    
    fig.subplots_adjust(left=.05,right=.99,top=.95,bottom=.045,hspace=.075)
    fig.suptitle(rwp_loc)
    #fig.suptitle('Temperature Anomalies and 1% most extreme cyclone statistics over Europe in DJF (30 Year Average) (1400 BC - 1800 AD)')
    fig.savefig(path+'figs/Timeseries_future_99p_track_{}_rwp_{}_{}.png'.format(track_loc,rwp_loc,period),dpi=300)
    #plt.close(fig)

    



    import matplotlib
    df_cov = pd.DataFrame({"RWP ENvelope [m/s]": RWP.values.squeeze(),
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
                ax.set_xlabel(columns[y])
                ax.set_ylabel(columns[x])
                ax.label_outer()
        fig.tight_layout()
        fig.savefig(path+'figs/cov_plot_past_track_{}_rwp_{}_{}.png'.format(track_loc,rwp_loc,period),dpi=300)
        return fig
            
    fig = cov_plot(df_cov_past)
    plt.close(fig)


# %%
