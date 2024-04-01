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
    full_time = np.around(np.arange(5,3515,1/12),3)
    PC1 = xr.open_dataset(path+'extracted/EOF_season/PC_mode_01.nc')
    time_arr = PC1.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    PC1 = PC1.assign_coords(time=time)
    PC1 = PC1.reindex(time=full_time)
    PC1 = PC1.pc.rolling(time=360,min_periods=340).mean(dim='time',skipna=True)
    PC1_std = float(PC1.sel(time=slice(5,3352)).std())
    PC1_plus_1std = PC1_std; PC1_min_1std = -PC1_std
    PC1_plus_2std = 2*PC1_std; PC1_min_2std = -2*PC1_std
    PC1_plus_3std = 3*PC1_std; T850_min_3std = -3*PC1_std
    
    PC2 = xr.open_dataset(path+'extracted/EOF_season/PC_mode_02.nc')
    time_arr = PC2.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    PC2 = PC2.assign_coords(time=time)
    PC2 = PC2.reindex(time=full_time)
    PC2 = PC2.pc.rolling(time=360,min_periods=340).mean(dim='time',skipna=True)
    PC2_std = float(PC2.sel(time=slice(5,3352)).std())
    PC2_plus_1std = PC2_std; PC2_min_1std = -PC2_std
    PC2_plus_2std = 2*PC2_std; PC2_min_2std = -2*PC2_std
    PC2_plus_3std = 3*PC2_std; T850_min_3std = -3*PC2_std
    
    PC3 = xr.open_dataset(path+'extracted/EOF_season/PC_mode_03.nc')
    time_arr = PC3.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    PC3 = PC3.assign_coords(time=time)
    PC3 = PC3.reindex(time=full_time)
    PC3 = PC3.pc.rolling(time=360,min_periods=340).mean(dim='time',skipna=True)
    PC3_std = float(PC3.sel(time=slice(5,3352)).std())
    PC3_plus_1std = PC3_std; PC3_min_1std = -PC3_std
    PC3_plus_2std = 2*PC3_std; PC3_min_2std = -2*PC3_std
    PC3_plus_3std = 3*PC3_std; T850_min_3std = -3*PC3_std
    
    PC4 = xr.open_dataset(path+'extracted/EOF_season/PC_mode_04.nc')
    time_arr = PC4.time.values
    time = np.around(np.array([x.year + x.month/12-1/12 for x in time_arr]),3)
    PC4 = PC4.assign_coords(time=time)
    PC4 = PC4.reindex(time=full_time)
    PC4 = PC4.pc.rolling(time=360,min_periods=340).mean(dim='time',skipna=True)
    PC4_std = float(PC4.sel(time=slice(5,3352)).std())
    PC4_plus_1std = PC4_std; PC4_min_1std = -PC4_std
    PC4_plus_2std = 2*PC4_std; PC4_min_2std = -2*PC4_std
    PC4_plus_3std = 3*PC4_std; T850_min_3std = -3*PC4_std
    
    print(full_time)
    
    seasons_dic = {1:'DJF',2:'DJF',3:'MAM',4:'MAM',5:'MAM',6:'JJA',
                7:'JJA',8:'JJA',9:'SON',10:'SON',11:'SON',12:'DJF'}
    names = ['Mediterranean [47 °N - 28 °N, 10 °W - 40 °E]',
            'Atlantic [30 °N - 30 °N, 65 °W - 10 °W]']
    path = "/storage/climatestor/PleioCEP/doensen/data/"
    fig,axz=plt.subplots(4,1,figsize=(15,5))
    plt.rcParams['axes.labelsize'] = 10
    index = PC1.time - 1502
    axz[0].set_ylabel('PC1\n(NAO-like)',fontsize=12)
    axz[0].fill_between(index, 0, PC1.values, where=PC1.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[0].fill_between(index, 0, PC1.values, where=PC1.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    axz[0].fill_between(index, 0, PC1_plus_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[0].fill_between(index, 0, PC1_min_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[0].axhline(y=0,color='k')

    axz[1].set_ylabel('PC2\n(EA-like)',fontsize=12)
    axz[1].fill_between(index, 0, PC2.values, where=PC2.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[1].fill_between(index, 0, PC2.values, where=PC2.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    axz[1].fill_between(index, 0, PC2_plus_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[1].fill_between(index, 0, PC2_min_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[1].axhline(y=0,color='k')
    
    axz[2].set_ylabel('PC3\n(EAWR-like)',fontsize=12)
    axz[2].fill_between(index, 0, PC3.values, where=PC3.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[2].fill_between(index, 0, PC3.values, where=PC3.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    axz[2].fill_between(index, 0, PC3_plus_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[2].fill_between(index, 0, PC3_min_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[2].axhline(y=0,color='k')
    
    axz[3].set_ylabel('PC4\n(SCAN-like)',fontsize=12)
    axz[3].fill_between(index, 0, PC4.values, where=PC4.values >= 0, facecolor='red', interpolate=True,zorder=2)
    axz[3].fill_between(index, 0, PC4.values, where=PC4.values <= 0, facecolor='blue', interpolate=True,zorder=2)
    axz[3].fill_between(index, 0, PC4_plus_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[3].fill_between(index, 0, PC4_min_2std, color='grey', alpha=0.5,edgecolor='none')
    axz[3].axhline(y=0,color='k')
    
    
    for ax in axz.flat:
        ax.grid()
        ax.set_xlabel('Year',fontsize=12)
        ax.label_outer()
        if period=='past':
            ax.set_xlim([-1500,1850])
        elif period=='fut':
            ax.set_xlim([1700,2100])
    
    fig.subplots_adjust(left=.07,right=.985,top=.92,bottom=.09,hspace=.075)
    if period=='past':
        fig.suptitle('30-year running average teleconnection modes',fontsize=16)
    else:
        fig.suptitle('Temperature Anomalies and Cyclone Features in the Mediterranean using a RCP8.5 scenario DJF',fontsize=16)
    #fig.suptitle('Temperature Anomalies and 1% most extreme cyclone statistics over Europe in DJF (30 Year Average) (1400 BC - 1800 AD)')
    fig.savefig(path+'figs/Toulouse/Timeseries_{}_{}_modes.png'.format(track_loc,period),dpi=500)
    #plt.close(fig)

    sys.exit()




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
    
    df_cov_past = df_cov[df_cov.index<=3514].dropna(how='all')
    
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
