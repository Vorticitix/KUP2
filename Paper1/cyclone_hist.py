#%%

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import sys
path='/storage/climatestor/PleioCEP/doensen/data/'
#%%
regions = ['med','natl']
file_era = 'fort_36_total_{}_ERA5.txt'
file_cesm = 'fort_36_total_{}.txt'

df_era = pd.read_csv(path+file_era.format(regions[0])).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
df_cesm = pd.read_csv(path+file_cesm.format(regions[0])).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
df_cesm['year'] = df_cesm['year'] - 1502
df_cesm = df_cesm[(df_cesm['year']>=1981)&(df_cesm['year']<=2010)]

#%%
data = [df_cesm,df_era]
data_grouped = []
seasons = [(3,5),(6,8),(9,11),(12,2)]
fig,axz =plt.subplots(4,5,figsize=(10,8))
xlabels = ['Total cyclone\n lifetime [days]',
           'Maximum cyclone\n depth [gpm]',
           'Minimum sea level\n pressure [hPa]',
           'Maximum precipitation\n rate [mm/day]',
           'Maximum\n wind speed [m/s]']
ylabels = ['MAM','JJA','SON','DJF']
           
for i in range(len(seasons)):
    for df in data:
        df_gb_max = df[['zrad','zdep','precmean','preclmean','wsmean','iic']].groupby('iic').max()
        df_gb_min = df[['slp','iic','agetot']].groupby('iic').min()
        df_gb_yas = df[['month','iic']].groupby('iic').first()
        df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
        df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                                    &(df_gb_iic['zrad']!=0)
                                    &(df_gb_iic['zdep']!=0)).dropna(how='all')
        data_grouped.append(df_gb_iic)
    df_cesm_grouped = data_grouped[0]
    df_era_grouped = data_grouped[1]
    if seasons[i][1]<seasons[i][0]:
        df_era_season = df_era_grouped[(df_era_grouped['month']>=seasons[i][0])|(df_era_grouped['month']<=seasons[i][1])]
        df_cesm_season = df_cesm_grouped[(df_cesm_grouped['month']>=seasons[i][0])|(df_cesm_grouped['month']<=seasons[i][1])]
    else:
        df_era_season = df_era_grouped[(df_era_grouped['month']>=seasons[i][0])&(df_era_grouped['month']<=seasons[i][1])]
        df_cesm_season = df_cesm_grouped[(df_cesm_grouped['month']>=seasons[i][0])&(df_cesm_grouped['month']<=seasons[i][1])]  
    #Agetot
    ax=axz[i,0]
    agetot=[df_cesm_season['agetot'].values/4,df_era_season['agetot'].values/4]
    ax.hist(agetot,bins=np.linspace(0,12,13),label=['CESM','ERA5'],density=True)
    ax.set_xlabel(xlabels[0],fontsize=10)
    ax.set_ylabel(ylabels[i])
    ax.grid()
    #zdep
    ax=axz[i,1]
    zdep=[df_cesm_season['zdep'].values,df_era_season['zdep'].values]
    ax.hist(zdep,bins=np.linspace(0,420,13),label=['CESM','ERA5'],density=True)
    ax.grid()
    ax.set_xlabel(xlabels[1],fontsize=10)
    #slp
    ax=axz[i,2]
    slp=[df_cesm_season['slp'].values,df_era_season['slp'].values*100]
    ax.hist(slp,bins=np.linspace(970,1030,13),label=['CESM','ERA5'],density=True)
    ax.grid()
    ax.set_xlabel(xlabels[2],fontsize=10)
    #precmean
    ax=axz[i,3]
    prec=[df_cesm_season['precmean'].values,df_era_season['precmean'].values]
    ax.hist(prec,bins=np.linspace(0,48,13),label=['CESM','ERA5'],density=True)
    ax.grid()
    ax.set_xlabel(xlabels[3],fontsize=10)
    #wsmean
    ax=axz[i,4]
    ws=[df_cesm_season['wsmean'].values,df_era_season['wsmean'].values]
    ax.hist(ws,bins=np.linspace(0,32,13),label=['CESM','ERA5'],density=True)
    ax.grid()
    ax.set_xlabel(xlabels[4],fontsize=10)
    if i == 0:
        ax.legend(loc='upper center', bbox_to_anchor=(0.43, 1.25),
           ncol=2)
for ax in axz.flat:
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.yaxis.set_major_locator(plt.NullLocator())
    
    ax.label_outer()

fig.suptitle('Cyclone features CESM vs ERA5 for the Mediterannean (1981-2010)                              ')
# handles, labels = ax.get_legend_handles_labels()
# fig.legend(handles, labels, loc='upper right')
fig.subplots_adjust(left=.05,right=.95,bottom=.1,top=.95,wspace=.1,hspace=.1)
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/figs/Paper1/hist_era5.png',dpi=300)

    
#%%
        
