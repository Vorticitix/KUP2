import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
from my_tools import *


#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone/'
file_PI = 'fort_34_total.txt'
file_RCP85 = 'RCP85/fort_34_RCP85_total.txt'
header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
             'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmax','precmean']
zeiger = 31000000

df_PI  = pd.read_csv(path+file_PI,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()
df_RCP85 = pd.read_csv(path+file_RCP85,delim_whitespace=True,names=header,skiprows=[0,1]).dropna()

df_RCP85['iic'] = df_RCP85['iic'] + zeiger
df_PI['date'] = [float(x.replace(x[:4],str(int(x[:4])-1502))) for x in df_PI['date'].astype(str)]
df_PI.loc[df_PI['date'] > 18600000, 'iic'] = df_PI['iic'] + 1000000
df = pd.concat([df_PI,df_RCP85])
df['lat']=ilat_to_lat(df['ilat'])
df['lon']=ilon_to_lon(df['ilon'])
df['lon']=((df['lon'] + 180) % 360) - 180
#df = df.where((df['lat']<=47)&(df['lat']>=28)&(df['lon']>=-10)&(df['lon']<=40)).dropna()
df = df.where(df['date']>1e7).dropna()
df['month'] = [int(x[4:6]) for x in df['date'].astype(str)]
df['year'] = [int(x[:4]) for x in df['date'].astype(str)]
df = df.where(df['year']>=1750).dropna().reset_index()
df_jja = df.where(df['month']==6).dropna()
#df_djf = df.where(np.isin(df['month'],[12,1,2])).dropna()
print(df_jja)
print(df_djf)
sys.exit()

#idxs = df[['iic','precmax']].groupby('iic')['precmax'].idxmax().values
df_med = df.where((df['lat']<=47)&(df['lat']>=28)&(df['lon']>=-10)&(df['lon']<=40)).dropna()
idxs_med = df_med[['iic','precmax']].groupby('iic')['precmax'].idxmax().values
df_atl = df.where((df['lat']<=60)&(df['lat']>=40)&(df['lon']>=-60)&(df['lon']<=-10)).dropna()
idxs_atl = df_atl[['iic','precmax']].groupby('iic')['precmax'].idxmax().values
# %%
df_atl_max = df_atl.loc[idxs_atl].dropna()
date = df_atl_max.groupby(['year','month']).min()['date']
date = [pd.to_datetime(x[:6]+'01',format='%Y%m%d') for x in date.astype(str)]

df_q50_5y_atl = df_atl_max.groupby(['year','month']).quantile(q=.5).rolling(window=12*3,center=True).mean()
#df_q90_5y = df_max.groupby(['year','month']).quantile(q=.9).rolling(window=12*3,center=True).mean()
df_q99_5y_atl = df_atl_max.groupby(['year','month']).quantile(q=.99).rolling(window=12*3,center=True).mean()

df_q50_30y_atl = df_atl_max.groupby(['year','month']).quantile(q=.5).rolling(window=30*12,center=True).mean()
#df_q90_30y = df_max.groupby(['year','month']).quantile(q=.9).rolling(window=30*12,center=True).mean()
df_q99_30y_atl = df_atl_max.groupby(['year','month']).quantile(q=.99).rolling(window=30*12,center=True).mean()

df_med_max = df_med.loc[idxs_med].dropna()
date_med = df_med_max.groupby(['year','month']).min()['date']
date_med = [pd.to_datetime(x[:6]+'01',format='%Y%m%d') for x in date_med.astype(str)]

df_q50_5y_med = df_med_max.groupby(['year','month']).quantile(q=.5).rolling(window=12*3,center=True).mean()
#df_q90_5y = df_max.groupby(['year','month']).quantile(q=.9).rolling(window=12*3,center=True).mean()
df_q99_5y_med = df_med_max.groupby(['year','month']).quantile(q=.99).rolling(window=12*3,center=True).mean()

df_q50_30y_med = df_med_max.groupby(['year','month']).quantile(q=.5).rolling(window=30*12,center=True).mean()
#df_q90_30y = df_max.groupby(['year','month']).quantile(q=.9).rolling(window=30*12,center=True).mean()
df_q99_30y_med = df_med_max.groupby(['year','month']).quantile(q=.99).rolling(window=30*12,center=True).mean()

#df_min = df.groupby('iic').min()


# %%

ticks = [pd.Timestamp('1750-01-01'),
          pd.Timestamp('1800-01-01'),
          pd.Timestamp('1850-01-01'),
          pd.Timestamp('1900-01-01'),
          pd.Timestamp('1950-01-01'),
          pd.Timestamp('2000-01-01'),
          pd.Timestamp('2050-01-01'),
          pd.Timestamp('2099-12-01')]
# ticks = [0,510,1020,1530,2040,2550,3060,3570,4079]
ticklabels = [1750,1800,1850,1900,1950,2000,2050,2100]

to_plot_30y = [df_q99_30y_atl,df_q50_30y_atl,df_q99_30y_med,df_q50_30y_med]
to_plot_5y = [df_q99_5y_atl,df_q50_5y_atl,df_q99_5y_med,df_q50_5y_med]
labels = ['p=99th','p=50th','p=99th','p=50th']
colors = ['r','orange','r','orange']
dates = [date,date,date_med,date_med]
lims = [[25,40],[9,12],[12,20],[4,8]]
yticks =  [[25,30,35,40],[9,10,11,12],[12,14,16,18,20],[4,5,6,7,8]]
titles = ['Northern Atlantic','Northern Atlantic','Mediterranean','Mediterranean']
fig,axz = plt.subplots(4,1,figsize=(4.5,7))

for i,ax in enumerate(axz.flat):

    ax.plot(dates[i],to_plot_30y[i]['precmax'],label=labels[i],
            color=colors[i],linewidth=3)
    ax.plot(dates[i],to_plot_5y[i]['precmax'],color=colors[i],linewidth=.75)
    ax.set_ylim(lims[i])
    if i<4:
        ax.legend(loc='upper left')
    else:
        ax.legend(loc='lower left')
    #ax.label_outer()
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.set_yticks(yticks[i])
    ax.grid()
    ax.margins(0.005)
    # ax2.set_ylabel('Total Daily Precipitation [mm]'")
    if i >=3:
        ax.set_xlabel('year')
    ax.set_title(titles[i])
fig.text(0.03, 0.5, 'Total Daily Precipitation [mm]', ha='center', va='center', rotation='vertical')
fig.subplots_adjust(left=0.125,bottom=0.075,top=0.96,right=0.95,
                    wspace=0.15,hspace=0.4)
fig.savefig('/storage/mirrored/homes/doensen/poster_figs/prec_timeseries.png',dpi=300)
    
    
    
