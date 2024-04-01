#%%
import xarray as xr
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
from scipy import stats
#%%

#Load data
path='/storage/climatestor/PleioCEP/doensen/data/'
file='fort_36_total_med_full.txt'

df=pd.read_csv(path+file).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
df = df.where(df['year']>=5).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
#df['yearseason']=(df['year'].astype(int)-1502).astype(str)+'_'+df['month'].map(seasons_dic)
df['yearmonth']=np.around(df['year'].values + df['month'].values/12 -1/12,3)
#df['cen'] = (df['year']-1502)//100
#df_groupby= df.groupby('yearseason')
df_gb_max = df[['zrad','zdep','precmean','preccmean','preclmean','wsmean','iic']].groupby('iic').max()
df_gb_min = df[['year','slp','iic','agetot']].groupby('iic').min()
df_gb_yas = df[['yearmonth','iic']].groupby('iic').first()
df_gb_iic = pd.concat([df_gb_max,df_gb_min,df_gb_yas],axis=1)
df_gb_iic = df_gb_iic.where((df_gb_iic['slp']!=0)
                            &(df_gb_iic['zrad']!=0)
                            &(df_gb_iic['zdep']!=0)).dropna(how='all')

#%%
#Define compound events. If cyclone passes specific percentile for both wind speed and precipitation, it is a compound event.

qs = [.7,.8,.9,.95]
bars = []; winds = []; precs = []; compounds = []; prec_flags = []; ws_flags = []


for q in qs:
    prec_flag = (df_gb_iic['precmean']>df_gb_iic['precmean'].quantile(q))
    ws_flag = (df_gb_iic['wsmean']>df_gb_iic['wsmean'].quantile(q))
    flag = prec_flag & ws_flag

    compound = df_gb_iic[flag]
    prec_df = df_gb_iic[prec_flag]
    ws_df = df_gb_iic[ws_flag]
    compound['dec'] = (compound['year']-1502)//10 
    prec_df['dec'] = (prec_df['year']-1502)//10 
    ws_df['dec'] = (ws_df['year']-1502)//10 
    prec_df['prec_flag']=prec_flag
    ws_df['ws_flag']=ws_flag
    bar = compound.groupby('dec').count().reindex(np.arange(-150,210),fill_value=0)['precmean']
    prec =  compound.groupby('dec').mean().reindex(np.arange(-150,210),fill_value=0)['precmean']
    wind =  compound.groupby('dec').mean().reindex(np.arange(-150,210),fill_value=0)['wsmean']
    prec_flag_sum = prec_df.groupby('dec').sum().reindex(np.arange(-150,210),fill_value=0)['prec_flag']
    ws_flag_sum = ws_df.groupby('dec').sum().reindex(np.arange(-150,210),fill_value=0)['ws_flag']
    
    bars.append(bar)
    winds.append(wind)
    precs.append(prec)
    compounds.append(compound)
    prec_flags.append(prec_flag_sum)
    ws_flags.append(ws_flag_sum)
    


#%%
"""
Compound events vs Singular extreme events
"""
#T
fig,ax = plt.subplots(figsize=(12,8))
sums = [x.sum() for x in bars]
x = np.array([1,2,3,4])
bar_width = 0.35
ax.bar(x-bar_width/2,sums,color='blue',width=bar_width)
ax.set_ylabel('Number of Compound Events',color='blue')
ax2=ax.twinx()
ax2.bar(x+bar_width/2,[x.sum() for x in prec_flags],color='red',width=bar_width)
ax2.set_ylabel('Number of cyclones 5% most extreme precipitation',color='red')
ax.set_xticks(x)
ax.set_xticklabels(('.7','.8','.9','.95'))
fig.suptitle('Number of compound events according to percentile threshold')


"""
Evolution of compound events in contrast to singular events.
"""
#%%
start_dec=180; end_dec=210
fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    n=bars[i][(bars[i].index.values>=start_dec)&(bars[i].index.values<end_dec)]
    ax.bar(np.arange(start_dec,end_dec),n)
    xticks = n.index
    xtick_labels = ['{}0-{}9'.format(x,x) for x in n.index]
    ax.set_xticks(xticks[0::2],fontsize=8)
    ax.set_xticklabels(xtick_labels[0::2],rotation=45,fontsize=8,ha='right')
    #ax.set_xlabel('Decades')
    ax.set_ylabel('Number of Compound Events')
    ax.set_title('q = {}'.format(qs[i]))
    ax.margins(x=0)
fig.subplots_adjust(bottom=0.1, right=0.9, top=0.95,left=0.1,hspace=.3)
#%%
start_dec=180; end_dec=210
fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    n=prec_flags[i][(prec_flags[i].index.values>=start_dec)&(prec_flags[i].index.values<end_dec)]
    ax.bar(np.arange(start_dec,end_dec),n)
    xticks = n.index
    xtick_labels = ['{}0-{}9'.format(x,x) for x in n.index]
    ax.set_xticks(xticks[0::2],fontsize=8)
    ax.set_xticklabels(xtick_labels[0::2],rotation=45,fontsize=8,ha='right')
    #ax.set_xlabel('Decades')
    ax.set_ylabel('Number of Extreme Precipitation Event')
    ax.set_title('q = {}'.format(qs[i]))
    ax.margins(x=0)
fig.subplots_adjust(bottom=0.1, right=0.9, top=0.95,left=0.1,hspace=.3)
#%%
fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    n=ws_flags[i][(ws_flags[i].index.values>=start_dec)&(ws_flags[i].index.values<end_dec)]
    ax.bar(np.arange(start_dec,end_dec),n)
    xticks = n.index
    xtick_labels = ['{}0-{}9'.format(x,x) for x in n.index]
    ax.set_xticks(xticks[0::2],fontsize=8)
    ax.set_xticklabels(xtick_labels[0::2],rotation=45,fontsize=8,ha='right')
    #ax.set_xlabel('Decades')
    ax.set_ylabel('Number of Extreme Wind Speed Event')
    ax.set_title('q = {}'.format(qs[i]))
    ax.margins(x=0)
fig.subplots_adjust(bottom=0.1, right=0.9, top=0.95,left=0.1,hspace=.3)

#%%
#Combined number precipitation, wind speed and compound
start_dec=180; end_dec=210

i=1
data_to_plot = [
    prec_flags[i][(prec_flags[i].index.values>=start_dec)&(prec_flags[i].index.values<end_dec)],
    ws_flags[i][(ws_flags[i].index.values>=start_dec)&(ws_flags[i].index.values<end_dec)],
    bars[i][(bars[i].index.values>=start_dec)&(bars[i].index.values<end_dec)],
    ]
ylabels=[
    'Number of extreme \nprecipitation events',
    'Number of extreme \nwind speed events',
    'Number of extreme \ncompound events']
titles = [
    'Precipitation Events only',
    'Wind Speed Events only',
    'Compound Events',
]
fig,axs=plt.subplots(1,3,figsize=(15,7))
for j,ax in enumerate(axs.flat):
    n=data_to_plot[j]
    ax.bar(np.arange(start_dec,end_dec),n)
    xticks = n.index
    xtick_labels = ['{}0-{}9'.format(x,x) for x in n.index]
    ax.set_xticks(xticks[0::2],fontsize=8)
    ax.set_xticklabels(xtick_labels[0::2],rotation=45,fontsize=10,ha='right')
    #ax.set_xlabel('Decades')
    ax.set_ylabel(ylabels[j],fontsize=12)
    ax.set_title(titles[j],fontsize=12)
    ax.margins(x=0)
    ax.set_ylim([0,70])
    #ax.label_outer()
fig.subplots_adjust(bottom=0.12, right=0.95, top=0.9,left=0.05,hspace=.3)
fig.suptitle('Precipitation, wind speed and compound events passing the 80th percentile threshold (1800 - 2100) DJF',fontsize=16)
fig.savefig('/storage/climatestor/PleioCEP/doensen/data/figs/Toulouse/compound.png',dpi=300)
#%%

fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    ax.bar(np.arange(350,360),winds[i])
    ax.set_xlabel('Decades')
    ax.set_ylabel('Mean compound event wind speed')
    ax.set_title('q = {}'.format(qs[i]))

#%%

fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    ax.bar(np.arange(350,360),precs[i])
    ax.set_xlabel('Decades')
    ax.set_ylabel('Mean compound event precipitation')
    ax.set_title('q = {}'.format(qs[i]))
    
#%%  
"""
Regression
"""


fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    compound_fut = compounds[i][compounds[i]['year']>=3502]
    ax.scatter(compound_fut['yearmonth']-1502,compound_fut['precmean'])
    slope, intercept, _, p_value, _ = stats.linregress(compound_fut['yearmonth']-1502,compound_fut['precmean'])
    trend_line = intercept + slope * (compound_fut['yearmonth']-1502)
    ax.plot(compound_fut['yearmonth']-1502, trend_line, color='red')
    ax.set_xlabel('Year')
    ax.set_ylabel('Compound event precipitation')
    ax.set_title('q = {} --- p = {}'.format(qs[i],p_value))
    
#%%
fig,axs = plt.subplots(2,2,figsize=(16,9))
for i,ax in enumerate(axs.flat):
    compound_fut = compounds[i][compounds[i]['year']>=3502]
    ax.scatter(compound_fut['yearmonth']-1502,compound_fut['wsmean'])
    slope, intercept, _, p_value, _ = stats.linregress(compound_fut['yearmonth']-1502,compound_fut['wsmean'])
    trend_line = intercept + slope * (compound_fut['yearmonth']-1502)
    ax.plot(compound_fut['yearmonth']-1502, trend_line, color='red')
    ax.set_xlabel('Year')
    ax.set_ylabel('Compound event wind speed')
    ax.set_title('q = {} --- p = {}'.format(qs[i],p_value))
