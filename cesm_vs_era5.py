#%%
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
path = '/storage/climatestor/PleioCEP/doensen/data/'
Flag=True
#%%
area = 'eur'
file_cesm = 'fort_36_total_{}.txt'
file_era5 = 'fort_36_total_{}_ERA5.txt'
file_rcp85 = 'fort_36_total_{}_RCP85.txt'

df_dic = {'cesm_med':pd.read_csv(path+file_cesm.format('med')),
          'era5_med':pd.read_csv(path+file_era5.format('med')),
          'rcp85_med':pd.read_csv(path+file_rcp85.format('med')),
          'cesm_satl':pd.read_csv(path+file_cesm.format('satl')),
          'era5_satl':pd.read_csv(path+file_era5.format('satl')),
          'rcp85_satl':pd.read_csv(path+file_rcp85.format('satl')),
          'cesm_natl':pd.read_csv(path+file_cesm.format('natl')),
          'era5_natl':pd.read_csv(path+file_era5.format('natl')),
          'rcp85_natl':pd.read_csv(path+file_rcp85.format('natl')),
          'cesm_eur':pd.read_csv(path+file_cesm.format('eur')),
          'era5_eur':pd.read_csv(path+file_era5.format('eur')),
          'rcp85_eur':pd.read_csv(path+file_rcp85.format('eur')),}

df_gb_dic = {}
regions = ['med','satl','natl','eur']
#%%
for key in df_dic.keys():
    df = df_dic[key]
    if 'rcp85' in key:
        df['year']=(df['year'].astype(int)-1502)
        df = df.where((df['year']>=2070)&(df['year']<=2099))
    elif 'cesm' in key:
        df['year']=(df['year'].astype(int)-1502)
        df = df.where((df['year']>=1981)&(df['year']<=2010))
    df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
    df['yearmonth']=df['year'] + df['month']/12 -1/12 
    df_gb_max = df[['agetot','zrad','zdep','precmax','precmean','wsmean','wsmax','iic']].groupby('iic').max()
    df_gb_min = df[['cz','iic']].groupby('iic').min()
    df_gb_yas = df[['yearmonth','iic']].groupby('iic').first()
    df_gb = pd.concat([df_gb_max,df_gb_yas,df_gb_min],axis=1)
    df_gb = df_gb.where((df_gb['zrad']!=0)&(df_gb['zdep']!=0)).dropna()
    df_gb_dic[key] = df_gb
    




# %%
#Plot
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    precmean_cesm = df_gb_dic['cesm_{}'.format(region)]['precmean']
    precmean_era5 = df_gb_dic['era5_{}'.format(region)]['precmean']
    precmean_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['precmean']
    ax.hist(precmean_cesm,bins=16,range=(0,40),alpha=.5,label='CESM 1981-2010 N={}'.format(len(precmean_cesm)),density=Flag)
    #ax.hist(precmean_era5,bins=16,range=(0,40),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(precmean_era5)),density=Flag)
    ax.hist(precmean_rcp85,bins=16,range=(0,40),alpha=.5,label='CESM 2070-2099 N={}'.format(len(precmean_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Mean Cyclone Related precipitation')
#fig.savefig(path+'figs/comparison/precmean_era5_density_{}.png'.format(str(Flag)))
fig.savefig(path+'figs/comparison/precmean_rcp85_density_{}.png'.format(str(Flag)))


# %%
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['precmax']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['precmax']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['precmax']
    ax.hist(var_cesm,bins=16,range=(0,40),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    #ax.hist(var_era5,bins=16,range=(0,40),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    ax.hist(var_rcp85,bins=16,range=(0,40),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Max Cyclone Related precipitation')
#fig.savefig(path+'figs/comparison/precmax_era5_density_{}.png'.format(str(Flag)))
fig.savefig(path+'figs/comparison/precmax_rcp85_density_{}.png'.format(str(Flag)))


# %%

fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['wsmean']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['wsmean']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['wsmean']
    ax.hist(var_cesm,bins=16,range=(0,40),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    #ax.hist(var_era5,bins=16,range=(0,40),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    ax.hist(var_rcp85,bins=16,range=(0,40),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Mean Cyclone Related Wind Speed')
#fig.savefig(path+'figs/comparison/wsmean_era5_density_{}.png'.format(str(Flag)))
fig.savefig(path+'figs/comparison/wsmean_rcp85_density_{}.png'.format(str(Flag)))

# %%
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['wsmax']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['wsmax']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['wsmax']
    ax.hist(var_cesm,bins=16,range=(0,40),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    #ax.hist(var_era5,bins=16,range=(0,40),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    ax.hist(var_rcp85,bins=16,range=(0,40),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Max Cyclone Related Wind Speed')
#fig.savefig(path+'figs/comparison/wsmax_era5_density_{}.png'.format(str(Flag)))
fig.savefig(path+'figs/comparison/wsmax_rcp85_density_{}.png'.format(str(Flag)))

# %%
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['cz']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['cz']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['cz']
    ax.hist(var_cesm,bins=40,range=(-500,300),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    ax.hist(var_era5,bins=40,range=(-500,300),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    #ax.hist(var_rcp85,bins=40,range=(-500,300),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('1000 hPa Geopotenital Height in Cyclone Center')
fig.savefig(path+'figs/comparison/cz_era5_density_{}.png'.format(str(Flag)))
#fig.savefig(path+'figs/comparison/cz_rcp85_density_{}.png'.format(str(Flag)))

# %%

fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['zdep']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['zdep']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['zdep']
    ax.hist(var_cesm,bins=25,range=(0,500),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    ax.hist(var_era5,bins=25,range=(0,500),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    #ax.hist(var_rcp85,bins=25,range=(0,500),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Cyclone Depth')
fig.savefig(path+'figs/comparison/zdep_era5_density_{}.png'.format(str(Flag)))
#fig.savefig(path+'figs/comparison/zdep_rcp85_density_{}.png'.format(str(Flag)))
# %%
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['agetot']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['agetot']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['agetot']
    ax.hist(var_cesm,bins=16,range=(0,40),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    ax.hist(var_era5,bins=16,range=(0,40),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    #ax.hist(var_rcp85,bins=16,range=(0,40),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Cyclone Lifetime')
fig.savefig(path+'figs/comparison/agetot_era5_density_{}.png'.format(str(Flag)))
#fig.savefig(path+'figs/comparison/agetot_rcp85_density_{}.png'.format(str(Flag)))
# %%
fig,axz = plt.subplots(2,2,figsize=(16,9))
for i,region in enumerate(regions):
    ax = axz.flat[i]
    var_cesm = df_gb_dic['cesm_{}'.format(region)]['zrad']
    var_era5 = df_gb_dic['era5_{}'.format(region)]['zrad']
    var_rcp85 = df_gb_dic['rcp85_{}'.format(region)]['zrad']
    ax.hist(var_cesm,bins=25,range=(0,0.5),alpha=.5,label='CESM 1981-2010 N={}'.format(len(var_cesm)),density=Flag)
    ax.hist(var_era5,bins=25,range=(0,0.5),alpha=.5,label='ERA5 1981-2010 N={}'.format(len(var_era5)),density=Flag)
    #ax.hist(var_rcp85,bins=25,range=(0,0.5),alpha=.5,label='CESM 2070-2099 N={}'.format(len(var_rcp85)),density=Flag)
    ax.legend()
    ax.set_title(region.upper())
fig.suptitle('Cyclone Radius')
fig.savefig(path+'figs/comparison/zrad_era5_density_{}.png'.format(str(Flag)))
#fig.savefig(path+'figs/comparison/zrad_rcp85_density_{}.png'.format(str(Flag)))
# %%
