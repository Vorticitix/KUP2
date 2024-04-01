#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
import pandas as pd
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
from my_tools import *
import cartopy.crs as ccrs
import xarray as xr
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, Normalize
# %%

# colorz = [(0,'#006600'),(0.125*2,'#33cc33'),
#           (0.25*2,'#ffff00'),(0.375*2,'#ff6600'),(0.5*2,'#ff0000')]
# cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom', colorz, N=31)

# def cmap_value(df,thres_up,thres_down):
#     y = np.zeros(df.values.shape)
#     for i,x in enumerate(df.values):
#         if x>thres_up:
#             y[i]=1
#         elif x<thres_down:
#             y[i]=0
#         else:
#             y[i] = (x-thres_down)/(thres_up - thres_down)
#     return y
# %%

path = "/storage/climatestor/PleioCEP/doensen/data/"
file = "fort_36_total_med.txt"

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep']#,'slp', 'precmax','precmean']

df = pd.read_csv(path+file).dropna(how='all')
df['lat']=ilat_to_lat(df['ilat'])
df['lon']=ilon_to_lon(df['ilon'])
df['year']=(df['year'].astype(int))
#df = df.where((df['year']>=2009)&(df['year']<=2009)).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
#df = df.where((df['agetot']>=32)&(df['lat']>=0)).dropna()
# thres_up = df['gz'].quantile(q=0.975)
# thres_down = df['gz'].quantile(q=.5)

#df = df.where((df['agetot']>=20)).dropna()
#df['month'] = [int(x[4:6]) for x in df['date'].astype(str)]
#df_red = df[['month','lat','lon','zdep','iic','gz']]
#cmap_value = [i-thres_down/(thres_up-thres_down) if i<=thres_ else i=1 for i in df['gz'].values]

# %%
years = np.linspace(5,245,41)[:-1]
cmap = matplotlib.cm.get_cmap('gist_rainbow')
for year in years:
    years_ = np.arange(year,year+6)
    fig,axz = plt.subplots(2,3,figsize=(17,8),subplot_kw={'projection':ccrs.PlateCarree()})
    for j,year_ in enumerate(years_):
        df_year = df.where((df['year']>=year_)&(df['year']<=year_)).dropna(how='all')
        ax = axz.flat[j]
        count = 0
        for i,iic in enumerate(df_year['iic'].unique()):
            df_sub = df_year.where(df_year['iic']==iic).dropna(how='all')
            boolean = np.argmin(df_sub['slp'].values)-4<0
            if boolean:
                continue
            count += 1
            lonz = df_sub['lon'].values; latz = df_sub['lat'].values
            lonz = [(x + 180) % 360 - 180  for x in lonz]
            # if ((lonz<300)&(lonz>60)).any():
            #     continue
            
            points = np.array([lonz, latz]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            #ax.gridlines()
            cval = i/len(df_year['iic'].unique())
            lc = LineCollection(segments,linewidths=2,color=cmap(cval))
            #lc.set_array(np.arange(len(df['iic'].unique())))
            #c_value = cmap_value(df_sub['gz'], thres_up,thres_down)
            # lc.set_array(c_value)
            # lc.set_linewidth(1)
            # lc.set_transform(ccrs.Geodetic())
            #ax.plot(lonz,latz,transform=ccrs.Geodetic(),linewidth=.5,color=cmap(c_value))
            line = ax.add_collection(lc)
        ax.coastlines()
        ax.set_global()
        ax.set_extent([-20, 50, 25, 55])
        ax.set_title('Year: {} Count: {}'.format(year_.astype(int),count))
            #ax.set_extent([-180,180,0,90])
    fig.tight_layout()
    fig.savefig('/storage/climatestor/PleioCEP/doensen/data/tracking/figs/tracks_med_{:04d}.png'.format(years_[0].astype(int)))
    plt.close(fig)
            #ax.set_title('Month = {}'.format(month))
        # fig.subplots_adjust(left=0.025,bottom=0.015,top=0.965,right=0.975,
        #                     wspace=0.02,hspace=0.11)
            #fig.savefig("/storage/climatestor/PleioCEP/doensen/data/tracking/figs/{}.png".format(iic),dpi=300)

# %%
