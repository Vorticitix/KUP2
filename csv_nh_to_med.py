#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:25:24 2022

@author: doensen
"""

import pandas as pd
from my_tools import *
import numpy as np
import sys

path = '/storage/climatestor/PleioCEP/doensen/data/'
file_30 = 'fort_30_total_RCP85.txt'
file_36 = 'fort_36_total_RCP85.txt'

file_30_med = 'fort_30_total_med_RCP85.txt'
file_36_med = 'fort_36_total_med_RCP85.txt'

# df_30 = pd.read_csv(path+file_30)
# df_36 = pd.read_csv(path+file_36)

header_36= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmean',
		 'precmax','preccmean','preclmean','wsmean','wsmax']
header_30 = ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad','zdep' ]
files = [file_30,file_36]
headers = [header_30,header_36]
files_to_write = [file_30_med,file_36_med]


for i,file in enumerate(files):
    if i == 0:
        continue
    df = pd.read_csv(path+file,delim_whitespace=True,names=headers[i],skiprows=[0,1]).dropna()
    df['date'] = df['date'].astype(np.int64)
    df['year'] = [int(x[:4]) for x in df['date'].astype(str).str.zfill(8)]
    df['month'] = [int(x[4:6]) for x in df['date'].astype(str).str.zfill(8)]
    df['dayofyear'] = [pd.to_datetime('2022-{}-{}'.format(x[4:6],x[6:8])).dayofyear\
                       for x in df['date'].astype(str).str.zfill(8)]
    #df = df.where(df['date']>1e7).dropna()
    df['lat']=ilat_to_lat(df['ilat'])
    df['lon']=ilon_to_lon(df['ilon'])
    df = df.where(((df['lat']>28)&(df['lat']<47))&((df['lon']>350)|(df['lon']<40)))\
        .dropna(how='any')
    #df = df.where((df['lon']>=240)|(df['lon']<=60)).dropna().reset_index()
    #df = df.where(df['date']>50000).dropna().reset_index()
    print(df)
    
    
    print('hello')
    #df = df.where((df['year']>=start_year)&(df['year']<=start_year + dy -1)).dropna()

    df.to_csv(path+files_to_write[i])
