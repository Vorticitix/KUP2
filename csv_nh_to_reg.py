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
#file_30 = 'fort_30_total.txt'
file_36 = 'fort_36_total_ERA5.txt'

#file_30_output = ['fort_30_total_med.txt','fort_30_total_satl.txt','fort_30_total_natl.txt','fort_30_total_eur.txt']
file_36_output = [
                   'fort_36_total_med_ERA5.txt',
                  'fort_36_total_satl_ERA5.txt',
                  'fort_36_total_natl_ERA5.txt',
                   'fort_36_total_eur_ERA5.txt',
                  'fort_36_total_whem_ERA5.txt',
                  'fort_36_total_wmed_ERA5.txt',
                  'fort_36_total_cmed_ERA5.txt',
                  'fort_36_total_emed_ERA5.txt']

coords = [
    [350,40,28,47],
    [290,350,28,47],
    [290,350,47,70],
    [350,40,47,70],
    [260,60,25,85],
    [350,2.5,30,47],
    [2.6,17.5,30,47],
    [17.6,40,28,42],
    ]

# df_30 = pd.read_csv(path+file_30)
# df_36 = pd.read_csv(path+file_36)

header_36= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep','slp', 'precmean',
		 'precmax','qmean','preclmean','wsmean','wsmax']
header_30 = ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad','zdep' ]
files = [file_36]
headers = [header_36]


for i,file in enumerate(files):
    df = pd.read_csv(path+file,delim_whitespace=True,names=headers[i],skiprows=[0,1]).dropna()
    df['date'] = df['date'].astype(np.int64)
    df['year'] = [int(x[:4]) for x in df['date'].astype(str).str.zfill(8)]
    df['month'] = [int(x[4:6]) for x in df['date'].astype(str).str.zfill(8)]
    df['day'] = [int(x[6:8]) for x in df['date'].astype(str).str.zfill(8)]
    df = df.where(~((df['day']==29)&(df['month']==2))).dropna(how='all')
    print(df['date'])
    df['dayofyear'] = [pd.to_datetime('2022-{}-{}'.format(int(x[4:6]),int(x[6:8]))).dayofyear\
                       for x in df['date'].astype(str).str.zfill(8)]
    #df = df.where(df['date']>1e7).dropna()
    print(df)
    df['lat']=ilat_to_lat(df['ilat'])
    df['lon']=ilon_to_lon(df['ilon'])
    print('yo')
    for j,output in enumerate(file_36_output):
        if coords[j][0]>coords[j][1]:
            df_output = df.where(((df['lat']>=coords[j][2])&(df['lat']<=coords[j][3]))
                                &((df['lon']>=coords[j][0])|(df['lon']<=coords[j][1])))\
                .dropna(how='any')
        else:
            df_output = df.where(((df['lat']>=coords[j][2])&(df['lat']<=coords[j][3]))
                                &((df['lon']>=coords[j][0])&(df['lon']<=coords[j][1])))\
                .dropna(how='any')
        #df = df.where((df['lon']>=240)|(df['lon']<=60)).dropna().reset_index()
        #df = df.where(df['date']>50000).dropna().reset_index()
        print(df_output)
        
        
                #df = df.where((df['year']>=start_year)&(df['year']<=start_year + dy -1)).dropna()
    
        df_output.to_csv(path+output)
