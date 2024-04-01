#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:00:02 2022

@author: doensen
"""

import pandas as pd
import sys
path = '/storage/climatestor/PleioCEP/doensen/data/'
#%%
regions = ['natl','satl','eur','med','wmed','cmed','emed','whem']

for i,region in enumerate(regions):

    file_BP = 'fort_36_total_{}.txt'.format(region)
    file_RCP85 = 'fort_36_total_{}_RCP85.txt'.format(region)

    df_BP = pd.read_csv(path+file_BP,header=0,index_col=0)
    df_BP = df_BP.where(df_BP['year']<=3514).dropna()
    df_RCP85 = pd.read_csv(path+file_RCP85,header=0,index_col=0)
    df_RCP85 = df_RCP85.where(df_RCP85['year']>3514).dropna()
    frames = [df_BP,df_RCP85]

    df_total = pd.concat(frames)
    df_total.to_csv(path+'fort_36_total_{}_full.txt'.format(region))
# %%
