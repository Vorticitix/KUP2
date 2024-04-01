#%%
import xarray as xr
import wrf
from netCDF4 import Dataset 
import numpy as np
import pandas as pd
import sys
from matplotlib import pyplot as plt

#%%
file = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/geo_em.d02.nc'
ds = xr.open_dataset(file).squeeze()
lons_raw = ds['XLONG_M'].values
lats_raw = ds['XLAT_M'].values
max_lon_diff = np.abs(np.diff(lons_raw,axis=1)).max()
max_lat_diff = np.abs(np.diff(lats_raw,axis=0)).max()


locations = {
               'emed':[(17.6,40),(28,42)],
               'cmed':[(2.6,17.5),(28,47)],
               #'wmed':[(-10,2.5),(28,40)],
               #'med':[(-10,40),(28,47)]
              }

events = ['precmax','wsmax','compoundmax']
 # Works as expected
# %%
path = '/storage/climatestor/PleioCEP/doensen/data/'
file = 'cyclone_intensity/precmax_ext_cmed.csv'

df = pd.read_csv(path+file,index_col=0)
df_iic = pd.read_csv(path+'fort_36_total_whem.txt',index_col=0)
df_iic = df_iic[['iic','lat','lon','date','idumi','slp']]
df_iic['lon'] = (df_iic['lon'] + 180) % 360 - 180


# %%
df_dates = pd.DataFrame(columns=['iic','year_start','month_start','day_start','hour_start',
                                    'year_end','month_end','day_end','hour_end'])
for region in locations:
    for event in events:
        df = pd.read_csv(path+'cyclone_intensity/{}_ext_{}.csv'.format(event,region),index_col=0)
        for rank,iic in enumerate(df['iic'].values):
            if rank > 100:
                continue
            print(rank)
            idx = np.where(df_iic['iic']==iic)[0]
            df_sel = df_iic.iloc[idx]
            df_sel['slp'] = df_sel['slp'].where((df_sel['lon']>=locations[region][0][0])&
                                (df_sel['lon']<=locations[region][0][1])&
                                (df_sel['lat']>=locations[region][1][0])&
                                (df_sel['lat']<=locations[region][1][1]))
            central_idx = df_sel['slp'].reset_index().idxmin().values[1]
            min_bnd = -central_idx
            max_bnd = len(df_sel['slp'])-central_idx
            df_sel.index = np.arange(min_bnd,max_bnd)
            lat = df_sel['lat']
            lon = df_sel['lon']
            date = df_sel['date']
            idumi = df_sel['idumi']
            lat_bool = [np.abs(x-lats_raw).min() for x in lat]<max_lat_diff
            lon_bool = [np.abs(x-lons_raw).min() for x in lon]<max_lon_diff
            total_bool = (lat_bool)&(lon_bool)
            df_sel['bool']=False
            idxs_lower = df_sel.index.where(total_bool==False)[df_sel.index<0]
            idxs_higher = df_sel.index.where(total_bool==False)[df_sel.index>0]
            if len(idxs_lower)>0:
                idx_lower = idxs_lower.max()
                df_sel['bool'].loc[idx_lower+1:0]=True
            else:
                idx_lower = df_sel.index.min()
                df_sel['bool'].loc[idx_lower:0]=True
            if len(idxs_higher)>0:
                idx_higher = idxs_higher.min()
                df_sel['bool'].loc[0:idx_higher]=True
            else:
                idx_higher = df_sel.index.max()
                df_sel['bool'].loc[0:idx_higher+1]=True
                
            df_sel['bool'].loc[idx_lower+1:idx_higher]=True
            date_bool = date[(lat_bool)&(lon_bool)].astype(int).astype(str).str.zfill(8)
            hour = (idumi[(lat_bool)&(lon_bool)]//10000).astype(int).astype(str).str.zfill(2).values
            year = [x[:4] for x in date_bool]
            month = [x[4:6] for x in date_bool]
            day = [x[6:] for x in date_bool]
            
            df_dates.loc[rank+1]=[int(iic),year[0],month[0],day[0],hour[0],year[-1],month[-1],day[-1],hour[-1]]
        sys.exit()
# %%
