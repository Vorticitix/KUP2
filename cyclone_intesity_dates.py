#%%
import xarray as xr
import wrf
from netCDF4 import Dataset 
import numpy as np
import pandas as pd
import sys
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")

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

for region in locations:
    for event in events:
        df_dates = pd.DataFrame(columns=['iic','year_start','month_start','day_start','hour_start',
                                    'year_end','month_end','day_end','hour_end'])
        df = pd.read_csv(path+'cyclone_intensity/{}_ext_{}.csv'.format(event,region),index_col=0)
        if event == 'precmax':
            df = df.iloc[:100]
        if event == 'wsmax':
            df_= pd.read_csv(path+'cyclone_intensity/{}_ext_{}.csv'.format('precmax',region),index_col=0)
            bool = ~np.isin(df['iic'].iloc[:100],df_['iic'].iloc[:100])
            df = df.iloc[:100].iloc[bool]
        if event == 'compoundmax':
            df_1 = pd.read_csv(path+'cyclone_intensity/{}_ext_{}.csv'.format('precmax',region),index_col=0)
            df_2 = pd.read_csv(path+'cyclone_intensity/{}_ext_{}.csv'.format('wsmax',region),index_col=0)
            bool = ~(np.isin(df['iic'].iloc[:100],df_1['iic'].iloc[:100]) | np.isin(df['iic'].iloc[:100],df_2['iic'].iloc[:100]))
            df = df.iloc[:100].iloc[bool]
        
        
        for rank,iic in enumerate(df['iic'].values):
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
            idxs= df_sel.index.where(total_bool==True)
            idx_min = idxs.min(); idx_max = idxs.max()
            
            df_sel['bool'].loc[idx_min:idx_max+1]=True
            date_bool = date[df_sel['bool']].astype(int).astype(str).str.zfill(8)
            idumi_bool = (idumi[df_sel['bool']].astype(int)//10000).astype(str).str.zfill(2)
            dummy_date = np.datetime64('1950-{:02d}-{:02d}T{:02d}:00:00'.format(int(date_bool.iloc[0][4:6]),int(date_bool.iloc[0][6:]),int(idumi_bool.iloc[0])))
            spinup_date = pd.Timestamp(dummy_date - np.timedelta64(6,'h'))
            if (spinup_date.month==12)&(spinup_date.day==31)&(spinup_date.hour==18):
                date_bool.loc[date_bool.index.min()-1]='{}{:02d}{:02d}'.format(int(date_bool.iloc[0][:4])-1,
                                                spinup_date.month,
                                                spinup_date.day)
            else:
                date_bool.loc[date_bool.index.min()-1]='{}{:02d}{:02d}'.format(date_bool.iloc[0][:4],
                                                spinup_date.month,
                                                spinup_date.day)
            idumi_bool.loc[idumi_bool.index.min()-1]='{:02d}'.format(spinup_date.hour)
            date_bool = date_bool.sort_index()
            idumi_bool = idumi_bool.sort_index()
            hour = idumi_bool.astype(int).values
            year = [int(x[:4]) for x in date_bool]
            month = [int(x[4:6]) for x in date_bool]
            day = [int(x[6:]) for x in date_bool]
            df_dates.loc[rank+1]=['{:08d}'.format(int(iic)),year[0],month[0],day[0],hour[0],year[-1],month[-1],day[-1],hour[-1]]
        df_dates.to_csv(path+'cyclone_intensity/dates/{}_ext_{}_dates_CSCS_test.csv'.format(event,region),
                        sep='\t',index=None,header=None)
# %%
