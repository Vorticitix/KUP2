#%%
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy
import seaborn
import cartopy.crs as ccrs

#%%
#Load cyclone script data
path = '/storage/climatestor/PleioCEP/doensen/data/'
file = 'fort_36_total_med.txt'

#Read data using pandas. Only use data from DJF and remove invalid data
df = pd.read_csv(path+file).drop(['Unnamed: 0'],axis=1).reset_index(drop=True)
df = df.where((df['year']>=5)&(df['year']<=3352)).dropna(how='all')
df['lon'] = (df['lon'] + 180) % 360 - 180
df = df.where((df['lat']<42)|(df['lon']<23))
#df = df.where(df['day']==1).dropna(how='all')
df = df.where((df['month']>=12) | (df['month']<=2)).dropna(how='all')
df['yearmonth']=np.around(df['year'].values + df['month'].values/12 -1/12,3)
df.index = np.arange(len(df['iic']))
df_ = df[:]
#%%
# #Set indices where different cyclone is 
# fig,axz=plt.subplots(2,3,figsize=(16,9))
# for i,lower_window in enumerate(np.arange(0,6)):
#     ax = axz.flat[i]
#     bnds = np.where(np.diff(df['iic']) != 0)[0] + 1
#     lwr_bnds = np.insert(bnds,0,0)
#     upp_bnds = np.insert(bnds,-1,len(df['iic']))
#     bnds = np.array(list(zip(lwr_bnds,upp_bnds)))[:-1]
#     #Set index where lowest slp is observed per cyclone 
#     idx = df.groupby('iic')['slp'].idxmin().values[:-1]
#     #Set time window over which precipitation has to be summed together
#     idx_prec = np.array([[x-lower_window,x+1] for x in idx])
#     idx_mask = np.array([[x-lower_window,x+lower_window+1] for x in idx])
#     #If the precipitation time windows exceeds the lifetime of the cyclone we ignore it
#     mask = (bnds[:,0]<idx_mask[:,0])&(bnds[:,1]>=idx_mask[:,1])
#     #Compute wind speed at t_slp_min and compute precipitation sum over time window
#     ws_tmin = df['wsmean'].iloc[idx[mask]].values
#     preccum_tmin = np.array([df['precmean'].iloc[x[0]:x[1]].sum() for x in idx_prec[mask]])
#     ws_tmin_norm = (ws_tmin -ws_tmin.mean())/ws_tmin.std()
#     preccum_tmin_norm = (preccum_tmin - preccum_tmin.mean())/preccum_tmin.std()

#     assert(len(ws_tmin)==len(preccum_tmin))

#     slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(ws_tmin, preccum_tmin)
#     seaborn.regplot(x=ws_tmin,y=preccum_tmin, line_kws={"color": "red"},ax=ax)
#     #ax.scatter(ws_tmin,preccum_tmin)
#     ax.set_title('R squared = {:.2f}\nNumber of cyclones = {}\nTime Window {} hours'.format(r_value**2,len(ws_tmin),i*6+6))
#     ax.set_xlabel('Wind Speed')
#     ax.set_ylabel('Accumulated precipitation')
    
# fig.subplots_adjust(left=.02,right=.98,bottom=0.05,top=.95,hspace=.02,wspace=.02)
# fig.tight_layout()
# fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/figs/corr_p_v.png',dpi=300)

# # %%
# fig,axz=plt.subplots(2,3,figsize=(16,9))
# for i,lower_window in enumerate(np.arange(0,6)):
#     ax = axz.flat[i]
#     bnds = np.where(np.diff(df['iic']) != 0)[0] + 1
#     lwr_bnds = np.insert(bnds,0,0)
#     upp_bnds = np.insert(bnds,-1,len(df['iic']))
#     bnds = np.array(list(zip(lwr_bnds,upp_bnds)))[:-1]
#     #Set index where lowest slp is observed per cyclone 
#     idx = df.groupby('iic')['slp'].idxmin().values[:-1]
#     #Set time window over which precipitation has to be summed together
#     idx_prec = np.array([[x-lower_window,x+1] for x in idx])
#     idx_mask = np.array([[x-lower_window,x+lower_window+1] for x in idx])
#     #If the precipitation time windows exceeds the lifetime of the cyclone we ignore it
#     mask = (bnds[:,0]<idx_mask[:,0])&(bnds[:,1]>=idx_mask[:,1])

#     #Compute wind speed at t_slp_min and compute precipitation sum over time window
#     ws_tmin = df['wsmean'].iloc[idx[mask]].values
#     q_tmin = df['qmean'].iloc[idx[mask]].values
#     preccum_tmin = np.array([df['precmean'].iloc[x[0]:x[1]].sum() for x in idx_prec[mask]])
#     qcum_tmin = np.array([df['qmean'].iloc[x[0]:x[1]].sum() for x in idx_prec[mask]])
#     flux_tmin = ws_tmin*qcum_tmin
#     flux_tmin_norm = (flux_tmin -flux_tmin.mean())/flux_tmin.std()
#     preccum_tmin_norm = (preccum_tmin - preccum_tmin.mean())/preccum_tmin.std()

#     assert(len(ws_tmin)==len(preccum_tmin))

#     slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(flux_tmin, preccum_tmin)
#     seaborn.regplot(x=flux_tmin,y=preccum_tmin, line_kws={"color": "red"},ax=ax)
#     #ax.scatter(flux_tmin,preccum_tmin)
#     ax.set_title('R squared = {:.2f}\nNumber of cyclones = {}\nTime Window {} hours'.format(r_value**2,len(flux_tmin),i*6+6))
#     ax.set_xlabel('Wind Speed * TCWV')
#     ax.set_ylabel('Accumulated precipitation')
    
# fig.subplots_adjust(left=.02,right=.98,bottom=0.05,top=.95,hspace=.02,wspace=.02)
# fig.tight_layout()
# fig.savefig('/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/figs/corr_p_flux.png',dpi=300)
#%%
fig,axz=plt.subplots(3, 2, figsize=(16,9),subplot_kw={'projection': ccrs.PlateCarree()})
name = 'wsmean'
labels = np.arange(6,37,6)

for i in np.arange(0,6):
    bnds = np.where(np.diff(df_['iic']) != 0)[0] + 1
    lwr_bnds = np.insert(bnds,0,0)
    upp_bnds = np.insert(bnds,-1,len(df_['iic']))
    bnds = np.array(list(zip(lwr_bnds,upp_bnds)))[:-1]
    #Set index where lowest slp is observed per cyclone 
    idx = df_.groupby('iic')['slp'].idxmin().values[:-1]
    #Set time window over which precipitation has to be summed together
    idx_prec = np.array([[x-i,x+1] for x in idx])
    mask = (bnds[:,0]<idx_prec[:,0])&(bnds[:,1]>=idx_prec[:,1])
    df = df_[['iic','year','wsmean','lat','lon']].iloc[idx[mask]]
    df.index = df['iic'].values
    preccum_tmin = np.array([df_['precmean'].iloc[x[0]:x[1]].sum() for x in idx_prec[mask]])
    df['preccum']=preccum_tmin
    df['wsnorm']=(df['wsmean']-df['wsmean'].mean())/df['wsmean'].std()
    df['preccum_norm']=(df['preccum']-df['preccum'].mean())/df['preccum'].std()
    
    #prec_ext = df.where(df['preccum']>=df['preccum'].quantile(q=0.1)).dropna().sort_values('preccum',ascending=False)
    #ws_ext = df.where(df['wsmean']>=df['wsmean'].quantile(q=0.1)).dropna().sort_values('wsmean',ascending=False)
    
    prec_ext = df.where(df['preccum']>=df['preccum'].quantile(q=0.9)).dropna().sort_values('preccum',ascending=False)
    ws_ext = df.where(df['wsmean']>=df['wsmean'].quantile(q=0.9)).dropna().sort_values('wsmean',ascending=False)
    cmpnd_ext = np.intersect1d(prec_ext['iic'],ws_ext['iic'])
    cmpnd_df = df.where(df['iic'].isin(cmpnd_ext)).dropna()
    varname = cmpnd_df

    
    da_tot = []
    for j,var_ in enumerate([varname]):
        da = xr.DataArray(dims=['lat','lon'],
                      coords={'lat':sorted(df_['lat'].unique(),reverse=True),
                              'lon':sorted(df_['lon'].unique())})
        var = var_.groupby(['lat','lon'])[name].count()
        midx = var.index
        for lat,lon in midx:
            da.loc[lat,lon]=var.loc[lat].loc[lon]
        da.coords['lon'] = (da.coords['lon'] + 180) % 360 - 180
        da = da.sortby(da.lon)
        da_tot.append(da)
    
    ax = axz.flat[i]
    da=da_tot[j]
    #if i ==0:
    im=ax.pcolormesh(da.lon,da.lat,da)
    #else:
    #    ax.pcolormesh(da.lon,da.lat,da)
    ax.coastlines()
    fig.colorbar(im,ax=ax)
    ax.set_title('N:{} {}h-window\ncompund:{}'.format(len(varname),labels[i],len(cmpnd_df)))
#cbar_ax = fig.add_axes([0.05, 0.065, 0.9, 0.015])
#fig.colorbar(im,cax=cbar_ax,orientation='horizontal')
fig.suptitle(name)
fig.subplots_adjust(left=.05,right=.95,bottom=.05,top=.93,wspace=.1,hspace=.2)

    
        

# %%