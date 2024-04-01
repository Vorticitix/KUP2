#%%
import pandas as pd
#from my_tools import *
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import cftime
import cartopy.crs as ccrs
import sys
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import matplotlib as mpl
#%%

"""
Quickly plot 99th percentile of decadal/centural variable to get first glance of
long-term trends in cyclone statistics
"""

path = '/storage/climatestor/PleioCEP/doensen/data/'
file = 'fort_36_total_med_full.txt'

df = pd.read_csv(path+file).drop('Unnamed: 0',axis=1).reset_index(drop=True)
df['dec'] = ((df['year']-5)//10)
df['cen'] = ((df['year']-5)//100)*10

#%%
df_group_year = df.groupby('dec')
#%%
cz = df.groupby('iic').min()
cz_dec = cz.groupby('dec').quantile(q=.01)['cz']
cz_cen = cz.groupby('cen').quantile(q=.01)['cz']
ax = cz_dec.plot(); cz_cen.plot()
ax.set_title("1000 hPa geopotential height")
ax.set_ylabel("m")
#%%
slp = df.groupby('iic').min()
slp_dec = slp.groupby('dec').quantile(q=.01)['slp']
slp_cen = slp.groupby('cen').quantile(q=.01)['slp']
ax=slp_dec.plot(); slp_cen.plot()
ax.set_title("Sea level pressure")
ax.set_ylabel("hPa")
# %%
precmean = df.groupby('iic').max()
precmean_dec = precmean.groupby('dec').quantile(q=.99)['precmean']
precmean_cen = precmean.groupby('cen').quantile(q=.99)['precmean']
ax=precmean_dec.plot(); precmean_cen.plot()

precmax = df.groupby('iic').max()
precmax_dec = precmax.groupby('dec').quantile(q=.99)['precmax']
precmax_cen = precmax.groupby('cen').quantile(q=.99)['precmax']
precmax_dec.plot(); precmax_cen.plot()
ax.set_title("Average and max cyclone precipitation")
ax.set_ylabel("mm")
# %%
zrad = df.groupby('iic').max()
zrad_dec = zrad.groupby('dec').quantile(q=.99)['zrad']
zrad_cen = zrad.groupby('cen').quantile(q=.99)['zrad']
ax=zrad_dec.plot(); zrad_cen.plot()
ax.set_title("Cyclone radius")
ax.set_ylabel("km")
 #%%
# wsmean = df.groupby('iic').max()
# wsmean_dec = wsmean.groupby('dec').quantile(q=.99)['wsmean']
# wsmean_cen = wsmean.groupby('cen').quantile(q=.99)['wsmean']
# ax=wsmean_dec.plot(); wsmean_cen.plot()

wsmax = df.groupby('iic').max()
wsmax_dec = wsmax.groupby('dec').quantile(q=.99)['wsmax']
wsmax_cen = wsmax.groupby('cen').quantile(q=.99)['wsmax']
wsmax_dec.plot(); wsmax_cen.plot()
ax.set_title("Average and max cyclone wind (850 hPa)")
ax.set_ylabel("m/s")

#%%
"""
Rank storms based on cyclone intensity from different variables
"""

df_cz = df[['iic','cz']].groupby('iic',as_index=False).min().sort_values('cz').reset_index(drop=True)
df_cz['norm'] = (df_cz['cz'] - df_cz['cz'].mean())/df_cz['cz'].std()*-1

df_gz = df[['iic','gz']].groupby('iic',as_index=False).max().sort_values('gz',ascending=False).reset_index(drop=True)
df_gz['norm'] = (df_gz['gz'] - df_gz['gz'].mean())/df_gz['gz'].std()

df_wsmean = df[['iic','wsmean']].groupby('iic',as_index=False).max().sort_values('wsmean',ascending=False).reset_index(drop=True)
df_wsmean['norm'] = (df_wsmean['wsmean'] - df_wsmean['wsmean'].mean())/df_wsmean['wsmean'].std()

df_precmean = df[['iic','precmean']].groupby('iic',as_index=False).max().sort_values('precmean',ascending=False).reset_index(drop=True)
df_precmean['norm'] = (df_precmean['precmean'] - df_precmean['precmean'].mean())/df_precmean['precmean'].std()

#%%
df_groupby_min = df.groupby('iic',as_index=False).min()
df_groupby_min['count'] = df.groupby('iic').count().iloc[:,0].values
df_groupby_min = df_groupby_min.where(df_groupby_min['count']>8).dropna()

df_groupby_max = df.groupby('iic',as_index=False).max()
df_groupby_max['count'] = df.groupby('iic').count().iloc[:,0].values
df_groupby_max = df_groupby_max.where(df_groupby_max['count']>8).dropna()


df_cz = df_groupby_min[['iic','cz']].sort_values('cz').reset_index(drop=True)
df_cz = df_cz.sort_values('iic')

df_gz = df_groupby_max[['iic','gz']].sort_values('gz',ascending=False).reset_index(drop=True)
df_gz = df_gz.sort_values('iic')

df_wsmean = df_groupby_max[['iic','wsmean']].sort_values('wsmean',ascending=False).reset_index(drop=True)
df_wsmean = df_wsmean.sort_values('iic')

df_precmean = df_groupby_max[['iic','precmean']].sort_values('precmean',ascending=False).reset_index(drop=True)
df_precmean = df_precmean.sort_values('iic')

df_ext = pd.DataFrame(index=df_cz['iic'],
                      data={'cz':df_cz.index.values,
                            'gz':df_gz.index.values,
                            'wsmean':df_wsmean.index.values,
                            'precmean':df_precmean.index.values})
df_ext['mean'] = df_ext[['cz','gz','wsmean','precmean']].mean(axis=1)
df_ext = df_ext.sort_values('mean')

# %%

iic_ext = df_ext.iloc[:1000].index.values
year_ext = df[np.isin(df['iic'].values,iic_ext)].groupby('iic').min()['date']//10000

fig,ax = plt.subplots(figsize=(16,5))
ax.hist(year_ext,range=(0,3601),bins=3601)
ax.set_ylim([0,4])
ax.set_yticks([0,1,2,3,4])
ax.set_xlabel('Year')
ax.set_ylabel('Count per year')
ax.set_title("Occurences 500 most extreme cyclones")
# %%
"""
Plotting synoptic situation in location of interest
"""
df_select = df.where(df['iic']==df_ext.index.values[83]).dropna()
date = df_select['date']
year = [int(str(int(x)).zfill(8)[:4]) for x in date.values]
month = [int(str(int(x)).zfill(8)[4:6]) for x in date.values]
day = [int(str(int(x)).zfill(8)[6:]) for x in date.values]
hour = [x//10000 for x in df_select['idumi'].values]

cf_dates = []
for i in np.arange(len(year)):
    cf_date = cftime.DatetimeNoLeap(year[i],month[i],day[i],hour[i])
    cf_dates.append(cf_date)

if year[0] <= 604:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_0005_0604/'
elif year[0] <= 1204:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_0605_1204/'
elif year[0] <= 1804:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_1205_1804/'
elif year[0] <= 2404:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_1805_2404/'
elif year[0] <= 3004:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_2405_3004/'
else:
    path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_3005_3514/'
    
years = np.arange(5,3516,10)

closest_year = years[np.abs(years-year[0]).argmin()]

if year[0]>=closest_year:
    yy=closest_year
    eyy=closest_year+9
else:
    eyy=closest_year-1
    yy=closest_year-10
    
file='TRANS.3501BP.cam.h1.{:04d}_{:04d}_01-01.sel.nc'.format(yy,eyy)
    
ds = xr.open_dataset(path+file).squeeze()
ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180)).sortby('lon')

# %%
#Plot 1000 hPa geopotential height

ds_sel = ds.sel(time=cf_dates).sel(lat=slice(60,25),lon=slice(-15,60))
fig,ax = plt.subplots(figsize=(12,5),subplot_kw={'projection': ccrs.PlateCarree()})
ax.coastlines(zorder=3)
ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00'.format(cf_dates[0].day,
                                    cf_dates[0].month,
                                    cf_dates[0].year-1502,
                                    cf_dates[0].hour,
                                    int(df_select.iloc[0]['age'])))
scat = ax.scatter(x=df_select['lon'].iloc[0],
           y=df_select['lat'].iloc[0],
           marker="X",s=200,color='k',
           transform=ccrs.PlateCarree(),zorder=4)
lat = ds_sel.lat.values; lon = ds_sel.lon.values
PSL = ds_sel.sel(time=cf_dates[0]).PSL
WS = ds_sel.sel(time=cf_dates[0]).WS
cont = ax.contour(lon,lat,PSL,transform=ccrs.PlateCarree(),
                vmin=95000,vmax=105000,cmap="binary")
im = ax.pcolormesh(lon,lat,WS,transform=ccrs.PlateCarree(),
                vmin=0,vmax=50,cmap="hot_r")


plt.colorbar(im)
def animate(i):
    global cont, im, scat
    PSL = ds_sel.sel(time=cf_dates[i]).PSL
    WS = ds_sel.sel(time=cf_dates[i]).WS
    ax.clear()
    ax.coastlines(zorder=3)
    ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00'.format(cf_dates[i].day,
                                        cf_dates[i].month,
                                        cf_dates[i].year-1502,
                                        cf_dates[i].hour,))
    scat = ax.scatter(x=df_select['lon'].iloc[i],
               y=df_select['lat'].iloc[i],
               marker="X",s=200,color='k',
               transform=ccrs.PlateCarree(),zorder=4)
    for c in cont.collections:
        c.remove()
    
    cont = ax.contour(lon,lat,PSL,transform=ccrs.PlateCarree(),
                    vmin=97000,vmax=105000,cmap="binary")
    
    im = ax.pcolormesh(lon,lat,WS,transform=ccrs.PlateCarree(),
                    vmin=0,vmax=30,cmap="hot_r")
    
    return cont,im, scat
    
#ax.coastlines()   
ani = animation.FuncAnimation(fig,animate,frames=len(cf_dates),interval=600)
#fig.colorbar(cont)
ani.save('/storage/climatestor/PleioCEP/doensen/data/figs/WS_{}_{}.mp4'.format(df_select.iloc[0]['iic'],int(df_select.iloc[0]['year'])),
          dpi=300)
plt.show()

#%%

path='/storage/climatestor/PleioCEP/doensen/data/'
times  = pd.date_range(pd.Timestamp('1990-01-17T06:00:00.000000000'),
                 pd.Timestamp('1990-01-19T11:00:00.000000000'),
                 freq='1h')
ws850 = xr.open_dataset(path+'ws850.nc').rename({'__xarray_dataarray_variable__':'ws'})
SLP = xr.open_dataset(path+'slp.nc')
ws = ws850.sel(Time=times.values)
slp = SLP.sel(Time=times.values)

#ax.set_ylim([0, 90])

#%%

fig,ax=plt.subplots(figsize=(12, 5),subplot_kw={'projection': ccrs.PlateCarree()})
#ax.set_global()
im=ws.ws[0,:,:].plot.pcolormesh(
    ax=ax, transform=ccrs.PlateCarree(), x="XLONG", y="XLAT", add_colorbar=True,
    cmap='hot_r',vmax=50,
    cbar_kwargs={'label': '850 hPa Wind Speed'})

cont=slp.slp[0,:,:].plot.contour(
    ax=ax, transform=ccrs.PlateCarree(), x="XLONG", y="XLAT",
    cmap='binary',vmin=950,vmax=1050,levels=np.arange(950,1051,5),
    add_colorbar=False)


ax.set_extent([-15, 60, 25, 60])
ax.coastlines(zorder=3)
ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00'.format(times[0].day,
                                        times[0].month,
                                        times[0].year,
                                        times[0].hour))
# im = ax.contourf(lon,lat,ws,transform=ccrs.PlateCarree(),
#                 vmin=0,vmax=30,cmap="hot_r")

def animate(i):
    global cont, im, scat
    ws_ = ws.sel(Time=times[i]).ws
    slp_ = slp.sel(Time=times[i]).slp
    ax.clear()
    ax.set_extent([-15, 60, 25, 60])
    ax.coastlines(zorder=3)

    # scat = ax.scatter(x=df_select['lon'].iloc[i],
    #            y=df_select['lat'].iloc[i],
    #            marker="X",s=200,color='k',
    #            transform=ccrs.PlateCarree(),zorder=4)
    for c in cont.collections:
        c.remove()
    
    cont = slp_.plot.contour(
    ax=ax, transform=ccrs.PlateCarree(), x="XLONG", y="XLAT",
    cmap='binary',vmin=950,vmax=1050,levels=np.arange(950,1051,5),
    add_colorbar=False)
    
    im = ws_.plot.pcolormesh(
    ax=ax, transform=ccrs.PlateCarree(), x="XLONG", y="XLAT",
    cmap='hot_r',vmax=50,add_colorbar=False)
    ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00'.format(times[i].day,
                                            times[i].month,
                                            times[i].year,
                                            times[i].hour))
    
    return cont,im
    
#ax.coastlines()   
ani = animation.FuncAnimation(fig,animate,frames=len(times),interval=100)
#fig.colorbar(cont)
# ani.save('/storage/climatestor/PleioCEP/doensen/data/figs/WRF_{}_{}.mp4'.format(df_select.iloc[0]['iic'],int(df_select.iloc[0]['year'])),
#           dpi=300)
ani.save('/storage/climatestor/PleioCEP/doensen/data/figs/test.mp4',
          dpi=300)
plt.show()

# %%
#Plot 1000 hPa geopotential height
ds_sel = ds.sel(time=cf_dates).sel(lat=slice(60,25),lon=slice(-15,60))
fig,ax = plt.subplots(figsize=(12,5),subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([-15,60,60,25], crs=ccrs.PlateCarree())
ax.coastlines(zorder=3)
ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00     AGE = {}'.format(cf_dates[0].day,
                                    cf_dates[0].month,
                                    cf_dates[0].year,
                                    cf_dates[0].hour,
                                    int(df_select.iloc[0]['age'])))
scat = ax.scatter(x=df_select['lon'].iloc[0],
           y=df_select['lat'].iloc[0],
           marker="X",s=200,color='r',
           transform=ccrs.PlateCarree(),zorder=4)
lat = ds_sel.lat.values; lon = ds_sel.lon.values
PSL = ds_sel.sel(time=cf_dates[0]).PSL
PRECT = ds_sel.sel(time=cf_dates[0]).PRECT*86400*1000
cont = ax.contour(lon,lat,PSL,transform=ccrs.PlateCarree(),
                vmin=97000,vmax=105000,cmap="binary")
im = ax.contourf(lon,lat,PRECT,transform=ccrs.PlateCarree(),
                levels=[0.01,5,10,15,20,25,30,35,40,45,50],
                vmin=0.01,vmax=50,cmap="ocean_r",extend="max")
cbar = plt.colorbar(im,ticks=np.linspace(0,50,6))


def animate(i):
    global cont, im, scat
    PRECT = ds_sel.sel(time=cf_dates[i]).PRECT*86400*1000
    PSL = ds_sel.sel(time=cf_dates[i]).PSL
    ax.clear()
    ax.coastlines(zorder=3)
    ax.set_extent([-15,60,60,25], crs=ccrs.PlateCarree())
    ax.set_title('{:02d}-{:02d}-{:04d} {:02d}:00     AGE = {}'.format(cf_dates[i].day,
                                        cf_dates[i].month,
                                        cf_dates[i].year,
                                        cf_dates[i].hour,
                                        int(df_select.iloc[i]['age'])))
    scat = ax.scatter(x=df_select['lon'].iloc[i],
               y=df_select['lat'].iloc[i],
               marker="X",s=200,color='r',
               transform=ccrs.PlateCarree(),zorder=4)
    for c in cont.collections:
        c.remove()
    
    cont = ax.contour(lon,lat,PSL,transform=ccrs.PlateCarree(),
                    vmin=97000,vmax=105000,cmap="binary")
    
    im = ax.contourf(lon,lat,PRECT,transform=ccrs.PlateCarree(),
                levels=[0.01,5,10,15,20,25,30,35,40,45,50],
                vmin=0.01,vmax=50,cmap="ocean_r",extend="max")
    
    return cont,im, scat
    
#ax.coastlines()   
ani = animation.FuncAnimation(fig,animate,frames=len(cf_dates),interval=500)
ani.save('/storage/climatestor/PleioCEP/doensen/tmp_figs/cyclones_3500/PRECT_{}_{}.mp4'.format(df_select.iloc[0]['iic'],int(df_select.iloc[0]['year'])),
         dpi=300)

#fig.colorbar(cont)
plt.show()
    
#%%

import pylab as plt
import numpy
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = r"C:\some_path\ffmpeg.exe"   # if necessary

# Generate data for plotting
Lx = Ly = 3
Nx = Ny = 11
Nt = 20
x = numpy.linspace(0, Lx, Nx)
y = numpy.linspace(0, Ly, Ny)
x,y = numpy.meshgrid(x,y)
z0 = numpy.exp(-(x-Lx/2)**2-(y-Ly/2)**2)   # 2 dimensional Gaussian

def some_data(i):   # function returns a 2D data array
    return z0 * (i/Nt)

fig = plt.figure()
ax = plt.axes(xlim=(0, Lx), ylim=(0, Ly), xlabel='x', ylabel='y')

cvals = numpy.linspace(0,1,Nt+1)      # set contour values 
cont = plt.contourf(x, y, some_data(0), cvals)    # first image on screen
plt.colorbar()

# animation function
def animate(i):
    global cont
    z = some_data(i)
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(x, y, z, cvals)
    plt.title('t = %i:  %.2f' % (i,z[5,5]))
    return cont

anim = animation.FuncAnimation(fig, animate, frames=Nt, repeat=False)
#anim.save('animation.mp4', writer=animation.FFMpegWriter())
plt.show()



    