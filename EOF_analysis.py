# %%
import xarray as xr
import numpy as np
from eofs.xarray import Eof
from sklearn.decomposition import PCA
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
from sklearn.preprocessing import normalize
import sys
path = '/storage/climatestor/PleioCEP/doensen/data/extracted/'
file_EOF = 'Z500_V300/EOF/Z500_EOF_ref.nc'


# %%
seasons = {'MAM': [3, 4, 5],
           'JJA': [6, 7, 8],
           'SON': [9, 10, 11],
           'DJF': [12, 1, 2]}
#months = np.arange(1,13)
ds_EOF = xr.open_dataset(path+file_EOF).squeeze()
ds_EOF.coords['lon'] = (ds_EOF.coords['lon'] + 180) % 360 - 180
ds_EOF = ds_EOF.sortby(ds_EOF.lon)

# %%
for key in seasons:
    da = ds_EOF.sel(lat=slice(80, 20), lon=slice(-90, 40)).Z
    da = da[(da.time.dt.year >= 105) & (da.time.dt.year <= 204)]
    da = da[da.time.dt.day == 1]
    #da = da[(da.time.dt.year>1000)&(da.time.dt.year<1200)]
    da = da[np.isin(da.time.dt.month, seasons[key])]
    # std = da.std(dim='time')
    # da = da/std
    time = da.time.values
    lat = da.lat.values
    lon = da.lon.values

    weights = np.cos(np.deg2rad(da.lat))
    weights = np.tile(weights.values.reshape(-1, 1), (1, len(da.lon)))
    solver = Eof(da, weights)
    eofs = solver.eofsAsCovariance(neofs=10)
    pcs = solver.pcs(npcs=10)
    exv = solver.varianceFraction(neigs=10)
    # eofs.to_dataset().to_netcdf(path+'EOF_season/EOFS_{}.nc'.format(key))
    # pcs.to_dataset().to_netcdf(path+'EOF_season/PCS_{}.nc'.format(key))

    fig, axz = plt.subplots(2, 3, subplot_kw={'projection': ccrs.Orthographic(-25, 50)},
                            figsize=(16, 10))
    for i, ax in enumerate(axz.flat):
        # ax.set_extent([-100,60,20,80])
        eof = eofs.sel(mode=i)
        # eof = eofs[i]
        ax.coastlines()
        ax.set_global()
        if i == 0:
            im = ax.pcolormesh(da.lon, da.lat, eof, vmin=-60, vmax=60,
                               transform=ccrs.PlateCarree(), cmap='bwr')
        else:
            ax.pcolormesh(da.lon, da.lat, eof, vmin=-60, vmax=60,
                          transform=ccrs.PlateCarree(), cmap='bwr')
        ax.set_title('Ex Var = {:.1f}%'.format(float(exv[i]*100)))
        #ax.set_title('Ex Var = {:.1f}%'.format(float(rat[i]*100)))
        # fig.colorbar(im)
    fig.subplots_adjust(left=0.05, right=0.95, hspace=0.1,
                        wspace=0.01, top=0.935, bottom=0.1)
    cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.025])
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    fig.suptitle('EOF {}'.format(key))
    continue
    fig.savefig(path+'EOF_season/figs/EOFS_{}.png'.format(key), dpi=300)
    # plt.close(fig)
    time = pcs.time.dt.year + pcs.time.dt.month/12 - 1/12
    fig, axz = plt.subplots(2, 3, figsize=(16, 9))
    for i, ax in enumerate(axz.flat):
        # ax.set_extent([-100,60,20,80])
        pc = pcs.sel(mode=i)
        ax.grid()
        ax.plot(time, pc)
        ax.set_ylim([-3, 3])
        ax.set_title('Mode {}'.format(i+1))
    fig.suptitle('Principal Component')
    fig.savefig(path+'EOF_season/figs/PCS_{}.png'.format(key), dpi=300)
    # plt.close(fig)

    # fig.colorbar(im)
    # fig.savefig('/storage/climatestor/PleioCEP/doensen/figs/pc.png',dpi=300)

# %%
"""
Project EOFS on longer dataset
"""
#Open seasonal avaeraged Z500 data for the EOF analysis for the CESM year 105 - 204
file_data = 'Z500_V300/EOF/Z500_EOF_ref.nc'
Z500_all = xr.open_dataset(
    path+file_data)
#Change longitudes to -180 to 180 degrees
Z500_all.coords['lon'] = (Z500_all.coords['lon'] + 180) % 360 - 180
Z500_all = Z500_all.sortby(Z500_all.lon)
#Select region relevant for North-Atlantic patterns
Z500_all = Z500_all.sel(lat=slice(80, 20), lon=slice(-90, 60)).squeeze().Z

#%%
#Select releavnt dates
Z500_sel = Z500_all[(Z500_all.time.dt.year >= 105) &
                    (Z500_all.time.dt.year <= 204)]
Z500_sel = Z500_sel[(Z500_sel.time.dt.month >= 12) |
                    (Z500_sel.time.dt.month <= 2)]
Z500_sel = Z500_sel[Z500_sel.time.dt.day == 1]
#Compute weights to adjust for latitude scaling
latitudes = Z500_sel.lat
weights = np.cos(np.deg2rad(latitudes)).values
weights_broadcasted = weights[np.newaxis, :, np.newaxis]
Z500_sel = weights_broadcasted * Z500_sel
#Flatten Z500 data so it becomes a 2D array with shape (time,space)
Z500 = Z500_sel.values.reshape(Z500_sel.values.shape[0], -1)
#Perform singular value decomposition
pc,std,eofs = np.linalg.svd(Z500, full_matrices=False)
var = std**2
#reshape eofs back to 3D array
eofs_reshaped = eofs.reshape((eofs.shape[0],) + Z500_sel.shape[1:])
for i in range(100):
    copy = Z500_sel.isel(time=0).squeeze().copy(deep=True)
    copy[:,:]=eofs_reshaped[0,:,:]
    #copy.to_netcdf(path+'EOF_season/EOF_mode_{:02d}.nc'.format(i+1))
    
# %%
#Plot leading EOFs 
fig, axz = plt.subplots(2,4,figsize=(16, 9), subplot_kw={'projection': ccrs.Orthographic(-25, 50)})
for i,ax in enumerate(axz.flat):
    ax = axz.flat[i]
    im = ax.pcolormesh(Z500_sel.lon, Z500_sel.lat,eofs_reshaped[i,:,:], vmin=-0.1, vmax=0.1,
                                transform=ccrs.PlateCarree(), cmap='bwr')
    ax.set_global()
    ax.coastlines()
    ex_var = (var[i]/var.sum())
    ax.set_title('Ex Var = {0:.1f}%'.format(ex_var*100))
fig.tight_layout()
fig.savefig(path+'/Z500_V300/EOF/EOF_DHF_105_204.png',dpi=300)



# %%
#Project EOF on longer monthly dataset

#Load full monthly Z500 data
Z500_full = xr.open_dataset(path+'Z500_V300/Z500_anom.nc').Z.squeeze()
Z500_full.coords['lon'] = (Z500_full.coords['lon'] + 180) % 360 - 180
Z500_full = Z500_full.sortby(Z500_full.lon)
Z500_full = Z500_full.sel(lat=slice(80,20),lon=slice(-90,40))
Z500_full = Z500_full[Z500_full.time.dt.day == 1]
Z500_full = Z500_full[(Z500_full.time.dt.year >= 5)]
latitudes = Z500_full.lat
weights = np.cos(np.deg2rad(latitudes)).values
weights_broadcasted = weights[np.newaxis, :, np.newaxis]
Z500_full_arr = Z500_full.values*weights_broadcasted
#%%
#project EOFs on data and adding them to lists
pcs = []; vars = []
for i in range(100):
    #select EOF
    eof = eofs_reshaped[i,:,:].reshape(-1)
    #Again reshape Z500 data to 2D array
    Z500_full_flat = Z500_full_arr.reshape(Z500_full_arr.shape[0],-1)
    #compute dot product to end up with PC component
    pc = np.dot(Z500_full_flat,eof)
    #normalize PC
    total_variance = np.sum(Z500_full_flat**2)
    exp_var = np.sum(pc**2)/total_variance
    mean = pc.mean(); std = pc.std()
    pc_norm = (pc-mean)/std
    #save PC to netcdf file
    ds = xr.Dataset(data_vars={'pc':('time',pc_norm)},
                    coords={'time':Z500_full.time.values},
                    attrs={'exp_var':exp_var})
    ds.to_netcdf(path+'EOF_season/PC_mode_{:02d}.nc'.format(i+1))
    print('{0:.1f}%'.format(exp_var*100))



#Z500_full_rs = Z500_full

# %%
#Project EOFs for other seasons
seasons = {'MAM': [3, 4, 5],
           'JJA': [6, 7, 8],
           'SON': [9, 10, 11],
           'DJF': [12, 1, 2]}

dic = {'MAM':[],
       'JJA':[],
       'SON':[],
       'DJF':[]}
for key in seasons:
    print(key)
    Z500_season = Z500_full[np.isin(Z500_full.time.dt.month,seasons[key])]
    for i in range(8):
        eof = eofs_reshaped[i,:,:].reshape(-1)
        Z500_season_flat = Z500_season.values.reshape(Z500_season.values.shape[0],-1)
        pc = np.dot(Z500_season_flat,eof)
        total_variance = np.sum(Z500_season_flat**2)
        exp_var = np.sum(pc**2)/total_variance
        mean = pc.mean(); std = pc.std()
        pc_norm = (pc-mean)/std
        dic[key].append(exp_var)
        # ds = xr.Dataset(data_vars={'pc':('time',pc_norm)},
        #                 coords={'time':Z500_season.time.values},
        #                 attrs={'exp_var':exp_var})
        # ds.to_netcdf(path+'EOF_season/PC_mode_{:02d}.nc'.format(i+1))
        print('{0:.1f}%'.format(exp_var*100))
    fig.savefig(path+'/Z500_V300/EOF/EOF_DHF_105_204.png',dpi=300)
#%%    
fig, axz = plt.subplots(2,4,figsize=(16, 9), subplot_kw={'projection': ccrs.Orthographic(-25, 50)})
for i,ax in enumerate(axz.flat):
    ax = axz.flat[i]
    im = ax.pcolormesh(Z500_sel.lon, Z500_sel.lat,eofs_reshaped[i,:,:], vmin=-0.1, vmax=0.1,
                                transform=ccrs.PlateCarree(), cmap='bwr')
    ax.set_global()
    ax.coastlines()
    ex_var = (var[i]/var.sum())
    MAM = dic['MAM'][i]; JJA = dic['JJA'][i]; SON = dic['SON'][i]; DJF = dic['DJF'][i]
    ax.set_title('MAM {:.1f}% JJA {:.1f}%\n SON {:.1f}% DJF {:.1f}%'\
        .format(MAM*100,JJA*100,SON*100,DJF*100))
fig.tight_layout()
sys.exit()
#fig.subplots_adjust(hspace=0.3)
fig.savefig(path+'/Z500_V300/EOF/EOF_season.png',dpi=300)
    
# %%
