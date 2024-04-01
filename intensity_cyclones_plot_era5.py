#%%
%matplotlib qt
import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import sys
import cartopy.crs as ccrs
import colormaps as cmaps
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from my_tools import FigureTemplate
from scipy.stats import ttest_ind
import matplotlib as mpl
import string
mpl.rcParams['hatch.linewidth']=.3
#%%
plotmax = True
if plotmax:
    maxstring = '_max'
else:   
    maxstring = ''
    
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/ERA5/'
file = '{}_{}_top_{:04d}_all_{}.nc'

regions = ['emed','cmed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
suptitles = ['Central\nMediterranean','Eastern\nMediterranean']
titles = ['Most extreme cyclones with respect to precipitation',
          'Most extreme cyclones with respect to wind speed',
          'Compound Events']
axis_labels = ['CESM','ERA5']
compare_med = False
nrows=2
ncols=3
ranks = [10,20]
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
letters = list(string.ascii_lowercase)


#stat_PRECT = ttest_ind(ds_wmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]    
for k,region in enumerate(regions):
    for rank in ranks:
        obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
        obj.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05,wspace=0,hspace=0)
        fig = obj.fig; axz = axz = obj.axz
        letter_index = 0
        for i in range(2):
            for j in range(3):
                ds_cesm = xr.open_dataset(path+file.format(region,events[j],rank,'CESM_modern')).isel(time=8)
                ds_era5 = xr.open_dataset(path+file.format(region,events[j],rank,'ERA5')).isel(time=8)
                ds_era5 = ds_era5.rename({'msl':'PSL','prec':'PRECT'})
                
                stat_WS = ttest_ind(ds_cesm['WS'].values,ds_era5['WS'].values,equal_var=False)[1]
                if i == 0:
                    ds = ds_cesm.mean(dim='rank')

                    lats=ds['lat'].values; lons = ds['lon'].values
                    ax = axz[i,j]
                    ax.set_xlabel('Degrees')
                    ax.set_ylabel(axis_labels[i])
                    ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                    im =ax.contourf(lons,lats,ds['WS'].values,levels=np.arange(0,30.1,2.5),
                                cmap=cmap1)
                    ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                    cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                                cmap='Greys')
                    ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                    if j + 1 == ncols:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        fig.colorbar(im, cax=cax, orientation='vertical',label='850 hPa Wind Speed [m/s]',
                                    ticks=np.arange(0,30.1,2.5))
                    ax.set_title(titles[j])
                
                else:
                    ds = ds_era5.mean(dim='rank')

                    lats=ds['lat'].values; lons = ds['lon'].values
                    ax = axz[i,j]
                    ax.set_xlabel('Degrees')
                    ax.set_ylabel(axis_labels[i])
                    ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                    im =ax.contourf(lons,lats,ds['WS'].values,levels=np.arange(0,30.1,2.5),
                                cmap=cmap1)
                    ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                    cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                                cmap='Greys')
                    ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                    if j + 1 == ncols:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        fig.colorbar(im, cax=cax, orientation='vertical',label='850 hPa Wind Speed [m/s]',
                                    ticks=np.arange(0,30.1,2.5))
                ax.label_outer()
                ax.text(-0.06, 1.025, letters[letter_index]+')', transform=ax.transAxes, 
                        size=16, weight='bold')
                letter_index += 1
        fig.suptitle('{} most extreme cyclones in the {} '.format(rank,suptitles[k]),fontsize=16)
        fig.tight_layout()
        if compare_med==True:
            fig.savefig(path+'figs/rank_{}_spatial_WS_comparison{}_{}.png'.format(rank,maxstring,region),dpi=400)
        elif compare_med==False:
            fig.savefig(path+'figs/rank_{}_spatial_WS_comparison_reg{}_{}.png'.format(rank,maxstring,region),dpi=400)

#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/ERA5/'
file = '{}_{}_top_{:04d}_all_{}.nc'

regions = ['emed','cmed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
suptitles = ['Central\nMediterranean','Eastern\nMediterranean']
titles = ['Most extreme cyclones with respect to precipitation',
          'Most extreme cyclones with respect to wind speed',
          'Compound Events']
axis_labels = ['CESM','ERA5']
compare_med = False
nrows=2
ncols=3
ranks = [10,20]
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
letters = list(string.ascii_lowercase)


#stat_PRECT = ttest_ind(ds_wmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]    
for k,region in enumerate(regions):
    for rank in ranks:
        obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
        obj.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05,wspace=0,hspace=0)
        fig = obj.fig; axz = axz = obj.axz
        letter_index = 0
        for i in range(2):
            for j in range(3):
                ds_cesm = xr.open_dataset(path+file.format(region,events[j],rank,'CESM_modern')).isel(time=8)
                ds_era5 = xr.open_dataset(path+file.format(region,events[j],rank,'ERA5')).isel(time=8)
                ds_era5 = ds_era5.rename({'msl':'PSL','prec':'PRECT'})
                #ds_era5['PRECT'] = #ds_era5['PRECT']*1000*21600
                
                stat_PRECT = ttest_ind(ds_cesm['PRECT'].values,ds_era5['PRECT'].values,equal_var=False)[1]
                if i == 0:
                    ds = ds_cesm.mean(dim='rank')

                    lats=ds['lat'].values; lons = ds['lon'].values
                    ax = axz[i,j]
                    ax.set_xlabel('Degrees')
                    ax.set_ylabel(axis_labels[i])
                    ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                    im = ax.contourf(lons,lats,(ds['PRECT']*1000*21600).values,
                                cmap=cmap2,levels=[1,2,3,4,5,7.5,10,12.5,15])
                    cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                                cmap='Greys')
                    ax.contourf(lons,lats,stat_PRECT<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                    ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                    if j + 1 == ncols:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        fig.colorbar(im, cax=cax, orientation='vertical',label='850 hPa Wind Speed [m/s]',
                                    ticks=np.arange(0,30.1,2.5))
                    ax.set_title(titles[j])
                
                else:
                    ds = ds_era5.mean(dim='rank')

                    lats=ds['lat'].values; lons = ds['lon'].values
                    ax = axz[i,j]
                    ax.set_xlabel('Degrees')
                    ax.set_ylabel(axis_labels[i])
                    ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                    im = ax.contourf(lons,lats,(ds['PRECT']*21600).values,
                                cmap=cmap2,levels=[1,2,3,4,5,7.5,10,12.5,15])
                    cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                                cmap='Greys')
                    ax.contourf(lons,lats,stat_PRECT<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                    ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                    if j + 1 == ncols:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        fig.colorbar(im, cax=cax, orientation='vertical',label='850 hPa Wind Speed [m/s]',
                                    ticks=np.arange(0,30.1,2.5))
                ax.label_outer()
                ax.text(-0.06, 1.025, letters[letter_index]+')', transform=ax.transAxes, 
                        size=16, weight='bold')
                letter_index += 1
        fig.suptitle('{} most extreme cyclones in the {} '.format(rank,suptitles[k]),fontsize=16)
        fig.tight_layout()
        if compare_med==True:
            fig.savefig(path+'figs/rank_{}_spatial_PRECT_comparison{}_{}.png'.format(rank,maxstring,region),dpi=400)
        elif compare_med==False:
            fig.savefig(path+'figs/rank_{}_spatial_PRECT_comparison_reg{}_{}.png'.format(rank,maxstring,region),dpi=400)
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
file = '{}_{}_top_{:04d}_all.nc'

regions = ['emed','cmed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
axis_labels = ['Central\nMediterranean','Eastern\nMediterranean']
titles = ['Most extreme cyclones with respect to precipitation',
          'Most extreme cyclones with respect to wind speed',
          'Compound Events']

nrows=2
ncols=3
ranks = [10,100,1000]
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
compare_med = False

#stat_PRECT = ttest_ind(ds_wmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]    
for rank in ranks:
    obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
    obj.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05,wspace=0.05,hspace=0.05)
    fig = obj.fig; axz = axz = obj.axz
    letter_index = 0

    if compare_med==False:
        stat_extra_emed = ttest_ind(xr.open_dataset(path+file.format('emed',events[0],rank)).isel(time=8)['PRECT'].values,
                                    xr.open_dataset(path+file.format('emed',events[1],rank)).isel(time=8)['PRECT'].values,
                                    equal_var=False)[1]
        stat_extra_cmed = ttest_ind(xr.open_dataset(path+file.format('cmed',events[0],rank)).isel(time=8)['PRECT'].values,
                                    xr.open_dataset(path+file.format('cmed',events[1],rank)).isel(time=8)['PRECT'].values,
                                    equal_var=False)[1]
    for i in range(2):
        for j in range(3):
            ds_cmed = xr.open_dataset(path+file.format('cmed',events[j],rank)).isel(time=8)
            ds_emed = xr.open_dataset(path+file.format('emed',events[j],rank)).isel(time=8)
            stat_PRECT = ttest_ind(ds_cmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]
            if i == 0:
                ds = ds_cmed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,(ds['PRECT']*1000*21600).values,
                            cmap=cmap2,levels=[1,2,3,4,5,7.5,10,12.5,15])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                if compare_med==True:
                    ax.contourf(lons,lats,stat_PRECT<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                elif compare_med==False and (j==0 or j==1):
                    ax.contourf(lons,lats,stat_extra_cmed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                if j +1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='Precipitation [mm/6h]',
                                ticks=[1,2,3,4,5,7.5,10,12.5,15])

                ax.set_title(titles[j])
            
            else:
                ds = ds_emed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,(ds['PRECT']*1000*21600).values,
                            cmap=cmap2,levels=[1,2,3,4,5,7.5,10,12.5,15])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                if compare_med==True:
                    ax.contourf(lons,lats,stat_PRECT<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                elif compare_med==False and (j==0 or j==1):
                    ax.contourf(lons,lats,stat_extra_emed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                ax.clabel(cs,fmt='%d',colors='k',fontsize=8)
                if j +1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='Precipitation [mm/6h]',
                                ticks=[1,2,3,4,5,7.5,10,12.5,15])
            ax.label_outer()
            ax.text(-0.06, 1.025, letters[letter_index]+')', transform=ax.transAxes, 
                    size=16, weight='bold')
            letter_index += 1
    fig.suptitle('{} most extreme cyclones in the Mediterranean'.format(rank),fontsize=16)
    fig.tight_layout()
    if compare_med==True:
        fig.savefig(path+'figs/rank_{}_spatial_PRECT_comparison{}.png'.format(rank,maxstring),dpi=400)
    else:
        fig.savefig(path+'figs/rank_{}_spatial_PRECT_comparison_reg{}.png'.format(rank,maxstring),dpi=400)

#%%
"""
Comparison Central vs Eastern Mediterranean Surface Wind speed
"""
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
file = '{}_{}_top_{:04d}_all.nc'
compare_med = True
regions = ['cmed','emed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
y_labels = ['Max. precipitation [mm/6h]','Max. surface wind speed [m/s]','Min. sea level pressure [hPa]']
titles = ['Most extreme cyclones\nwith respect to precipitation','Most extreme cyclones\nwith respect to wind speed','Compound Events']
ranks = [10,100,1000]
vars = ['PRECT','WS','PSL']

nrows=3
ncols=3
#fig,axz = plt.subplots(nrows,ncols,figsize=(16,10))
obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
obj.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05,wspace=0.05,hspace=0.05)
letter_index = 0

fig = obj.fig; axz = axz = obj.axz
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
for i in range(nrows):
    for j in range(ncols):
        ax = axz[i,j]
        if compare_med==False:
            ds_cmed_10 = xr.open_dataset(path+file.format('cmed',events[j],10)).sel(time=slice(-5,5))
            ds_cmed_100 = xr.open_dataset(path+file.format('cmed',events[j],100)).sel(time=slice(-5,5))
            ds_cmed_1000 = xr.open_dataset(path+file.format('cmed',events[j],1000)).sel(time=slice(-5,5))
            ds_emed_10 = xr.open_dataset(path+file.format('emed',events[j],10)).sel(time=slice(-5,5))
            ds_emed_100 = xr.open_dataset(path+file.format('emed',events[j],100)).sel(time=slice(-5,5))
            ds_emed_1000 = xr.open_dataset(path+file.format('emed',events[j],1000)).sel(time=slice(-5,5))
            if vars[i]=='PSL':
                p_value_cmed_10 = ttest_ind(ds_cmed_10[vars[i]].min(['lat','lon']).values,ds_cmed_100[vars[i]].min(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_emed_1000 = ttest_ind(ds_cmed_100[vars[i]].min(['lat','lon']).values,ds_cmed_1000[vars[i]].min(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_emed_10 = ttest_ind(ds_emed_10[vars[i]].min(['lat','lon']).values,ds_emed_100[vars[i]].min(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_cmed_1000 = ttest_ind(ds_emed_100[vars[i]].min(['lat','lon']).values,ds_emed_1000[vars[i]].min(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05 
            else:
                p_value_cmed_10 = ttest_ind(ds_cmed_10[vars[i]].max(['lat','lon']).values,ds_cmed_100[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_emed_1000 = ttest_ind(ds_cmed_100[vars[i]].max(['lat','lon']).values,ds_cmed_1000[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_emed_10 = ttest_ind(ds_emed_10[vars[i]].max(['lat','lon']).values,ds_emed_100[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
                p_value_cmed_1000 = ttest_ind(ds_emed_100[vars[i]].max(['lat','lon']).values,ds_emed_1000[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05 

        for rank in ranks:
            ds_cmed = xr.open_dataset(path+file.format('cmed',events[j],rank)).sel(time=slice(-5,5))
            ds_emed = xr.open_dataset(path+file.format('emed',events[j],rank)).sel(time=slice(-5,5))
            ds_cmed_mean = ds_cmed.mean(dim='rank')
            ds_emed_mean = ds_emed.mean(dim='rank')
            time = ds_cmed['time']*6
            # if j == 0:
            #     ds = ds_cmed.mean(dim='rank')
            # else:
            #     ds = ds_emed.mean(dim='rank')
            if vars[i] == 'PSL':
                p_value = ttest_ind(ds_cmed[vars[i]].min(['lat','lon']).values,ds_emed[vars[i]].min(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
            else:
                p_value = ttest_ind(ds_cmed[vars[i]].max(['lat','lon']).values,ds_emed[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
            if vars[i] == 'PSL':
                var_cmed = ds_cmed_mean['PSL'].min(['lat','lon'])/100
                var_emed = ds_emed_mean['PSL'].min(['lat','lon'])/100
                ax.set_ylim([980,1015])
            elif vars[i] == 'PRECT':
                var_cmed = ds_cmed_mean['PRECT'].max(['lat','lon'])*1000*21600
                var_emed = ds_emed_mean['PRECT'].max(['lat','lon'])*1000*21600
                ax.set_ylim([0,15])
            else:
                var_cmed = ds_cmed_mean['WS'].max(['lat','lon'])
                var_emed = ds_emed_mean['WS'].max(['lat','lon'])
                ax.set_ylim([10,30])
            if rank == 100:
                linestyle= '-'
                ax.plot(time.values,var_cmed.values,color='k',label='Central Med',linestyle=linestyle)
                ax.plot(time.values,var_emed.values,color='r',label='Eastern Med',linestyle=linestyle)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
            elif rank == 10:
                linestyle = '--'
                ax.plot(time.values,var_cmed.values,color='k',linestyle=linestyle,linewidth=0.7)
                ax.plot(time.values,var_emed.values,color='r',linestyle=linestyle,linewidth=0.7)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
                else:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value_cmed_10[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value_emed_10[i]]

            else:
                linestyle = ':'
                ax.plot(time.values,var_cmed.values,color='k',linestyle=linestyle,linewidth=0.7)
                ax.plot(time.values,var_emed.values,color='r',linestyle=linestyle,linewidth=0.7)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
                else:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value_cmed_1000[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value_emed_1000[i]]
            #lats=ds['lat'].values; lons = ds['lon'].values

        ax.set_xlabel('Hours')
        ax.set_ylabel(y_labels[i])
        ax.margins(0)
        if i==0:
            ax.set_title(titles[j])
        
        # ax2.set_ylabel('Sea level pressure [hPa]')
        # ax.set_ylim([0,35])
        # ax2.set_ylim([980,1015])
        ax.set_xticks(np.arange(-30,30.1,6))
        # ax2.set_xticks(np.arange(-48,48.1,6))

        if (i==0)&(j==0):
            ax.legend()
        ax.axvline(0,color='k')
        ax.label_outer()
        ax.grid()
        ax.text(-0.06, 1.075, letters[letter_index]+')', transform=ax.transAxes, 
                    size=16, weight='bold')
        letter_index += 1

        #ax.coastlines()
        # xr.plot.contour(ds['WS'],ax=ax)
        
#fig.suptitle('100 most extreme cyclones',fontsize=16)
fig.tight_layout()
if compare_med==True:
    fig.savefig(path+'figs/surface_evolution_compare_med{}.png'.format(maxstring),dpi=400)
else:
    fig.savefig(path+'figs/surface_evolution_compare_number{}.png'.format(maxstring),dpi=400)
#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
file = '{}_{}_top_{:04d}_all.nc'

regions = ['emed','cmed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
axis_labels = ['Central\nMediterranean','Eastern\nMediterranean']
titles = ['Most extreme cyclones with respect to precipitation',
          'Most extreme cyclones with respect to wind speed',
          'Compound Events']

nrows=2
ncols=3
ranks = [10,100,1000]
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
compare_med = True

#stat_PRECT = ttest_ind(ds_wmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]    
for rank in ranks:
    obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
    obj.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05,wspace=0.05,hspace=0.05)
    letter_index = 0
    fig = obj.fig; axz = axz = obj.axz
    if compare_med==False:
        stat_extra_emed = ttest_ind(xr.open_dataset(path+file.format('emed',events[0],rank)).isel(time=8)['WS300'].values,
                                    xr.open_dataset(path+file.format('emed',events[1],rank)).isel(time=8)['WS300'].values,
                                    equal_var=False)[1]
        stat_extra_cmed = ttest_ind(xr.open_dataset(path+file.format('cmed',events[0],rank)).isel(time=8)['WS300'].values,
                                    xr.open_dataset(path+file.format('cmed',events[1],rank)).isel(time=8)['WS300'].values,
                                    equal_var=False)[1]
    for i in range(2):
        for j in range(3):
            ds_cmed = xr.open_dataset(path+file.format('cmed',events[j],rank)).isel(time=8)
            ds_emed = xr.open_dataset(path+file.format('emed',events[j],rank)).isel(time=8)
            stat_WS = ttest_ind(ds_cmed['WS300'].values,ds_emed['WS300'].values,equal_var=False)[1]
            if i == 0:
                ds = ds_cmed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,ds['WS300'].values,levels=np.arange(35,50.1,2.5),
                                    cmap='Greens',vmin=20,vmax=50)
                # if compare_med==True:
                #     ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                # elif compare_med==False and (j==0 or j==1):
                #     ax.contourf(lons,lats,stat_extra_cmed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                cs2 = ax.contour(lons,lats,ds['T'].values,levels=[-8,-6,-4,-2,2,4,6,8],
                                    cmap='seismic',linewidths=2,linestyles='dashed')
                ax.clabel(cs2,fmt='%d',colors='k',fontsize=8)
                if j + 1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='300 hPa Wind Speed [m/s]',
                                ticks=np.arange(30,50.1,2.5))
                ax.set_title(titles[j])
            
            else:
                ds = ds_emed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,ds['WS300'].values,levels=np.arange(35,50.1,2.5),
                                    cmap='Greens',vmin=20,vmax=50)
                # if compare_med==True:
                #     ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                # elif compare_med==False and (j==0 or j==1):
                #     ax.contourf(lons,lats,stat_extra_emed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                cs2 = ax.contour(lons,lats,ds['T'].values,levels=[-8,-6,-4,-2,2,4,6,8],
                                    cmap='seismic',linewidths=2,linestyles='dashed')
                ax.clabel(cs2,fmt='%d',colors='k',fontsize=8)
                if j + 1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='300 hPa Wind Speed [m/s]',
                                ticks=np.arange(30,50.1,2.5))
                ax.set_title(titles[j])
            ax.label_outer()
            ax.text(-0.06, 1.025, letters[letter_index]+')', transform=ax.transAxes, 
                    size=16, weight='bold')
            letter_index += 1

    fig.suptitle('{} most extreme cyclones in the Mediterranean'.format(rank),fontsize=16)
    fig.tight_layout()
    if compare_med==True:
        fig.savefig(path+'figs/rank_{}_spatial_V300_T850_comparison{}.png'.format(rank,maxstring),dpi=400)
    else:
        fig.savefig(path+'figs/rank_{}_spatial_V300_T850_comparison_reg{}.png'.format(rank,maxstring),dpi=400)

#%%
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
file = '{}_{}_top_{:04d}_all.nc'

regions = ['emed','cmed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
axis_labels = ['Central\nMediterranean','Eastern\nMediterranean']
titles = ['Most extreme cyclones with respect to precipitation',
          'Most extreme cyclones with respect to wind speed',
          'Compound Events']

nrows=2
ncols=3
ranks = [10,100,1000]
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
compare_med = False

#stat_PRECT = ttest_ind(ds_wmed['PRECT'].values,ds_emed['PRECT'].values,equal_var=False)[1]    
for rank in ranks:
    obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
    obj.subplots_adjust(left=0.05,right=0.95,top=0.92,bottom=0.05,wspace=0.05,hspace=0.05)
    letter_index = 0
    fig = obj.fig; axz = axz = obj.axz
    if compare_med==False:
        stat_extra_emed = ttest_ind(xr.open_dataset(path+file.format('emed',events[0],rank)).isel(time=8)['RWP'].values,
                                    xr.open_dataset(path+file.format('emed',events[1],rank)).isel(time=8)['RWP'].values,
                                    equal_var=False)[1]
        stat_extra_cmed = ttest_ind(xr.open_dataset(path+file.format('cmed',events[0],rank)).isel(time=8)['RWP'].values,
                                    xr.open_dataset(path+file.format('cmed',events[1],rank)).isel(time=8)['RWP'].values,
                                    equal_var=False)[1]
    for i in range(2):
        for j in range(3):
            ds_cmed = xr.open_dataset(path+file.format('cmed',events[j],rank)).isel(time=8)
            ds_emed = xr.open_dataset(path+file.format('emed',events[j],rank)).isel(time=8)
            stat_WS = ttest_ind(ds_cmed['RWP'].values,ds_emed['RWP'].values,equal_var=False)[1]
            if i == 0:
                ds = ds_cmed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,ds['Z'].values,levels=np.arange(-250,250.1,25),
                                    cmap='bwr',vmin=-250,vmax=250)
                # if compare_med==True:
                #     ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                # elif compare_med==False and (j==0 or j==1):
                #     ax.contourf(lons,lats,stat_extra_cmed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                cs2 = ax.contour(lons,lats,ds['RWP'].values,levels=np.arange(24,40.1,2),
                                    cmap='Purples',linewidths=2,linestyles='dashed')
                ax.clabel(cs2,fmt='%d',colors='k',fontsize=8)
                if j + 1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='Z500 anomalies [gpm]',
                                ticks=np.arange(-250,250.1,50))
                ax.set_title(titles[j])
            
            else:
                ds = ds_emed.mean(dim='rank')

                lats=ds['lat'].values; lons = ds['lon'].values
                ax = axz[i,j]
                ax.set_xlabel('Degrees')
                ax.set_ylabel(axis_labels[i])
                ax.scatter(0,0,marker='X',zorder=3,s=80,color='k')
                im = ax.contourf(lons,lats,ds['Z'].values,levels=np.arange(-250,250.1,25),
                                    cmap='bwr',vmin=-250,vmax=250)
                # if compare_med==True:
                #     ax.contourf(lons,lats,stat_WS<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                # elif compare_med==False and (j==0 or j==1):
                #     ax.contourf(lons,lats,stat_extra_emed<0.05,levels=[0,0.05,1],colors='none',hatches=['','.'])
                cs = ax.contour(lons,lats,(ds['PSL']/100).values,levels=np.arange(980,1041,5),
                            cmap='Greys')
                cs2 = ax.contour(lons,lats,ds['RWP'].values,levels=np.arange(24,40.1,2),
                                    cmap='Purples',linewidths=2,linestyles='dashed')
                ax.clabel(cs2,fmt='%d',colors='k',fontsize=8)
                if j + 1 == ncols:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    fig.colorbar(im, cax=cax, orientation='vertical',label='Z500 anomalies [gpm]',
                                ticks=np.arange(-250,250.1,50))
                ax.set_title(titles[j])
            ax.label_outer()
            ax.text(-0.06, 1.025, letters[letter_index]+')', transform=ax.transAxes, 
                    size=16, weight='bold')
            letter_index += 1
    fig.suptitle('{} most extreme cyclones in the Mediterranean'.format(rank),fontsize=16)
    fig.tight_layout()
    if compare_med==True:
        fig.savefig(path+'figs/rank_{}_spatial_RWP_Z500_comparison{}.png'.format(rank,maxstring),dpi=400)
    else:
        fig.savefig(path+'figs/rank_{}_spatial_RWP_Z500_comparison_reg{}.png'.format(rank,maxstring),dpi=400)
        
#%%

#%%
"""
Comparison Central vs Eastern Mediterranean Surface Wind speed
"""
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'
file = '{}_{}_top_{:04d}_all.nc'
compare_med = True
regions = ['cmed','emed']
tops = [10,100,1000]
if plotmax==True:
    events = ['precmax','wsmax','compoundmax']
else:
    events = ['prec','ws','compound']
y_labels = ['300 hPa wind speed [m/s]','Rossby wave packet amplitude (m/s)']
titles = ['Most extreme cyclones\nwith respect to precipitation','Most extreme cyclones\nwith respect to wind speed','Compound Events']
ranks = [10,100,1000]
vars = ['WS300','RWP']

nrows=2
ncols=3
#fig,axz = plt.subplots(nrows,ncols,figsize=(16,10))
obj = FigureTemplate(nrows=nrows,ncols=ncols,figsize=(16,9))
obj.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05,wspace=0.05,hspace=0.05)
letter_index = 0
fig = obj.fig; axz = axz = obj.axz
manual_positions =  [(1, 1)]
cmap1=cmaps.NCV_jet
cmap2=cmaps.GMT_cool
for i in range(nrows):
    for j in range(ncols):
        ax = axz[i,j]
        if compare_med==False:
            ds_cmed_10 = xr.open_dataset(path+file.format('cmed',events[j],10)).sel(time=slice(-5,5))
            ds_cmed_100 = xr.open_dataset(path+file.format('cmed',events[j],100)).sel(time=slice(-5,5))
            ds_cmed_1000 = xr.open_dataset(path+file.format('cmed',events[j],1000)).sel(time=slice(-5,5))
            ds_emed_10 = xr.open_dataset(path+file.format('emed',events[j],10)).sel(time=slice(-5,5))
            ds_emed_100 = xr.open_dataset(path+file.format('emed',events[j],100)).sel(time=slice(-5,5))
            ds_emed_1000 = xr.open_dataset(path+file.format('emed',events[j],1000)).sel(time=slice(-5,5))

            p_value_cmed_10 = ttest_ind(ds_cmed_10[vars[i]].max(['lat','lon']).values,ds_cmed_100[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
            p_value_emed_1000 = ttest_ind(ds_cmed_100[vars[i]].max(['lat','lon']).values,ds_cmed_1000[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
            p_value_emed_10 = ttest_ind(ds_emed_10[vars[i]].max(['lat','lon']).values,ds_emed_100[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05
            p_value_cmed_1000 = ttest_ind(ds_emed_100[vars[i]].max(['lat','lon']).values,ds_emed_1000[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05 

        for rank in ranks:
            ds_cmed = xr.open_dataset(path+file.format('cmed',events[j],rank)).sel(time=slice(-5,5))
            ds_emed = xr.open_dataset(path+file.format('emed',events[j],rank)).sel(time=slice(-5,5))
            ds_cmed_mean = ds_cmed.mean(dim='rank')
            ds_emed_mean = ds_emed.mean(dim='rank')
            times = ds_cmed['time'].values
            # if j == 0:
            #     ds = ds_cmed.mean(dim='rank')
            # else:
            #     ds = ds_emed.mean(dim='rank')
            p_value = ttest_ind(ds_cmed[vars[i]].max(['lat','lon']).values,ds_emed[vars[i]].max(['lat','lon']).values,equal_var=False,nan_policy='omit')[1]<0.05

            if vars[i] == 'RWP':
                var_cmed = ds_cmed_mean['RWP'].max(['lat','lon'])
                var_emed = ds_emed_mean['RWP'].max(['lat','lon'])
                ax.set_ylim([20,50])
            elif vars[i] == 'WS300':
                var_cmed = ds_cmed_mean['WS300'].max(['lat','lon'])
                var_emed = ds_emed_mean['WS300'].max(['lat','lon'])
                ax.set_ylim([30,60])
            time = ds_cmed['time']*6
            if rank == 100:
                linestyle= '-'
                ax.plot(time.values,var_cmed.values,color='k',label='Central Med',linestyle=linestyle)
                ax.plot(time.values,var_emed.values,color='r',label='Eastern Med',linestyle=linestyle)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
            elif rank == 10:
                linestyle = '--'
                ax.plot(time.values,var_cmed.values,color='k',linestyle=linestyle,linewidth=0.7)
                ax.plot(time.values,var_emed.values,color='r',linestyle=linestyle,linewidth=0.7)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
                else:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value_cmed_10[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value_emed_10[i]]

            else:
                linestyle = ':'
                ax.plot(time.values,var_cmed.values,color='k',linestyle=linestyle,linewidth=0.7)
                ax.plot(time.values,var_emed.values,color='r',linestyle=linestyle,linewidth=0.7)
                if compare_med==True:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value[i]]
                else:
                    [ax.scatter(time.values[i],var_cmed.values[i],marker='X',zorder=3,color='k',s=20) for i in range(len(time)) if p_value_cmed_1000[i]]
                    [ax.scatter(time.values[i],var_emed.values[i],marker='X',zorder=3,color='r',s=20) for i in range(len(time)) if p_value_emed_1000[i]]
            #lats=ds['lat'].values; lons = ds['lon'].values

        ax.set_xlabel('Hours')
        ax.set_ylabel(y_labels[i])
        if i==0:
            ax.set_title(titles[j])
        
        # ax2.set_ylabel('Sea level pressure [hPa]')
        # ax.set_ylim([0,35])
        # ax2.set_ylim([980,1015])
        ax.margins(0)
        ax.set_xticks(np.arange(-30,30.1,6))
        # ax2.set_xticks(np.arange(-48,48.1,6))

        if (i==0)&(j==0):
            ax.legend()
        ax.axvline(0,color='k')
        ax.label_outer()
        ax.grid()
        ax.text(-0.06, 1.05, letters[letter_index]+')', transform=ax.transAxes, 
                    size=16, weight='bold')
        letter_index += 1

        #ax.coastlines()
        # xr.plot.contour(ds['WS'],ax=ax)
        
#fig.suptitle('100 most extreme cyclones',fontsize=16)
fig.tight_layout()
if compare_med==True:
    fig.savefig(path+'figs/upperatm_evolution_compare_med{}{}.png'.format(rank,maxstring),dpi=400)
else:
    fig.savefig(path+'figs/upperatm_evolution_compare_number{}.png'.format(rank,maxstring),dpi=400)


#%%
regions = ['cmed','emed']
events = ['prec']
path = '/storage/climatestor/PleioCEP/doensen/data/cyclone_intensity/'

for region in regions:
    for event in events:
        df = pd.read_csv(path+'{}/{}_check_latlon.csv'.format(region,event),index_col=0)
        df.columns = np.arange(-8,9)
        counts = df.apply(pd.value_counts)
        counts = counts.fillna(0)
        obj = FigureTemplate(nrows=5,ncols=4,figsize=(20,10),cartopy_projection=True)
        obj.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05,wspace=0.05,hspace=0.15)
        fig = obj.fig; axz = obj.axz
        for i,time in enumerate(counts):
            ax = axz.flat[i]
            ax.set_title(time)
            column = counts[time].iloc[:-1]
            column.index = pd.MultiIndex.from_tuples(column.index.map(eval), names=['latitude', 'longitude'])
            lats,lons = zip(*column.index)
            lats = list(lats); lons = list(lons)
            lon_grid,lat_grid = np.meshgrid(np.arange(min(lons),max(lons)+1,2.5),np.unique(lats)[::-1])
            count_grid = np.zeros_like(lat_grid)
            for lat,lon,count in zip(lats,lons,column.values):
                lat_idx = np.where(lat_grid[:,0]==lat)[0][0]
                lon_idx = np.where(lon_grid[0]==lon)[0][0]
                count_grid[lat_idx,lon_idx] = count
            count_grid[count_grid==0]=np.nan
            ax.pcolormesh(lon_grid[0],lat_grid[:,0],count_grid,vmin=0,vmax=20)
            ax.coastlines('110m')
        fig.suptitle('Region: {} | Event: {}'.format(region,event))
        fig.savefig(path+'figs/{}_{}{}latlon.png'.format(region,event,maxstring),dpi=400)





# %%
