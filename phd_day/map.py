#%%
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
path='/storage/climatestor/PleioCEP/doensen/data/figs/Toulouse/'
# %%
fig,ax = plt.subplots(figsize=(15,5),subplot_kw={'projection':ccrs.PlateCarree()})
ax.set_extent([-30, 60, 20, 55], crs=ccrs.PlateCarree())
ax.coastlines(resolution='110m') # add map
ax.plot([-10, 40, 40, -10, -10], 
        [28, 28, 47, 47, 28],
         color='black', linewidth=2,zorder=3,
         #transform=ccrs.Geodetic(), #remove this line to get straight lines
         )
gl=ax.gridlines(draw_labels=True)
gl.xlabels_top = False
gl.ylabels_right = False
ax.text(15,37.5,'MED',va='center',ha='center',
        fontsize=30,color='darkred',weight='bold')

# ax.plot([-10, 2.5, 2.5, -10, -10], 
#         [30, 30, 47, 47, 30],
#          color='darkred', linewidth=2,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )

# ax.plot([2.5, 17.5, 17.5, 2.5, 2.5], 
#         [30, 30, 47, 47, 30],
#          color='darkred', linewidth=2,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )

# ax.text(10,38.5,'CMED',va='center',ha='center',
#         fontsize=16,color='darkred',weight='bold')

# ax.plot([17.5, 40, 40, 17.5, 17.5], 
#         [28, 28, 42, 42, 28],
#          color='darkred', linewidth=2,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )


# ax.text(28.75,35,'EMED',va='center',ha='center',
#         fontsize=16,color='darkred',weight='bold')

# ax.plot([-10, 40, 40, -10, -10], 
#         [47, 47, 70, 70, 47],
#          color='darkred', linewidth=1,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )
# ax.text(15,58.5,'EUR',va='center',ha='center',
#         fontsize=16,color='darkred',weight='bold')
# ax.plot([-70, -10, -10, -70, -70], 
#         [28, 28, 47, 47, 28],
#          color='darkred', linewidth=1,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )
# ax.text(-40,37.5,'SATL',va='center',ha='center',
#         fontsize=16,color='darkred',weight='bold')
# ax.plot([-70, -10, -10, -70, -70], 
#         [47, 47, 70, 70, 47],
#          color='darkred', linewidth=1,
#          #transform=ccrs.Geodetic(), #remove this line to get straight lines
#          )
# ax.text(-40,58.5,'NATL',va='center',ha='center',
#         fontsize=16,color='darkred',weight='bold')
fig.savefig(path+'map_med.png',dpi=300)
# %%
