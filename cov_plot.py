import pandas as pd
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt

path = "/storage/climatestor/PleioCEP/doensen/data/cyclone/"
file = "fort_34_3491_3500"

header= ['iic', 't0', 'age', 'ilon', 'ilat', 'cat', 'date' , 'cz' , 'di' ,
         'realc', 'gz', 'agetot' , 'idumi', 'zrad' , 'zdep'
         , 'slp', 'precmax','precmean']

df = pd.read_csv(path+file,delim_whitespace=True,names=header,skiprows=[0,1])

# %%
colorz = [(0,'#006600'),(0.125,'#33cc33'),
          (0.25,'#ffff00'),(0.375,'#ff6600'),(0.5,'#ff0000'),(0.625,'#ff6600'),(0.75,'#ffff00'),
          (0.875,'#33cc33'),(1,'#006600')]   
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('custom', colorz, N=1024)

def cov_plot(df):
    columns = df.columns
    fig,axz = plt.subplots(len(columns),len(columns),figsize=(18,13))
    for y in range(len(columns)):
        for x in range(len(columns)):
            ax=axz[x,y]
            if y==x:
                df[columns[x]].plot.hist(ax=ax,bins=30)
                ax.set_ylabel([])
            if y<x:
                ax.scatter(df[columns[y]],df[columns[x]],
                           marker='.',s=0.2)
    
            if x<y:
                corr = df[columns[[y,x]]].corr().values[1,0]
                ax.text(0.5 , 0.5 , '{:.3f}'.format(corr),
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes,
                fontsize=16)
                ax.set_facecolor(cmap((corr+1)/2))
            ax.set_xlabel(columns[y])
            ax.set_ylabel(columns[x])
            ax.label_outer()
    fig = fig.tight_layout()
    return fig
        
            
# %%

df_iic_mean = df[['iic','cz','gz','zrad','zdep','slp','precmax','precmean']]\
    .groupby('iic').mean()

fig = cov_plot(df_iic_mean)

# %%

df_iic_min = df[['iic','cz','slp']].groupby('iic').min()
df_iic_max = df[['iic','gz','zrad','zdep','precmax','precmean']].groupby('iic').max()
df_iic_ext = df_iic_min.join(df_iic_max)\
    .reindex(columns=['cz','gz','zrad','zdep','slp','precmax','precmean'])

fig = cov_plot(df_iic_ext)

# %%
df_gz_75p = df_iic_ext.where(df_iic_ext['gz']>162.8).dropna()

fig = cov_plot(df_gz_75p)
# %%

df_slp = df_iic_ext.where((df_iic_ext['slp']>1005)&(df_iic_ext['slp']<1010)).dropna()

fig = cov_plot(df_slp)
    

    

    

