import pandas as pd
#from pylab import *
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

class FigureTemplate:
    def __init__(self, nrows=1, ncols=1,figsize=(16,9),cartopy_projection=False):
        self.cartopy_projection = cartopy_projection
        if cartopy_projection:
            self.fig, self.axz = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                                              subplot_kw={'projection': ccrs.PlateCarree()})
        else:
            self.fig, self.axz = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    def subplots_adjust(self,left=0.05,right=0.95,top=0.95,bottom=0.05,wspace=0.05,hspace=0.05):
        self.fig.subplots_adjust(left=left,right=right,top=top,bottom=bottom,wspace=wspace,hspace=hspace)
            

def ilon_to_lon(ilon,nlon=144):
    lon = (ilon-1)*(360/nlon)
    return lon



def ilat_to_lat(ilat,nlat=96):
    constant = 0.9473684210526301
    lat = -180*((ilat-1.5)/(nlat-1)-0.5)
    lat = lat-constant
    return lat
    
class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""

    name = 'nlcmap'

    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels/ self.levels.max()
        self.levmax = self.levels.max()
        self.levmin = self.levels.min()
        self._y = np.linspace(self.levmin, self.levmax, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi/self.levmax, alpha)