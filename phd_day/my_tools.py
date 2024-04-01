import pandas as pd
#from pylab import *
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

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