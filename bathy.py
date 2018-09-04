from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from LambertProjector import *

import matplotlib.pyplot as pl
import numpy as np

proj = LambertProjector(boundinglat=70,resolution=80.)
d    = Dataset('/mnt/climstorage/cian/Bathymetry_etopo2-cian.nc','r')

lon = d.variables['lon'][:]
lat = d.variables['lat'][:]
b   = d.variables['Bathymetry'][:].squeeze()
b   = proj(b,lon,lat)

cseq = np.arange(-600+5,5)
cf   = pl.contourf(proj.x,proj.y,b,cseq,cmap=pl.cm.coolwarm,extend='both')
cbar = pl.colorbar(cf)
proj.m.drawcoastlines()
pl.show()

