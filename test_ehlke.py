from netCDF4 import *
from LambertProjector import *

import matplotlib.pyplot as pl
import numpy as np



d = Dataset('/mnt/climstorage/obs/ERAInt/monthly/surface_analysis/tcw.mon.mean.nc','r')
s = d.variables['tcw'][:]
lat = d.variables['latitude'][:]
lon = d.variables['longitude'][:]

latx = np.argmin((lat--0)**2)

sJ,sA,sS = s[::6,latx:,:].mean(axis=0),s[::7,latx:,:].mean(axis=0),s[::8,latx:,:].mean(axis=0)
sJAS     = (sJ+sA+sS)/3

pl.contourf(lon,lat[latx:],sJAS)
pl.colorbar()
pl.show()
