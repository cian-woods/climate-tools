from netCDF4 import *
from ReanalysisMonthlyMeanDataServer import *
from LambertProjector import *

import sys
import matplotlib.pyplot as pl
import numpy as np

key = str(sys.argv[1])

if key not in ['slhf','sshf','ssr','str']: sys.exit()

# NetCDF in file
FileName_in = '/mnt/climstorage/cian/%s.mon.mean.nc' % (key)
d_in   = Dataset(FileName_in,'r')
Var_bi = d_in.variables[key][:]
time   = d_in.variables['time'][:]

n0,n1,n2 = Var_bi.shape
Var_mo   = Var_bi.reshape((-1,2,n1,n2)).sum(axis=1)/(48*60*60.)
time     = time[0::2]
n0,n1,n2 = Var_mo.shape


FileName_out = '%s.new.nc' % (FileName_in[0:-3])
print 'Creating %s ...' % (FileName_out)

# Create file
File = Dataset(FileName_out,'w',format='NETCDF3_CLASSIC')
# Define some global attribs
File.Conventions='COARDS'
# Time is record dimension
File.createDimension('time',n0)
var = File.createVariable('time','d',('time',))
var.long_name = d_in.variables['time'].long_name
var.units     = d_in.variables['time'].units
var[:]        = time

# axes
File.createDimension('latitude',n1)
var           = File.createVariable('latitude','f',('latitude',))
var.long_name = d_in.variables['latitude'].long_name
var.units     = d_in.variables['latitude'].units
var[:]        = d_in.variables['latitude'][:].astype('f')
File.createDimension('longitude',n2)
var = File.createVariable('longitude','f',('longitude',))
var.long_name = d_in.variables['longitude'].long_name
var.units     = d_in.variables['longitude'].units
var[:]        = d_in.variables['longitude'][:].astype('f')

# create variables
var           = File.createVariable(key,'f',('time','latitude','longitude'))
var.long_name = d_in.variables[key].long_name
var.units     = 'W m**-2'
var[:]        = Var_mo

File.close()
d_in.close()


# Test climatology
varnames    = ['slhf','sshf','ssr','str']
servernames = ['slhf','sshf','fsns','fls']
Field       = servernames[varnames.index(key)]

proj = LambertProjector(boundinglat=70,resolution=80.)
d    = DataServer(Field=Field,LevType='surface_forecast')
s    = proj(np.array([d.getSeason(year,Season='DJF') for year in range(1980,2016+1,1)]).mean(axis=0),d.lon,d.lat)

cseq = np.arange(-80,80+10,10)
cf   = pl.contourf(proj.x,proj.y,s,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('%s [%s]' % (d.long_name,d.units))
proj.m.drawcoastlines()
pl.show()
