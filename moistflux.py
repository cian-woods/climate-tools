from ReanalysisDataServer import DataServer as reDataServer
from netCDF4 import Dataset

import numpy as np
import os

def makeFile(FileName,F,hours,lons,lats,key,longname):

	print 'Creating %s ...' % (FileName)

	ntime,nlat,nlon = F.shape

	# Create file
	File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
	# Define some global attribs
	File.Conventions='COARDS'
	# Time is record dimension
	File.createDimension('time',ntime)
	var = File.createVariable('time','d',('time',))
	var.long_name = 'time'
	var.units = 'hours since 1900-01-01 00:00:0.0'
	var[:] = hours

	# axes
	File.createDimension(latnm,nlat)
	var = File.createVariable(latnm,'f',(latnm,))
	var.long_name = 'latitude'
	var.units = 'degrees north'
	var[:] = lats.astype('f')
	File.createDimension(lonnm,nlon)
	var = File.createVariable(lonnm,'f',(lonnm,))
	var.long_name = 'longitude'
	var.units = 'degrees east'
	var[:] = lons.astype('f')

	# create variables
	var = File.createVariable(key,'f',('time',latnm,lonnm))
	var.long_name = longname
	var.units = 'kg s**-1 m**-1'
	var[:] = F

	File.close()


du = reDataServer(Field='U',LevType='plev',LatRange=(50,90))
dv = reDataServer(Field='V',LevType='plev',LatRange=(50,90))
dq = reDataServer(Field='q',LevType='plev',LatRange=(50,90))

levs = du.lev
lats = du.lat
lons = du.lon

# Mass of girdbox on each level
dP_ = []
dP  = np.diff(levs)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665

latnm,lonnm  = 'latitude','longitude'
years,months = range(1979,2018+1,1),range(1,12+1,1)

for year in years:
	for month in months:
		FileName_vq = '/mnt/climstorage/cian/moistflux/vq/vq_%s_%02d.nc' % (year,month)
		FileName_uq = '/mnt/climstorage/cian/moistflux/uq/uq_%s_%02d.nc' % (year,month)
		if (os.path.isfile(FileName_vq)==False) or (os.path.isfile(FileName_uq)==False):
			uq,vq,q = du.getDataSnaps(Year=year,Month=month),dv.getDataSnaps(Year=year,Month=month),dq.getDataSnaps(Year=year,Month=month)
			uq,vq   = uq*q,vq*q
			# Vertically integrate
			uq,vq   = (uq*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1),(vq*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)
			hours    = du.getHourList(Year=year,Month=month) 
			# Make NetCDF file
			if not os.path.isfile(FileName_uq): makeFile(FileName_uq,uq,hours,lons,lats,'uq','Zonal moisture flux')
			if not os.path.isfile(FileName_vq): makeFile(FileName_vq,vq,hours,lons,lats,'vq','Meridional moisture flux')
