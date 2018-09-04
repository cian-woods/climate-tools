from cmipDataServer import DataServer as cmipDataServer
from scipy import interpolate
from netCDF4 import Dataset as Dataset

import numpy as np
import sys,os

def fill(data,invalid=None):
                if invalid is None: invalid = np.isnan(data)
                ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
                return data[tuple(ind)]

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

def makeNetCDF(field,times,time_units,time_calendar,lats,lons,FileName):
        print 'Creating %s ...' % (FileName)

	ntime,nlon,nlat = len(times),len(lons),len(lats)

        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'time'
        var.units     = time_units
	var.calendar  = time_calendar
        var[:]        = times
        # Horizontal axes
        File.createDimension('lat',nlat)
        var           = File.createVariable('lat','f',('lat',))
        var.long_name = 'latitude'
        var.units     = 'degrees_north'
        var[:]        = lats.astype('f')
        File.createDimension('lon',nlon)
        var           = File.createVariable('lon','f',('lon',))
        var.long_name = 'longitude'
        var.units     = 'degrees_east'
        var[:]        = lons.astype('f')
        # Create Variables
        var           = File.createVariable('prw','f',('time','lat','lon'))
        var.long_name = 'Precipitable water'
        var.units     = 'kg m**-2'
        var[:]        = field
	# Close file
        File.close()

# DataServers
Source = str(sys.argv[1])
Exp    = str(sys.argv[2])
qds    = cmipDataServer(Field='hus',LevType='plev',LevRange=(0,1000),LatRange=(-90,90),Source=Source,ExpType=Exp,DataFreq='day')
lats   = qds.lat
lons   = qds.lon

# Mass of girdbox on each level
plevs = np.linspace(30,1000,100)
dP    = []
dp    = np.diff(plevs)/2.
dP.append(dp[0])
for i in range(len(dp)-1):
        dP.append(dp[i]+dp[i+1])
dP.append(dp[-1])
dP = np.array(dP)
dM = 100.*dP/9.80665                             # Mass per unit area between pressure levels of data (kg/m^2)

years,months = range(2074,2100+1,1),[12,1,2]
#years,months = range(1981,2005+1,1),[12,1,2]
Dir          = '/mnt/climstorage/cian/%s/%s/day/surface/prw' % (Exp,Source)
os.system('mkdir -p %s' % (Dir))
# Make data
for year in years:
	for month in months:
		fname = '%s/prw_%s_%02d.nc' % (Dir,year,month)
		if not os.path.isfile(fname):
			q     = qds.getMonth(year,month)	
			if (type(q)==np.ma.core.MaskedArray) and (q.mask.shape == q.data.shape):
				q  = fill(q.data,invalid=q.mask)
			hours = [qds.getDays(*date) for date in qds.getDateList(year,month)]
			q     = interpolateND(q,qds.lev,plevs,axis=1)
			q    = (q*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)	# Take product and vertially integrate
			# Reorder axes and data
			if qds.InvertLatAxis: q,lats = q[:,::-1,:],qds.lat[::-1]
			if qds.InvertLonAxis: q,lons = q[:,:,::-1],qds.lon[::-1]
			makeNetCDF(q,hours,qds.time_units,qds.calendar,lats,lons,fname)
