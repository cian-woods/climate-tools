from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from scipy import interpolate

import numpy as np
import sys

def wave1(v):
	# Takes in (time)x(lat)x(lon) field and returns the wave1 pattern
	# Divide by n becasue discrete transform and extract wave1 (i0 index after shift)
	ftn = v.shape[-1]
	i0  = int(ftn/2) + 1
	Fk  = np.fft.fft(v,axis=2)/self.ftn
	Fk  = np.fft.fftshift(Fk,axes=(2,))[:,:,self.i0]
	# Calcutae amplitide and phase of wave1
	amp = 2*np.abs(Fk)
	ang = np.arctan2(Fk.imag,Fk.real)
	return amp,ang

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

def makeNetCDF(field,times,lats,lons,FileName,case):
        print 'Creating %s ...' % (FileName)

	if case == 'vq': orr = 'meridional'
	if case == 'uq': orr = 'zonal'
	ntime,nlon,nlat = len(times),len(lons),len(lats)

        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'Time'
        var.units     = 'hours since 1800-01-01 00:00:0.0'
        var[:]        = times
        # Horizontal axes
        File.createDimension('lat',nlat)
        var           = File.createVariable('lat','f',('lat',))
        var.long_name = 'Latitude'
        var.units     = 'degrees_north'
        var[:]        = lats.astype('f')
        File.createDimension('lon',nlon)
        var           = File.createVariable('lon','f',('lon',))
        var.long_name = 'Longitude'
        var.units     = 'degrees_east'
        var[:]        = lons.astype('f')
        # Create Variables
        var           = File.createVariable(case,'f',('time','lat','lon'))
        var.long_name = '4xDaily vertically integrated %s moisture flux' % (orr)
        var.units     = 'kg s**-1 m**-1'
        var[:]        = field
	# Close file
        File.close()

case = str(sys.argv[1])
# DataServers
if case == 'vq': vds = DataServer(Field='vplev',LevType='plev',LevRange=(0,1000),LatRange=(20,90),Source='NCEP1')
if case == 'uq': vds = DataServer(Field='uplev',LevType='plev',LevRange=(0,1000),LatRange=(20,90),Source='NCEP1')
qds = DataServer(Field='shum',LevType='plev',LevRange=(0,1000),LatRange=(20,90),Source='NCEP1')

# Mass of girdbox on each level
plevs = np.linspace(300,1000,100)
dP    = []
dp    = np.diff(plevs)/2.
dP.append(dp[0])
for i in range(len(dp)-1):
        dP.append(dp[i]+dp[i+1])
dP.append(dp[-1])
dP = np.array(dP)
dM = 100.*dP/9.80665                             # Mass per unit area between pressure levels of data (kg/m^2)

years,months = range(1950,2016+1,1),range(1,12+1,1)
Dir          = '/mnt/climstorage/cian/NCAR1/%s' % (case)
# Make data
for year in years:
	for month in months:
		fname = '%s/%s_%s_%02d.nc' % (Dir,case,year,month)
		if not os.path.isfile(fname):
			v,q   = vds.getMonth(year,month),qds.getMonth(year,month)	
			hours = [vds.getHours(*date) for date in vds.getDateList(year,month)]
			v     = interpolateND(v,vds.lev,plevs,axis=1)
			q     = interpolateND(q,qds.lev,plevs,axis=1)
			vq    = ((v*q)*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)	# Take product and vertially integrate
			# Reorder axes and data
			vq,lats = vq[:,::-1,:],vds.lat[::-1]
			makeNetCDF(vq,hours,lats,vds.lon,fname,case)
