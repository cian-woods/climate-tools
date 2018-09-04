#from ReanalysisMonthlyMeanDataServer import *
from ReanalysisDataServer import *
from toPick import *
import os

def meanField(	Field 	   =  'strd',
		YearRange  =   (1979,2012),
		Months	   =   arange(12)+1,
		Source	   =  'ERAInt',
		LevType	   =  'surface_forecast',
		LevRange   =   (0,1000),
		LatRange   =   (-90,90),
		LonRange   =   (0,360)		):

	d = DataServer(Field=Field,LevType=LevType,Source=Source,\
			LevRange=LevRange,LatRange=LatRange,LonRange=LonRange)

	var_name  = d.long_name
	var_units = d.units

	lev = d.lev
	if LevType == 'plev':
		lev_longname = 'pressure_level'
		lev_units    = 'hPa'
	elif LevType == 'modlev':
		lev_longname = 'model_level'
		lev_units    = 'int'
	else:
		lev = None

	S,R = [],[]
	for year in range(YearRange[0],YearRange[1]+1,1):

		for month in Months:

			print year,month

			try:
				data = d.getMonth(year,month)
				hour = d.getHours(year,month,1,0)

				S.append(mean(data,axis=0))
				R.append(hour)
			except:
				pass

	S = array(S)
	R = array(R)
	if Field == 'ci':
		S[where(S<0)] = 0

	key  = d.FieldName
	Var  = [S,var_name,var_units]
	Time = [R,'time',d.time_units]
	Lon  = [d.lon,'longitude','degrees_east']
	Lat  = [d.lat,'latitude','degrees_north']
	if lev != None:
		Lev = [lev,lev_longname,lev_units]
	else:
		Lev = []

	return key,Var,Time,Lon,Lat,Lev,d.FieldNames['lat'],d.FieldNames['lon']

if __name__ == "__main__":

#	import matplotlib.pyplot as pl

	# ------------- Input ------------- #

	Field 	    =   'uq'
	YearRange   =   (1979,2018)
	Months	    =   arange(12)+1
	Source	    =   'ERAInt'
	LevType	    =   'surface_analysis'
	LevRange    =   (0,1000)
	LatRange    =   (50,90)
	LonRange    =   (0,360)

	# -------------------------------- #


	key,Var,Time,Lon,Lat,Lev,latnm,lonnm = meanField(Field=Field,YearRange=YearRange,Months=Months,Source=Source,LevType=LevType,\
					                 LevRange=LevRange,LatRange=LatRange,LonRange=LonRange)

	ntime,nlon,nlat = len(Time[0]),len(Lon[0]),len(Lat[0])

	print shape(Var[0])


	FileName = '/mnt/climstorage/cian/%s.mon.mean.nc' % (key)
	print 'Creating %s ...' % (FileName)

	# Create file
	File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
	# Define some global attribs
	File.Conventions='COARDS'
        # Time is record dimension
   	File.createDimension('time',ntime)
    	var = File.createVariable('time','d',('time',))
    	var.long_name = Time[1]
    	var.units = Time[2]
	var[:] = Time[0]

        # Create z dimension if it exists
        if Lev != []:
                File.createDimension('level',len(Lev[0]))
                var = File.createVariable('level','f',('level',))
                var.long_name = Lev[1]
                var.units = Lev[2]
                var[:] = Lev[0].astype('f')

    	# axes
    	File.createDimension(latnm,nlat)
    	var = File.createVariable(latnm,'f',(latnm,))
    	var.long_name = Lat[1]
    	var.units = Lat[2]
    	var[:] = Lat[0].astype('f')
    	File.createDimension(lonnm,nlon)
    	var = File.createVariable(lonnm,'f',(lonnm,))
    	var.long_name = Lon[1]
    	var.units = Lon[2]
    	var[:] = Lon[0].astype('f')

    	# create variables
	if Lev != []:
   	 	var = File.createVariable(key,'f',('time','level',latnm,lonnm))
	else:
		var = File.createVariable(key,'f',('time',latnm,lonnm))
    	var.long_name = Var[1]
    	var.units = Var[2]
	var[:] = Var[0]
	
	File.close()

