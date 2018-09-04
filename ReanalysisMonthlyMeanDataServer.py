#!/usr/bin/env python

import os,sys
from numpy import *
from netCDF4 import Dataset as NetCDFFile
from progress import ProgressMeter
import time,datetime,calendar

HomeDir = '/mnt/climstorage'


# list of all data holdings in ERA
FieldNamesERA = {}
FieldNamesERA['time']  	=  ['time','coord']
FieldNamesERA['lev']   	=  ['level','coord']
FieldNamesERA['lat']   	=  ['latitude','coord']
FieldNamesERA['lon']   	=  ['longitude','coord']
FieldNamesERA['U']     	=  ['u','plev']
FieldNamesERA['V']     	=  ['v','plev']
FieldNamesERA['W']     	=  ['w','plev']
FieldNamesERA['Z']     	=  ['z','plev']
FieldNamesERA['T']     	=  ['t','plev']
FieldNamesERA['q']     	=  ['q','plev']
FieldNamesERA['PV']     =  ['pv','plev']
FieldNamesERA['PVth']   =  ['pv','thetalev']
FieldNamesERA['lnsp']  	=  ['lnsp','modlev']
FieldNamesERA['slp']   	=  ['msl','surface_analysis']
FieldNamesERA['skt']   	=  ['skt','surface_analysis']
FieldNamesERA['sst']   	=  ['sst','surface_analysis']
FieldNamesERA['T2']    	=  ['t2m','surface_analysis']
FieldNamesERA['Ts']    	=  ['skt','surface_analysis']
FieldNamesERA['U10']   	=  ['u10','surface_analysis']
FieldNamesERA['V10']   	=  ['v10','surface_analysis']
FieldNamesERA['pw']    	=  ['tcw','surface_analysis']
FieldNamesERA['pwv']   	=  ['tcwv','surface_analysis']
FieldNamesERA['cld']   	=  ['cld','surface_analysis']
FieldNamesERA['ci']    	=  ['ci','surface_analysis']
FieldNamesERA['ie']     =  ['ie','surface_analysis']
FieldNamesERA['precc'] 	=  ['cp','surface_forecast']
FieldNamesERA['precl'] 	=  ['lsp','surface_forecast']
FieldNamesERA['flds']  	=  ['strd','surface_forecast']
FieldNamesERA['fsns']   =  ['ssr','surface_forecast']
FieldNamesERA['fls']  	=  ['str','surface_forecast']
FieldNamesERA['sshf']	=  ['sshf','surface_forecast']
FieldNamesERA['th2PV']	=  ['pt','PVlev']
#FieldNamesERA['ishf']	=  ['ishf','surface_analysis']
FieldNamesERA['slhf']   =  ['slhf','surface_forecast']
FieldNamesERA['cp']     =  ['cp','surface_forecast']
FieldNamesERA['ttr']    =  ['ttr','surface_forecast']
FieldNamesERA['vq']     =  ['vq','surface_analysis']
FieldNamesERA['uq']     =  ['uq','surface_analysis']
# List of all data holdings in CAM
FieldNamesCAM         = {}
FieldNamesCAM['time'] = ['time','coord']
FieldNamesCAM['lat']  = ['lat','coord']
FieldNamesCAM['lon']  = ['lon','coord']
FieldNamesCAM['lev']  = ['lev','coord']
FieldNamesCAM['Ts']   = ['TS','surface']
FieldNamesCAM['U']    = ['U','plev']
FieldNamesCAM['V']    = ['V','plev']
FieldNamesCAM['q']    = ['Q','plev']
FieldNamesCAM['T']    = ['T','plev']
FieldNamesCAM['W']    = ['OMEGA','plev']
FieldNamesCAM['slp']  = ['PSL','surface']
FieldNamesCAM['pw']   = ['TMQ','surface']
FieldNamesCAM['ttr']  = ['FLNT','surface']
FieldNamesCAM['tsr']  = ['FSNT','surface']
FieldNamesCAM['sshf'] = ['SHFLX','surface']
FieldNamesCAM['slhf'] = ['LHFLX','surface']
FieldNamesCAM['fls'] = ['FLNS','surface']
FieldNamesCAM['Z'] = ['Z3','plev']
# List of all data holdings in NCAR20C_V2c
FieldNamesNCAR20C         = {}
FieldNamesNCAR20C['time'] = ['time','coord']
FieldNamesNCAR20C['lat']  = ['lat','coord']
FieldNamesNCAR20C['lon']  = ['lon','coord']
FieldNamesNCAR20C['lev']  = ['lev','coord']
FieldNamesNCAR20C['T2']   = ['air','surface']
FieldNamesNCAR20C['slp']  = ['prmsl','surface']
FieldNamesNCAR20C['ci']   = ['icec','surface']
FieldNamesNCAR20C['U10']  = ['uwnd','surface']
FieldNamesNCAR20C['V10']  = ['vwnd','surface']

class DataServer:
    def __init__(self,
                 LevType  = None,
                 Field    = 'U',
                 Source	  ='ERAInt',
                 LevRange = (0.,1000.), 
                 LatRange = (-90.,90.), 
                 LonRange = (0.,360.)   ):

	self.Field  = Field
	self.dayref = 5				# dayref is for (year,month,dayref,0). Used for finding month data.
	self.leap   = True			# Some data uses different time of month to define month eg (1998,1,1,0) vs (1998,1,31,18)
        if Source[0:3] == 'ERA':
            self.FieldNames = FieldNamesERA
            if Source == 'ERA40':
                self.FieldNames['lat']    = 'lat'
                self.FieldNames['lon']    = 'lon'
	    self.LevType = self.FieldNames[Field][1]
            # set up field name and directory
            self.FieldName = self.FieldNames[Field][0]
	    # Access any pre-specified netcdf files not provided by Rodge (typically monthly climatologies)
	    if Field == 'flds':
		Dir = '/mnt/climstorage/cian/strd.mon.mean.nc'
            elif Field == 'ci':
                Dir = '/mnt/climstorage/cian/ci.mon.mean.nc'
            elif Field == 'U':
                Dir = '/mnt/climstorage/cian/u.mon.mean.nc'
            elif Field == 'V':
                Dir = '/mnt/climstorage/cian/v.mon.mean.nc'
            elif Field == 'Z':
                Dir = '/mnt/climstorage/cian/z.mon.mean.nc'
            elif Field == 'U10':
                Dir = '/mnt/climstorage/cian/u10.mon.mean.nc'
            elif Field == 'V10':
                Dir = '/mnt/climstorage/cian/v10.mon.mean.nc'
            elif Field == 'T':
                Dir = '/mnt/climstorage/cian/t.mon.mean.nc'
            elif Field == 'slp':
                Dir = '/mnt/climstorage/cian/msl.mon.mean.nc'
            elif Field == 'pw':
                Dir = '/mnt/climstorage/cian/tcw.mon.mean.nc'
            elif Field == 'sst':
                Dir = '/mnt/climstorage/cian/sst.mon.mean.nc'
            elif Field == 'slhf':
                Dir = '/mnt/climstorage/cian/slhf.mon.mean.new.nc'
            elif Field == 'sshf':
                Dir = '/mnt/climstorage/cian/sshf.mon.mean.new.nc'
	    elif Field == 'fls':
                Dir = '/mnt/climstorage/cian/str.mon.mean.new.nc'
            elif Field == 'fsns':
                Dir = '/mnt/climstorage/cian/ssr.mon.mean.new.nc'
	    elif LevType == 'modlev':
		Dir = '/qnap/cian/modlev/%s.mon.mean.nc' % (self.FieldName)
	    elif Field == 'ttr':
		Dir = '/mnt/climstorage/cian/%s.mon.mean.nc' % (self.FieldName)
            elif Field == 'T2':
                Dir = '/mnt/climstorage/cian/%s.mon.mean.nc' % (self.FieldName)
	    elif Field == 'vq':
                Dir = '/mnt/climstorage/cian/%s.mon.mean.nc' % (self.FieldName)
            elif Field == 'uq':
                Dir = '/mnt/climstorage/cian/%s.mon.mean.nc' % (self.FieldName)
	    else:
		Dir = HomeDir+'/obs/%s/monthly/%s/%s.mon.mean.nc' % (Source,self.LevType,self.FieldName)
	# CAM output data
	elif Source[0:3] == 'CAM':
	    self.dayref     = 25
	    self.leap       = False
            self.FieldNames = FieldNamesCAM
            self.LevType    = self.FieldNames[Field][1]
            # set up field name and directory
            self.FieldName  = self.FieldNames[Field][0]
	    Dir             = '%s/cian/EXPS/Expt_%s/h0/CAM_h0.nc' % (HomeDir,Source[3:])
        # NCAR20C_V2c output data
        elif Source[0:7] == 'NCAR20C':
            self.dayref     = 15
            self.leap       = True
            self.FieldNames = FieldNamesNCAR20C
            self.LevType    = self.FieldNames[Field][1]
            # set directory
            self.FieldName  = self.FieldNames[Field][0]
            Dir             = '%s/cian/%s/monthly/%s/%s.mon.mean.nc' % (HomeDir,Source,self.LevType,self.FieldName)
	self.Dir = Dir
	# Scale factors
	if Field == 'slhf':
	    self.sf = 1#2260e03
	else:
	    self.sf = 1.
	# NetCDF file
	print self.Dir
	self.File = NetCDFFile(Dir,'r')
	# Initialise data
        self.units = self.File.variables[self.FieldName].units
        self.long_name = self.File.variables[self.FieldName].long_name
        # base year for time computation
        self.Year0 = self._getBaseYear(self.File)
        # time axis
        self.time = self.File.variables['time'][:]
        self.FirstYear = self.getDate(self.time[0])[0]
        self.LastYear  = self.getDate(self.time[-1])[0]
        self.Years = range(self.FirstYear,self.LastYear+1)
        self.nyears = len(self.Years)
        # Print out info
        print '-------------'
        print 'Instantiated Reanalysis Monthly Mean Data Server'
#        print 'Data from: %s' % (HomeDir+'/obs/%s/monthly/%s/%s.mon.mean.nc'\
#                              %(Source,self.LevType,self.FieldName))
	print 'Data from: %s' % (Dir)
        print 'Field: %s' % Field
        print 'Year range: %s-%s (%s years)' %\
              (self.FirstYear,self.LastYear,self.nyears)


        # Initialize coord axes
        self.InvertLatAxis,self.InvertLonAxis = False,False
        self.j0,self.j1,self.i0,self.i1       = None,None,None,None
        try:
        	latnames = ['latitude','lat']
        	latname  = list(set(latnames) & set(self.File.variables.keys()))[0]
        	lat = array(self.File.variables[latname][:])*1.
        	(self.lat, self.j0, self.j1, self.InvertLatAxis) = self._setupAxis(lat,LatRange)
        	self.nlat = len(self.lat)
        except:
                pass
        try:
        	lonnames = ['longitude','lon']
        	lonname  = list(set(lonnames) & set(self.File.variables.keys()))[0]
        	lon = array(self.File.variables[lonname][:])*1.
        	(self.lon, self.i0, self.i1, self.InvertLonAxis) = self._setupLonAxis(lon,LonRange)
        	self.nlon = len(self.lon)
        except:
                pass
        # Level axis if it exists
        try:
        	levnames = ['lev','levelist','level']
        	levname  = list(set(levnames) & set(self.File.variables.keys()))[0]
        	lev = array(self.File.variables[levname][:])*1.
        	(self.lev, self.k0, self.k1, self.InvertLevAxis) = self._setupAxis(lev,LevRange)
        	self.nlev = len(self.lev)
        except:
                self.lev = None

	# Land Sea mask from ERAInt
        lsm_file      = NetCDFFile('/home/rca/obs/ERAInt/ERAInterimLandSeaMask_0.75x0.75.nc','r')
        self.lsm_lat  = lsm_file.variables['latitude'][:]
        self.lsm_lon  = lsm_file.variables['longitude'][:]
        self.lsm      = lsm_file.variables['lsm'][:].squeeze()
	if self.lsm_lat[1]<self.lsm_lat[0]: self.lsm_lat,self.lsm = self.lsm_lat[::-1],self.lsm[::-1,:]
	if self.lsm_lon[1]<self.lsm_lon[0]: self.lsm_lon,self.lsm = self.lsm_lon[::-1],self.lsm[:,::-1]
        #if self.InvertLatAxis: self.lsm = self.lsm[::-1,:]
        #if self.InvertLonAxis: self.lsm = self.lsm[:,::-1]
        #self.lsm = self.lsm[self.j0:self.j1, self.i0:self.i1]

    def _setupAxis(self,axis,range):
        if axis[0] > axis[-1]:
            axis = axis[::-1]
            invert = True
        else: invert = False
        i0  = argmin(abs(axis-range[0]))
        i1  = argmin(abs(axis-range[1]))+1
        axis = axis[i0:i1]
        return axis,i0,i1,invert

    def _setupLonAxis(self,axis,range):
        if axis[0] > axis[-1]:
            axis = axis[::-1]
            invert = True
        else: invert = False
        if range[0] < 0:
            axis = where(axis>180.,axis-360.,axis)
            self.shift = len(axis) - argmin(abs(axis-range[0]))
            axis = roll(axis,self.shift)
            i0 = 0
            i1  = argmin(abs(axis-range[1]))+1
            self.ShiftLon = True
        else:
            i0  = argmin(abs(axis-range[0]))
            i1  = argmin(abs(axis-range[1]))+1
            self.ShiftLon = False
        axis = axis[i0:i1]
        return axis,i0,i1,invert

    """
    def getDate(self,Hours,nd=1):
        # Given hours elapsed since midnight on 1 January Year0,
        # returns date (year, month, day, hour) 
        if self.leap == False:
                Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days=nd*Hours/24.)
                m      = 1
                S      = 0
                while m-S > 0:
                        if Date.month > 2:
                                m = calendar.leapdays(self.Year0,Date.year+1)
                        elif Date.month < 2:
                                m = calendar.leapdays(self.Year0,Date.year)
                        elif (Date.month == 2) and (Date.day < 29):
                                m = calendar.leapdays(self.Year0,Date.year)
                        else:
                                m = calendar.leapdays(self.Year0,Date.year+1)
                        if m-S == 0:
                                break
                        else:
                                Date  = Date + datetime.timedelta(days=float(m-S))
                                S = m
        elif self.leap == True:
        	Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(Hours/24.)
        return Date.year,Date.month,Date.day,Date.hour
    """

    def leapDaysBetweenDates(self,date1):

        if (date1[1]==2 and date1[2]<=28) or (date1[1]==1):
                year1,dn1 = date1[0],0
        elif (date1[1]==2) and (28 < date1[2]):
                year1,dn1 = date1[0],(date1[3]+6)/24.
        else:
                year1,dn1 = date1[0]+1,0

        n = calendar.leapdays(self.Year0,year1) + dn1
        return n

    def getDate(self,Hours,nd=1,Year0=None):
        # Given hours elapsed since midnight on 1 January Year0,
        # returns date (year, month, day, hour)
	if Year0 == None: Year0 = self.Year0
        if self.leap == False:
                Date = datetime.datetime(Year0,1,1,0) + datetime.timedelta(days=nd*Hours/24.)
                date = (Date.year,Date.month,Date.day,Date.hour)
                n    = self.leapDaysBetweenDates(date)
                N    = 0
                while n-N > 0:
                        N    = n
                        Date = datetime.datetime(Year0,1,1,0) + datetime.timedelta(days= n + nd*Hours/24.)
                        date = (Date.year,Date.month,Date.day,Date.hour)
                        n    = self.leapDaysBetweenDates(date)
        elif self.leap == True:
                Date = datetime.datetime(Year0,1,1,0) + datetime.timedelta(nd*Hours/24.)
        return Date.year,Date.month,Date.day,Date.hour

    def getDays(self,Year,Month,Day,Hour):
                hours = self.getHours(Year,Month,Day,Hour)
                days  = hours/24.
                return days

    def getHours(self,Year,Month,Day,Hour):
        # Given date (year, month, day, hour),
        # returns hours elapsed since midnight on 1 January Year0
        if self.leap == False:
                Days = datetime.datetime(int(Year),int(Month),int(Day),int(Hour)) \
                       - datetime.datetime(self.Year0,1,1,0)
                if Month > 2:
                        m = calendar.leapdays(self.Year0,Year+1)
                elif Month < 2:
                        m = calendar.leapdays(self.Year0,Year)
                elif (Month == 2) and (Day < 29):
                        m = calendar.leapdays(self.Year0,Year)
                else:
                        m = calendar.leapdays(self.Year0,Year+1)
                Hours = Days.days*24 + Days.seconds/3600 - m*24
        elif self.leap == True:
        	Days = datetime.datetime(Year,Month,Day,Hour) \
        	       - datetime.datetime(self.Year0,1,1,0)
        	Hours = Days.days*24 + Days.seconds/3600
        return Hours

    def getMonth(self, Year, Month):
        # Return 1 month of data.
        assert self.FirstYear <= Year <= self.LastYear, \
               'Year %s not in dataset !!' % Year
        h = self.getHours(Year,Month,self.dayref,0)
        i = argmin(abs(self.time-h))
        ## try:
        ##     scale = self.Field.scale_factor
        ##     offset = self.Field.add_offset
        ##     x = self.Field[i]*scale+offset
        ## except:
        f = self.File.variables[self.FieldName][i]

        # swap axes if necessary and select slab
        if len(f.shape) == 2:
            if self.InvertLatAxis: f = f[::-1,:]
            if self.InvertLonAxis: f = f[:,::-1]
            if self.ShiftLon: f = roll(f,self.shift,axis=-1)            
            f = f[self.j0:self.j1, self.i0:self.i1]
        if len(f.shape) == 3:
            if self.InvertLatAxis: f = f[:,::-1,:]
            if self.InvertLonAxis: f = f[:,:,::-1]
            try:
                if self.InvertLevAxis: f = f[::-1,:,:]
            except:
                pass
            if self.ShiftLon: f = roll(f,self.shift,axis=-1)            
            f = f[self.k0:self.k1, self.j0:self.j1, self.i0:self.i1]

        return self.sf*f

    def getMonthlyClimatology(self,Month):
        try:
            f = ma.masked_all((self.nyears,self.nlev,self.nlat,self.nlon))
        except:
            f = ma.masked_all((self.nyears,self.nlat,self.nlon))            
        for i in range(self.nyears):
            f[i] = self.getMonth(self.Years[i],Month)
        return f.mean(axis=0)

    def getSeason(self, Year, Season, mean=True):
        # Return 1 season of data.
        assert Season in ['Annual','Annual_wc','Annual_MF','MAMJJASO','MAMJJASO_l1','D','J','F','ON','DJ','FM','ONDJ','NDJ','DJF','JF','JFM','MAM','JJA','SON','NDJFM','M','S','OND','A','JFD','ONDJFM','DJFM','NDJF','ONDJF','AMJJAS','DJFMA'],\
               "Season must be one of 'DJF','MAM','JJA','SON'"
	if Season == 'Annual'      : Months = range(1,13)
	if Season == 'Annual_wc'   : Months = range(7,12+1,1) + range(1,6+1,1)
	if Season == 'Annual_MF'   : Months = range(3,12+1,1) + range(1,2+1,1)
	if Season == 'MAMJJASO'    : Months = [3,4,5,6,7,8,9,10]
	if Season == 'MAMJJASO_l1' : Months = [10,9,8,7,6,5,4,3]
	if Season ==  'ONDJ'       : Months = [10,11,12,1]
	if Season ==   'OND'       : Months = [10,11,12]
	if Season == 'NDJFM'       : Months = [11,12,1,2,3]
	if Season == 'ONDJFM'   : Months = [10,11,12,1,2,3]
	if Season == 'AMJJAS'   : Months = [4,5,6,7,8,9]
	if Season == 'ONDJF'    : Months = [10,11,12,1,2]
	if Season == 'DJFM'     : Months = [12,1,2,3]
	if Season == 'DJFMA'    : Months = [12,1,2,3,4]
	if Season == 'NDJF'     : Months = [11,12,1,2]
	if Season ==   'ON'     : Months = [10,11]
	if Season ==   'DJ'     : Months = [12,1]
	if Season ==   'FM'     : Months = [2,3]
	if Season ==  'NDJ'     : Months = [11,12,1]
        if Season ==  'DJF'     : Months = [12,1,2]
	if Season ==  'JFD'     : Months = [1,2,12]
	if Season ==   'JF'     : Months = [1,2]
	if Season ==  'JFM'     : Months = [1,2,3]
	if Season ==  'AMJ'     : Months = [4,5,6]
	if Season ==  'JAS'     : Months = [7,8,9]
	if Season ==  'OND'     : Months = [10,11,12]
        if Season ==  'MAM'     : Months = [3,4,5]
        if Season ==  'JJA'     : Months = [6,7,8]
        if Season ==  'SON'     : Months = [9,10,11]
	if Season ==    'D'     : Months = [12]
	if Season ==    'J'     : Months = [1]
	if Season ==    'F'     : Months = [2]
	if Season ==    'M'     : Months = [3]
	if Season ==    'S'     : Months = [9]
	if Season ==    'A'     : Months = [8]
        if len(self.File[self.FieldName][0].shape) >= 3:
            f = ma.masked_all((len(Months),self.nlev,self.nlat,self.nlon))
        else:
            f = ma.masked_all((len(Months),self.nlat,self.nlon))

	if Months[0] > Months[-1]:
		i = 0
		while Months[i]>Months[-1]:
			x    = self.getMonth(Year-1,Months[i])
			f[i] = x
			i    = i + 1
		while i <= len(Months)-1:
			x    = self.getMonth(Year,Months[i])
			f[i] = x
			i    = i + 1
	else:
        	for i in range(len(Months)):
		    x = self.getMonth(Year,Months[i])
                    f[i] = x
	if mean == True:
        	return f.mean(axis=0)
	else:
		return f

    def getSeasonalClimatology(self,Season):
        if Season == 'DJF':
            FirstYear = self.FirstYear+1
        else:
            FirstYear = self.FirstYear
        Years = range(FirstYear,self.LastYear+1)
        try:
            f = ma.masked_all((len(Years),self.nlev,self.nlat,self.nlon))
        except:
            f = ma.masked_all((len(Years),self.nlat,self.nlon))            
        for i in range(len(Years)):
            f[i] = self.getSeason(Years[i],Season)
        return ma.array(f).mean(axis=0)

    def _getBaseYear(self,File):
        # base year for time computation
        self.time_units = File.variables[self.FieldNames['time'][0]].units
        Year0 = self.time_units.split()[2].split('-')[0]
        if len(Year0) != 4:
            Year0 = self.time_units.split()[2].split('-')[-1]
        return int(Year0)


def CreateOutputFile(FileName,data):
    # Create file
    File = NetCDFFile(FileName,'w',format='NETCDF3_CLASSIC')
    # Define some global attribs
    File.Conventions='COARDS'
    # Time is record dimension
    File.createDimension('time',None)
    var = File.createVariable('time','d',('time',))
    var.long_name = 'time'
    var.units = 'year'
    # axes
    try:
        File.createDimension('level',data.nlev)
        var = File.createVariable('level','f',('level',))
        var.long_name = 'pressure level'
        var.units = 'mb'
        var[:] = data.File.variables['level'][:].astype('f')
    except: pass
    File.createDimension('lat',data.nlat)
    var = File.createVariable('lat','f',('lat',))
    var.long_name = 'latitude'
    var.units = 'degrees_north'
    var[:] = data.lat.astype('f')
    File.createDimension('lon',data.nlon)
    var = File.createVariable('lon','f',('lon',))
    var.long_name = 'longitude'
    var.units = 'degrees_east'
    var[:] = data.lon.astype('f')
    # create variables
    try:
        var = File.createVariable(Field,'f',('time','level','lat','lon'))
    except:
        var = File.createVariable(Field,'f',('time','lat','lon'))
    var.long_name = data.long_name
    var.units = data.units
    return File

def Seasonal(Field='U', Season='DJF', Source='ERA40', \
             YearStart=None, YearStop=None):
    # instatiate data server
    data = DataServer(Field=Field,Source=Source)
    if YearStart is None: YearStart = data.FirstYear
    if YearStop is None: YearStop = data.LastYear
    assert YearStart >= data.FirstYear,\
                       '\nFirst year in dataset is %s' % data.FirstYear
    assert YearStop <= data.LastYear,\
                       '\nLast year in dataset is %s' % data.LastYear
    # create output file
    FileName = '%s.%s.%s.%s-%s.nc' % (Field,Season,Source,YearStart,YearStop)
    File = CreateOutputFile(FileName,data)
    print 'Creating %s'%FileName
    TimeIndex = 0
    meter = ProgressMeter(total=YearStop-YearStart+1)
    for Year in range(YearStart,YearStop+1):
        meter.update(1)
        # get 1 season of data
        SeasonData = data.getSeason(Year,Season)
        File.variables['time'][TimeIndex]  = float(Year)
        File.variables[Field][TimeIndex] = SeasonData.astype('f')
        TimeIndex += 1
    File.close()

if __name__ == "__main__":

        from LambertProjector import *
        import matplotlib.pyplot as pl

	d = DataServer(Field='fsns')
	s = d.getMonth(2016,10)

	proj  = LambertProjector(boundinglat=60,resolution=100.)
	s     = proj(s,d.lon,d.lat)

	pl.contourf(proj.x,proj.y,s)
	proj.m.drawcoastlines()
	pl.colorbar()
	pl.show()

	print s.min(),s.max()

	"""
	LatRange = (70,82)
	LonRange = (20,90)
	years1   = range(1975,1998+1,1)
	years2   = range(1998,2014+1,1)

	proj  = LambertProjector(boundinglat=60,resolution=100.)
	mask1 = (proj.lat<LatRange[1])&(proj.lat>LatRange[0])&(proj.lon<LonRange[1])&(proj.lon>LonRange[0])==False
	d1    = DataServer(Field='ci',Source='NCAR20C_V2c')
	s1    = proj(np.array([d1.getSeason(Year=year,Season='DJF') for year in years1]),d1.lon,d1.lat)
	s2    = proj(np.array([d1.getSeason(Year=year,Season='DJF') for year in years2]),d1.lon,d1.lat)
	s1i   = np.ma.masked_array(s1,mask=np.tile(mask1,(len(years1),1,1))).mean(axis=1).mean(axis=1)
	s2i   = np.ma.masked_array(s2,mask=np.tile(mask1,(len(years2),1,1))).mean(axis=1).mean(axis=1)

	pl.figure(1)
	pl.plot(years1,s1i,'b',linewidth=2,alpha=0.6)
	pl.plot(years2,s2i,'r',linewidth=2,alpha=0.6)
	pl.grid()
	pl.xlabel('Year')
	pl.ylabel('%s [%s]' % (d1.long_name,d1.units))

	pl.figure(2)
	cseq = np.arange(-0.4,0.4+0.05,0.05)
	s    = s1.mean(axis=0)-s2.mean(axis=0)
	cf   = pl.contourf(proj.x,proj.y,s,cseq,cmap=pl.cm.coolwarm,extend='both')
	cbar = pl.colorbar(cf)
	cl1  = pl.contour(proj.x,proj.y,s1.mean(axis=0),[0.15],colors='b',linewidths=1.5,alpha=0.35)
	cl2  = pl.contour(proj.x,proj.y,s2.mean(axis=0),[0.15],colors='r',linewidths=1.5,alpha=0.35)
	cbar.set_label('%s [%s]' % (d1.long_name,d1.units))
	proj.m.drawparallels([70,80],latmax=90)
	proj.m.drawcoastlines()

	pl.show()
	"""

	"""
	from LambertProjector import *
	from mpl_toolkits.basemap import Basemap as Basemap
	import matplotlib.pyplot as pl

	proj = LambertProjector(boundinglat=40,resolution=100.)

	d4L  = DataServer(Field='Ts',Source='CAM4xCO2')
	d8L  = DataServer(Field='Ts',Source='CAM8xCO2')
	d16L = DataServer(Field='Ts',Source='CAM16xCO2')
	d4H  = DataServer(Field='Ts',Source='CAM4xCO2_high_continent')
	d8H  = DataServer(Field='Ts',Source='CAM8xCO2_high_continent')
	d16H = DataServer(Field='Ts',Source='CAM16xCO2_high_continent')

	years = range(1901,1919+1,1)
	s4L  = proj(array([d4L.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d4L.lon,d4L.lat)
	s8L  = proj(array([d8L.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d8L.lon,d8L.lat)
	s16L = proj(array([d16L.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d16L.lon,d16L.lat)
	s4H  = proj(array([d4H.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d4H.lon,d4H.lat)
	s8H  = proj(array([d8H.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d8H.lon,d8H.lat)
	s16H = proj(array([d16H.getSeason(Year=year,Season='NDJFM') for year in years]).mean(axis=0),d16H.lon,d16H.lat)

	cseq = np.arange(0,22+2,2)
	pl.figure(1)
	cf   = pl.contourf(proj.x,proj.y,s8L-s4L,cseq,cmap=pl.cm.OrRd,extend='both')
	cbar = pl.colorbar(cf)
	proj.m.drawparallels([70,80],latmax=90)
	pl.title('Low: 8xCO2 - 4xCO2')

        pl.figure(2)
        cf   = pl.contourf(proj.x,proj.y,s8H-s4H,cseq,cmap=pl.cm.OrRd,extend='both')
	cbar = pl.colorbar(cf)
        proj.m.drawparallels([70,80],latmax=90)
	pl.title('High: 8xCO2 - 4xCO2')

        pl.figure(3)
        cf   = pl.contourf(proj.x,proj.y,s16L-s8L,cseq,cmap=pl.cm.OrRd,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawparallels([70,80],latmax=90)
        pl.title('Low: 16xCO2 - 8xCO2')

        pl.figure(4)
        cf   = pl.contourf(proj.x,proj.y,s16H-s8H,cseq,cmap=pl.cm.OrRd,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawparallels([70,80],latmax=90)
        pl.title('High: 16xCO2 - 8xCO2')

	pl.show()
	"""
