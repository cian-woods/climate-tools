#!/usr/bin/env python

import os,sys,glob
from numpy import *
from netCDF4 import Dataset
import time,datetime
from progress import ProgressMeter
from interpol import vinterpol
import calendar
from UnPickle import *
from toPick import *

#HomeDir = os.getenv('HOME')
HomeDir = '/mnt/climstorage'

class DataServer:
    def __init__(	self, 
                 	Field       = 'U', 
                 	LevType     = 'plev',
			HybToPress  =  False,
			HybLevRange = (30,1000),
                 	Source      = 'ERAInt',
		 	DataFreq    = '6hourly',
                 	LevRange    = (0.,1000.), 
                 	LatRange    = (-90.,90.), 
                 	LonRange    = (0.,360.), 
		 	DailyMean   =  False,
			Year0       =  None		):
	# Attributes
        self.x,self.v,self.y = None,None,None
	self.LevType         = LevType
        self.Type            = LevType.split('_')[-1]
        if self.Type == 'forecast':
            self.shift = 6
        else:
            self.shift = 0
        self.FieldNames = {}
	self.leap       = True
	if Source[0:3] == 'ERA':
        	self.FieldNames['time'] 	= 'time'
        	self.FieldNames['lev'] 		= 'levelist'
        	self.FieldNames['lat'] 		= 'latitude'
        	self.FieldNames['lon'] 		= 'longitude'
        	self.FieldNames['slp'] 		= 'msl'
		self.FieldNames['msl_d'] 	= 'msl'
        	self.FieldNames['U10'] 		= 'u10'
		self.FieldNames['V10']          = 'v10'
        	self.FieldNames['U'] 		= 'u'
        	self.FieldNames['V'] 		= 'v'
        	self.FieldNames['W'] 		= 'w'
        	self.FieldNames['Z'] 		= 'z'
        	self.FieldNames['T'] 		= 't'
        	self.FieldNames['T2'] 		= 't2m'
        	self.FieldNames['Ts'] 		= 'skt'
        	self.FieldNames['pw'] 		= 'tcw'
		self.FieldNames['pv']           = 'pv'
		self.FieldNames['pwv'] 		= 'tcwv'
        	self.FieldNames['q'] 		= 'q'
        	self.FieldNames['ci'] 		= 'ci'
		self.FieldNames['lnsp']         = 'lnsp'
        	self.FieldNames['precc'] 	= 'cp'
        	self.FieldNames['precl'] 	= 'lsp'
		self.FieldNames['pottemp'] 	= 'pt'
		self.FieldNames['flds'] 	= 'strd'
		self.FieldNames['sshf'] 	= 'sshf'
		self.FieldNames['slhf'] 	= 'slhf'
		self.FieldNames['fls'] 		= 'str'
		self.FieldNames['fsns']         = 'ssr'
		self.FieldNames['z_filtered'] 	= 'z_filtered'
		self.FieldNames['msl_filtered'] = 'msl_filtered'
		self.FieldNames['shum'] 	= 'shum'
		self.FieldNames['ttr']		= 'ttr'
		self.FieldNames['tcc']		= 'cld'
		self.FieldNames['Fdiv']         = 'Fdiv'
		self.FieldNames['Fphi']         = 'Fphi'
		self.FieldNames['Fp']           = 'Fp'
		self.FieldNames['vq']           = 'vq'
		self.FieldNames['uq']           = 'uq'
		if Source == 'ERAInt':
			if (Field == 'ci') and (DataFreq == 'monthly'):
				Dir = '/qnap/cian/ci'
			elif LevType == 'modlev':
				Dir = HomeDir+'/obs/%s/%s/%s/60-90N/%s' % (Source,DataFreq,LevType,self.FieldNames[Field])
			elif (Field == 'Fdiv') or (Field == 'Fphi') or (Field == 'Fp'):
				Dir = '%s/cian/EPFD' % (HomeDir)
			elif Field == 'fsns':
				Dir = HomeDir+'/obs/%s/%s/%s/%s' % (Source,DataFreq,LevType,'fsns')
			elif (Field == 'vq') or (Field == 'uq'):
				Dir = '/mnt/climstorage/cian/moistflux/%s/' % (Field)
			else:
				Dir = HomeDir+'/obs/%s/%s/%s/%s' % (Source,DataFreq,LevType,self.FieldNames[Field])
		if Source == 'ERA20c':
			Dir = HomeDir+'/obs/%s/%s' % (Source,self.FieldNames[Field])
		FileNames = glob.glob(Dir+'/*.nc')
		if Field == 'fsns':
			FileNames.pop(FileNames.index('%s/fsns_1994.nc' % (Dir)))
			FileNames.append('%s/cian/fsns_1994.nc' % (HomeDir))	
        elif Source[0:4] == 'NCEP':
                self.FieldNames['lat']   = 'lat'
                self.FieldNames['lon']   = 'lon'
                self.FieldNames['shum']  = 'shum'
		self.FieldNames['rhum']  = 'rhum'
		self.FieldNames['vplev'] = 'vwnd'
		self.FieldNames['air_t'] = 'air'
		self.FieldNames['omega'] = 'omega'
		self.FieldNames['uplev'] = 'uwnd'
		self.FieldNames['t2m']   = 'air'
		self.FieldNames['dlwrf'] = 'dlwrf'
		self.FieldNames['vq']    = 'vq'
		self.FieldNames['uq']    = 'uq'
		self.FieldNames['nswrs'] = 'nswrs'
		self.FieldNames['nlwrs'] = 'nlwrs'
		self.FieldNames['slp']   = 'slp'
		Dir = '/mnt/climstorage/gabriele/obs/%s/%s/%s/%s' % (Source,DataFreq,LevType,Field)
		if (Source[-1] == '2') and (Field == 'shum'):
			Dir = '/mnt/climstorage/cian/NCAR2/shum'
		if (Field == 'vq') or (Field == 'uq'):
			Dir ='/mnt/climstorage/cian/NCAR1/%s' % (Field)
		elif DataFreq == 'monthly':
			Dir = '/mnt/climstorage/cian/NCAR1/monthly/%s/%s' % (LevType,Field)
		FileNames = glob.glob(Dir+'/*.nc')
	elif Source[0:5] == 'MERRA':
                self.FieldNames['lat']   = 'lat'
                self.FieldNames['lon']   = 'lon'
                self.FieldNames['QV']    = 'QV'
                self.FieldNames['V']     = 'V'
		Dir                      = '/mnt/climstorage/cian/MERRA/M%s/%s' % (Source[-1],Source)
		FileNames                = glob.glob(Dir+'/*')
	elif Source[0:3] == 'CAM':
		self.ExpType            = Source[3:]
		self.leap               = False
		self.FieldNames['lat']  = 'lat'
		self.FieldNames['lon']  = 'lon'
		self.FieldNames['lev']  = 'lev'
		self.FieldNames['Ts']   = 'TS'
		self.FieldNames['U']    = 'U'
                self.FieldNames['V']    = 'V'
		self.FieldNames['q']    = 'Q'
                self.FieldNames['T']    = 'T'
		self.FieldNames['W']    = 'OMEGA'
                self.FieldNames['slp']  = 'PSL'
		self.FieldNames['pw']   = 'TMQ'
		self.FieldNames['fls']  = 'FLNS'
		self.FieldNames['sshf'] = 'SHFLX'
		self.FieldNames['slhf'] = 'LHFLX'
		Dir                     = '/mnt/climstorage/cian/EXPS/Expt_%s/h1/%s' % (self.ExpType,self.FieldNames[Field])
		FileNames               = glob.glob(Dir+'/*.nc')
	print Dir
        # Dictionary of files
	Handles      = [Dataset(Name) for Name in FileNames]
	self.Handles = Handles[:]
	if Source[0:3] == 'CAM':
		tname = '/mnt/climstorage/cian/scripts/times/%s.%s.h1.p' % (Source,Field)
		if not os.path.isfile(tname):
			Times = [Handle.variables['time'][:] for Handle in Handles]
			toPick(Times,tname)
		print tname
		Times = unpick(tname)
	else:
		Times = [Handle.variables['time'][:] for Handle in Handles]
	#Times = []
	#for ii in range(len(Handles)):
	#	print FileNames[ii]
	#	Times.append(Handles[ii].variables['time'][:])
        self.Files = dict(zip(range(len(Handles)),Handles))
	self.Names = dict(zip(range(len(Handles)),FileNames))
        self.Times = dict(zip(range(len(Times)),Times))
        # Base year for time computation
        self.time_units    = Handles[0].variables['time'].units
	try:
		self.time_calendar = Handles[0].variables['time'].calendar
	except:
		self.time_calendar = 'N/A'
	if Year0 != None:
		self.Year0 = Year0
	else:
        	Year0 = self.time_units.split()[2]
        	self.Year0 = Year0.split('-')[0]
        	if len(self.Year0) != 4: self.Year0 = Year0.split('-')[-1]
        	self.Year0 = int(self.Year0)
	print self.Year0
        # Period covered
        self.MinHour = array([Time[0] for Time in Times]).min()
        self.MaxHour = array([Time[-1] for Time in Times]).max()
        self.MinDate = self.getDate(self.MinHour)
        self.MaxDate = self.getDate(self.MaxHour)
        print 'Data from ',self.MinDate,' to ',self.MaxDate
        # Initialize field 
        self.Field = Field
        self.FieldName = self.FieldNames[Field]
        if Field == 'tcc': self.FieldName = 'tcc'
	self.units = Handles[0].variables[self.FieldName].units
	try:
        	self.long_name = Handles[0].variables[self.FieldName].long_name
	except:
		self.long_name = 'N/A'
        # Initialize coord axes
	self.InvertLatAxis,self.InvertLonAxis = False,False
	self.j0,self.j1,self.i0,self.i1       = None,None,None,None
	try:
		latnames = ['latitude','lat']
		latname  = list(set(latnames) & set(Handles[0].variables.keys()))[0]
        	lat = array(Handles[0].variables[latname][:])*1.
        	(self.lat, self.j0, self.j1, self.InvertLatAxis) = self._setupAxis(lat,LatRange)
        	self.nlat = len(self.lat)
	except:
		pass
	try:
                lonnames = ['longitude','lon']
                lonname  = list(set(lonnames) & set(Handles[0].variables.keys()))[0]
        	lon = array(Handles[0].variables[lonname][:])*1.
        	(self.lon, self.i0, self.i1, self.InvertLonAxis) = self._setupAxis(lon,LonRange)
        	self.nlon = len(self.lon)
	except:
		pass
	# Level axis if it exists
	try:
		levnames = ['lev','levelist','level']
		levname  = list(set(levnames) & set(Handles[0].variables.keys()))[0]
		lev = array(Handles[0].variables[levname][:])*1.
		(self.lev, self.k0, self.k1, self.InvertLevAxis) = self._setupAxis(lev,LevRange)
		self.nlev = len(self.lev)	
	except:
		self.lev = None
	# If EP data, alternate the slice indexes
	if (Field == 'Fdiv') or (Field == 'Fphi') or (Field == 'Fp'):
		self.i0,self.i1   = self.j0,self.j1
		self.j0,self.j1   = self.k0,self.k1
	# Deal with modlev data
	self.Interpolate = False
	if (self.LevType == 'modlev') and (Field != 'lnsp'):
	    d = Dataset('/mnt/climstorage/obs/ERAInt/ERA-Interim_coordvars.nc','r')
	    self.hyam = d.variables['a_model_alt'][self.k0:self.k1]
	    self.hybm = d.variables['b_model_alt'][self.k0:self.k1]
	    d.close()
	    self.FillValue = 1.e20
            if HybToPress:
                self.Interpolate = True
                self.ps_data     = DataServer(Field='lnsp',Source=Source,LevType='modlev',HybToPress=False)
		self.p           = linspace(HybLevRange[0],HybLevRange[1],300)
		self.dP          = diff(self.p)[0]*100
	    else:
		self.p = (self.hyam+self.hybm*105000.)/100.
	# Land sea mask for ERAInt
	if Source == 'ERAInt':
		lsm_file = Dataset('/mnt/climstorage/obs/ERAInt/ERAInterimLandSeaMask_0.75x0.75.nc','r')
		lsm_lat  = lsm_file.variables[self.FieldNames['lat']][:]
		lsm_lon  = lsm_file.variables[self.FieldNames['lon']][:]
		self.lsm = lsm_file.variables['lsm'][:].squeeze()
		if self.InvertLatAxis: self.lsm = self.lsm[::-1,:]
		if self.InvertLonAxis: self.lsm = self.lsm[:,::-1]
	self.closeFiles()

    def _setupAxis(self,axis,range):
        if axis[0] > axis[-1]:
            axis = axis[::-1]
            invert = True
        else: invert = False
        i0  = argmin(abs(axis-range[0]))
        i1  = argmin(abs(axis-range[1]))+1
        axis = axis[i0:i1]
        return axis,i0,i1,invert


    def closeFiles(self):
        for File in self.Handles: File.close()

    def setPlevs(self,p1=None,p2=None):
    # Only 1 pressure level: Use p1 and p2 = None. If using 2 levels: p1 < p2.
        P1 = None
        P2 = None
        if p1 != None:
                x = argmin((self.lev-p1)**2)
                P1 = self.lev[x]
                if p2 == None:
                        v = x + 1
                        if v < len(self.lev):
                                P2 = self.lev[v]
                        else:
                                P2 = self.lev[x]
 #                      y = 0
                else:
                        v = argmin((self.lev-p2)**2)
                        P2 = self.lev[v]
                        v = v+1
        else:
#               x,v,y = None,None,None
                x,v = None,None
#       self.x,self.v,self.y = x,v,y
        self.x,self.v = x,v
#       print x,v,P1,P2
        lp = self.lev[x:v]
        return P1,lp

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
        	Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(nd*Hours/24.)
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

    def getDate(self,Hours,nd=1):
        # Given hours elapsed since midnight on 1 January Year0,
        # returns date (year, month, day, hour)
        if self.leap == False:
                Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days=nd*Hours/24.)
		date = (Date.year,Date.month,Date.day,Date.hour)
		n    = self.leapDaysBetweenDates(date)
		N    = 0
		while n-N > 0:
			N    = n
			Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days= n + nd*Hours/24.)
			date = (Date.year,Date.month,Date.day,Date.hour)
			n    = self.leapDaysBetweenDates(date)
        elif self.leap == True:
                Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(nd*Hours/24.)
        return Date.year,Date.month,Date.day,Date.hour

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
        	Days = datetime.datetime(int(Year),int(Month),int(Day),int(Hour)) \
        	       - datetime.datetime(self.Year0,1,1,0)
        	Hours = Days.days*24. + Days.seconds/3600.
        return Hours

    def getDays(self,Year,Month,Day,Hour):
		hours = self.getHours(Year,Month,Day,Hour)
		days  = hours/24.
		return days
    
    def getNdaysInMonth(self,Year,Month):
        if Month==12:
            Ndays = 31
        else:
            h1 = self.getHours(Year,Month  ,1,0)
            h2 = self.getHours(Year,Month+1,1,0)
            Ndays = int((h2-h1)/24)
        return Ndays

    def getDateList(self,Year=None, Month=None, Day=None, Hour=None, \
                    Season=None, step=6.):
        if Year is not None:
            if Month is not None:
                if Day is not None:
                    if Hour is not None:
                        Nsnapshots = 1
                        h = self.getHours(Year,Month,Day,Hour)
                    else: 
                        Nsnapshots = 4
                        h = self.getHours(Year,Month,Day,0)
                else:
                    Nsnapshots = self.getNdaysInMonth(Year,Month)*4
                    h = self.getHours(Year,Month,1,0)
            if Season is not None:
                if Season == 'DJF':
                    Months = [12,1,2]
                    h = self.getHours(Year-1,12,1,0)
                if Season == 'NDJF':
                    Months = [11,12,1,2]
                    h = self.getHours(Year-1,11,1,0)
                if Season == 'ONDJF':
                    Months = [10,11,12,1,2]
                    h = self.getHours(Year-1,10,1,0)
                if Season == 'ONDJFM':
                    Months = [10,11,12,1,2,3]
                    h = self.getHours(Year-1,10,1,0)
                if Season == 'JFM':
                    Months = [1,2,3]
                    h = self.getHours(Year,1,1,0)
                if Season == 'JFMA':
                    Months = [1,2,3,4]
                    h = self.getHours(Year,1,1,0)
                if Season == 'DJ':
                    Months = [12,1]
                    h = self.getHours(Year-1,12,1,0)
                if Season == 'FM':
                    Months = [2,3]
                    h = self.getHours(Year,2,1,0)
                if Season == 'NDJ':
                    Months = [11,12,1]
                    h = self.getHours(Year-1,11,1,0)
		if Season == 'ONDJ':
		    Months = [10,11,12,1]
		    h = self.getHours(Year-1,10,1,0)
                if Season == 'NDJFM':
                    Months = [11,12,1,2,3]
                    h = self.getHours(Year-1,11,1,0)
                if Season == 'ND':
                    Months = [11,12]
                    h = self.getHours(Year,11,1,0)
                if Season == 'MAM':
                    Months = [3,4,5]
                    h = self.getHours(Year,3,1,0)
                if Season == 'JJA':
                    Months = [6,7,8]
                    h = self.getHours(Year,6,1,0)
                if Season == 'SON':
                    Months = [9,10,11]
                    h = self.getHours(Year,9,1,0)
                if Season == 'OND':
                    Months = [10,11,12]
                    h = self.getHours(Year,10,1,0)
                if Season == 'Annual':
                    Months = [1,2,3,4,5,6,7,8,9,10,11,12]
                    h = self.getHours(Year,1,1,0)
                if Season == 'Annual_wc':
                    Months = [7,8,9,10,11,12,1,2,3,4,5,6]
                    h = self.getHours(Year-1,1,1,0)
                if Season == 'F':
                    Months = [2]
                    h = self.getHours(Year,2,1,0)
                if Season == 'D':
                    Months = [12]
                    h = self.getHours(Year-1,12,1,0)
                if Season == 'M':
                    Months = [3]
                    h = self.getHours(Year,3,1,0)
                if Season == 'J':
                    Months = [1]
                    h = self.getHours(Year,1,1,0)
                Nsnapshots = 0
                for Month in Months:
                    Nsnapshots += self.getNdaysInMonth(Year,Month)*4

	Nsnapshots = int(Nsnapshots*(6./step))
        Dates = []
        for i in range(Nsnapshots):
            Dates.append( self.getDate(h) )
            h += step
        return Dates

    def getHourList(self,Year=None,Month=None,Day=None,Hour=None,Season=None,step=6.):
	datelist = self.getDateList(Year,Month,Day,Hour,Season,step)
	return [self.getHours(*date) for date in datelist]

    def snapshot(self, Year=0, Month=1, Day=1, Hour=0):
        """
        # Extracts a single snapshot of Field.
        # Note that netCDF4 will automatically
        # - scale/offset values
        # - output a masked array if nc file has fill value specify
        """
        # Select file and time index
        now  = self.getHours(Year,Month,Day,Hour)
	now  = now + self.shift
        if now < self.MinHour or now > self.MaxHour:
            raise ValueError('Date ',self.getDate(now),' not in dataset!!')
        for key in self.Names:
            if now in self.Times[key]:
		File = Dataset(self.Names[key],'r')
                l    = argmin(abs(self.Times[key]-now))
                break
        # Retrieve variable
        f = File.variables[self.FieldName][l][self.x:self.v]
        f = squeeze(f)
        if self.Type == 'forecast':
            hour = self.getDate(self.Times[key][l])[-1]
            if (hour == 0) or (hour == 12):
                f = (f - squeeze(File.variables[self.FieldName][l-1][self.x:self.v]))/21600.
            else:
                f = f/21600. # divied by 6 hours in seconds
        # Swap axes if necessary and select slab
        if len(f.shape) == 2:
            if self.InvertLatAxis: f = f[::-1,:]
            if self.InvertLonAxis: f = f[:,::-1]
            f = f[self.j0:self.j1, self.i0:self.i1]
        if len(f.shape) == 3:
            if self.InvertLatAxis: f = f[:,::-1,:]
            if self.InvertLonAxis: f = f[:,:,::-1]
            try:
                if self.InvertLevAxis: f = f[::-1,:,:]
            except:
                pass
            f = f[self.k0:self.k1, self.j0:self.j1, self.i0:self.i1]
        # Interpolate from model to pressure coords if required
        if self.Interpolate:
            lnsp = self.ps_data.snapshot(Year,Month,Day,Hour) 
            ps   = exp(lnsp)
            f    = self.interpolate(f,ps,self.p)
	f = squeeze(f)
	File.close()
        return f

    def interpolate(self,xold,ps,p):
        pnew = p[:,None,None] + ps[None,:,:]*0.
        pold = (self.hyam[:,None,None]+self.hybm[:,None,None]*ps[None,:,:])/100.
        xnew = vinterpol(pold,xold,pnew,Extrapolate=False,FillValue=self.FillValue)
        return xnew

    def getDataSnaps(self, Year=None, Month=None, Day=None, Hour=None, Season=None, datelist=None, dailymean=False, step=6.):
	if dailymean == False:
	        if datelist == None:
			datelist = self.getDateList(Year,Month,Day,Hour,Season,step=step)
     	  	f = ma.array([self.snapshot(*date) for date in datelist])
      		f = f.squeeze()
	elif dailymean == True:
		f = self.getDailyMeanSnaps(Year,Month,Day,Hour,Season,datelist)
        return f

    def getDailyMeanSnaps(self, Year=None, Month=None, Day=None, Hour=None, Season=None, datelist=None):
	if datelist == None:
		datelist = self.getDateList(Year,Month,Day,Hour=None,Season=Season,step=24.)
       	f = array([self.dailyMeanData(date[0],date[1],date[2]) for date in datelist])
        f = f.squeeze()
        return f

    def dailyMeanData(self, Year=1990, Month=1, Day=1):
	datelist = self.getDateList(Year,Month,Day)
	f = array([self.snapshot(*date) for date in datelist])
	f = mean(f,axis=0)
	return f

    def getDay(self, Year=1958, Month=1, Day=1, Daily=None, TimeZoneOffset=0):
        # Return 1 day of data. TimeZoneOffset = -XX means time zone
        # is XX hours behind (earlier) than UTC
        f = []
        if TimeZoneOffset > 0 and (Year,Month,Day) == (1958,1,1): \
           TimeZoneOffset = 0
        if TimeZoneOffset < 0 and (Year,Month,Day) == (2001,12,31): \
           TimeZoneOffset = 0 
        dates = self.getDateList(Year,Month,Day)
        for date in dates:
            h    = self.getHours(*date) - TimeZoneOffset
            date = self.getDate(h)
            f.append(self.snapshot(*date))
        f = array(f)
        if   Daily is None:  return f
        elif Daily == 'max': return f.max(axis=0)
        elif Daily == 'min': return f.min(axis=0)
        elif Daily == 'avg': return f.mean(axis=0)
        else: raise ValueError,'operation %s not recognized' % Daily

    def getMonth(self, Year=1958, Month=1, Daily=None, TimeZoneOffset=0):
        # Return 1 month of data.
        # Keep iterating over day of month until exception
        # occurs, marking end of month
        print 'Getting %s %s %s' %(self.Field,Year,Month)
        Ndays = self.getNdaysInMonth(Year,Month)
        meter = ProgressMeter(total=Ndays)        
        f = []
        for Day in range(1,Ndays+1):
            meter.update(1)
            x = self.getDay(Year,Month,Day,Daily,TimeZoneOffset)
            if Daily is None: f.extend( x.tolist() )
            else: f.append(x)
        return array(f)

    def getSeason(self, Year=1959, Season='DJF'):
        if Season == 'DJF':
            DateStart = (Year-1,12,1,0)
            DateEnd   = (Year,2,28,18)
        if Season == 'MAM': 
            DateStart = (Year,3,1,0)
            DateEnd   = (Year,5,31,18)
        if Season == 'JJA': 
            DateStart = (Year,6,1,0)
            DateEnd   = (Year,8,31,18)
        if Season == 'SON': 
            DateStart = (Year,9,1,0)
            DateEnd   = (Year,11,30,18)
	if Season == 'M':
	    DateStart = (Year,3,1,0)
	    DateEnd   = (Year,3,31,18)
        return self.getTimeSlice(DateStart,DateEnd)        

    def getData(self, Year=1958, Month=1, Day=None, Hour=None, Season=None,\
                TimeZoneOffset=0, Daily=None):
        if Season is not None:
            return self.getSeason(Year,Season,Daily,TimeZoneOffset)
        if Hour is None:
            if Day is None:
                return self.getMonth(Year,Month,Daily,TimeZoneOffset)
            else:
                return self.getDay(Year,Month,Day,Daily=Daily)
        else:
            return self.snapshot(Year,Month,Day,Hour)
        
    def getTimeSlice(self, DateStart = (1958,1,1,0), DateEnd = (1958,12,31,18) ):
        print ' -- Getting timeslice %s to %s' % (DateStart,DateEnd)
        h0 = self.getHours(*DateStart)
        h1 = self.getHours(*DateEnd)        
        N = int((h1-h0)/6+1)
        f = self.snapshot(*self.getDate(h0))
        shape = (N,) + f.shape
        if hasattr(f,'mask'):
            f = ma.zeros(shape,dtype=float)
        else:
            f = zeros(shape,dtype=float)
        meter = ProgressMeter(total=N)
        for l in range(N):
            meter.update(1)
            f[l] = self.snapshot(*self.getDate(h0)) 
            h0 += 6
        return f

if __name__ == '__main__':

	import matplotlib.pyplot as pl
	import numpy as np
	from LambertProjector import *
	from UnPickle import *

	df     = DataServer(Field='T',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
	years  = range(1979,2016+1,1)
	months = range(1,12+1,1)
	for year in years:
		for month in months:
			print year,month
			F = df.getDataSnaps(Year=year,Month=month,dailymean=True)
			print F.shape

	"""
	# DataServer and get data
	d1   = DataServer(Field='Fdiv',LatRange=(60,90),LonRange=(0,360),LevRange=(10,500))
	d2   = DataServer(Field='Fphi',LatRange=(60,90),LonRange=(0,360),LevRange=(10,500))
	d3   = DataServer(Field='Fp',LatRange=(60,90),LonRange=(0,360),LevRange=(10,500))


	# Clims
	Years    = range(1991,2016+1,1)
	divclim  = np.array([d1.getDataSnaps(Year,Season='NDJFM').mean(axis=0) for Year in Years]).mean(axis=0)
	Fphiclim = np.array([d2.getDataSnaps(Year,Season='NDJFM').mean(axis=0) for Year in Years]).mean(axis=0)
	Fpclim   = np.array([d3.getDataSnaps(Year,Season='NDJFM').mean(axis=0) for Year in Years]).mean(axis=0)

	days  = np.arange(-10,0+1,1)
	dates = unpick('../WarmArctic/dates.cold.50.ERAInt.p')
	div,Fphi,Fp = [],[],[]
	for date0 in dates:
		hour0    = d1.getHours(*date0)
		datelist = [d1.getDate(hour0 + ii*24) for ii in days]
		try:
			div.append([d1.snapshot(*date) for date in datelist])
			Fphi.append([d2.snapshot(*date) for date in datelist])
			Fp.append([d3.snapshot(*date) for date in datelist])
		except:
			pass
	div,Fphi,Fp = np.array(div).mean(axis=0),np.array(Fphi).mean(axis=0),np.array(Fp).mean(axis=0)
	div,Fphi,Fp = div-divclim[np.newaxis,:,:],Fphi-Fphiclim[np.newaxis,:,:],Fp-Fpclim[np.newaxis,:,:]
	print div.shape

        sid    = 24*3600
        lat    = d1.lat
        lev    = d1.lev
        loglev = np.exp(np.linspace(np.log(10),np.log(500),20))
        xs     = [np.argmin((lev-loglev[ii])**2) for ii in range(len(loglev))]
        N      = np.ones(len(lev))
        N[np.where(lev<100)] = 3
        N      = N[np.newaxis,:,np.newaxis]

	div  = -1*sid*div
	Fphi = Fphi*N
	Fp   = -1*Fp*N
	cseq = np.arange(-14,14+2,2)
	for ii in range(len(div)):
	        cf   = pl.contourf(lat,lev,div[ii],cseq,cmap=pl.cm.RdBu_r,extend='both')
	        cbar = pl.colorbar()
	        cbar.set_label('Zonal wind acceleration anomaly [m s$^{-1}$ day$^{-1}$]')

	        Q  = pl.quiver(lat[::5],lev[xs],Fphi[ii,xs,::5],Fp[ii,xs,::5],units='inches',scale=100,\
	                        scale_units='inches',headwidth=2.5,headlength=2.5,headaxislength=2.5,pivot='tail',alpha=0.6)
	        qk = pl.quiverkey(Q, 0.2, 1.02, 20, '%s' % (20), labelpos='W',fontproperties={'weight': 'bold'})

	        pl.ylabel('Pressure [hPa]')
	        pl.xlabel(r'Latitude [$^{\circ}$ N]')
	        pl.ylim(1000,10)
	        pl.yscale('log')
	        pl.xlim(-90,90)
	        pl.title('day %s' % (days[ii]))
	        pl.show()
	"""
