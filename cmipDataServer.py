#!/usr/bin/env python

import os,sys,glob
from UnPickle import *
from toPick import *
from netCDF4 import Dataset
import time,datetime
from progress import ProgressMeter
import calendar
import numpy as np


# ------------------------------------------------------------- #
#								#
# getDays and getDate have been rewritten to have no leap days	#
#								#
# --------------------------------------------------------------#

class DataServer:
    def __init__(	self, 
			Field    = 'prw', 
			LevType  = 'surface',
			Source   = 'CMCC-CESM',
			ExpType  = 'historical',
			DataFreq = 'mon',
			LevRange = (0,1000.), 
			LatRange = (-90.,90.), 
			LonRange = (0.,360.),
			close    = False		):

		if 'ERA' == 'ERA':
		    self.x,self.v,self.y = None,None,None 
		    self.FieldNames = {}
		    self.FieldNames['time']         = 'time'
		    self.FieldNames['plev']         = 'plev'
		    self.FieldNames['lat']          = 'lat'
		    self.FieldNames['lon']          = 'lon'
		    self.FieldNames['psl']          = 'psl'
		    self.FieldNames['ua']           = 'ua'
		    self.FieldNames['va']           = 'va'
		    self.FieldNames['W']            = 'w'
		    self.FieldNames['zg']           = 'zg'
		    self.FieldNames['T']            = 't'
		    self.FieldNames['T2']           = 't2m'
		    self.FieldNames['tas']          = 'tas'
		    self.FieldNames['sic']          = 'sic'
		    self.FieldNames['pw']           = 'tcw'
		    self.FieldNames['hus']          = 'hus'
		    self.FieldNames['rlds']         = 'rlds'
		    self.FieldNames['rlus']         = 'rlus'
		    self.FieldNames['hfss']         = 'hfss'
		    self.FieldNames['hfls']         = 'hfls'
		    self.FieldNames['ci']           = 'ci'
		    self.FieldNames['precc']        = 'cp'
		    self.FieldNames['precl']        = 'lsp'
		    self.FieldNames['pottemp']      = 'pt'
		    self.FieldNames['zg_filtered']  = 'zg_filtered'
		    self.FieldNames['psl_filtered'] = 'psl_filtered'
		    self.FieldNames['ta']           = 'ta'
		    self.FieldNames['rldscs']       = 'rldscs'
		    self.FieldNames['prw']          = 'prw'
		    self.FieldNames['wap']          = 'wap'
		    if ExpType == 'rcp85':
			homedir = '/mnt/climstorage/cian'
		    if ExpType == 'historical':
			homedir = '/mnt/climstorage/cian'#'/qnap/cian/cmip'
			if (Field == 'prw') and (DataFreq == 'day'):
				homedir = '/mnt/climstorage/cian'
		    if ExpType == 'piControl':
			homedir = '/mnt/climstorage/cian'
		    Dir = '%s/%s/%s/%s/%s/%s' % (homedir,ExpType,Source,DataFreq,LevType,self.FieldNames[Field])
			
		    print Dir
		    # dictionary of files
		    Names         = glob.glob(Dir+'/*')
		    Handles       = [Dataset(Name) for Name in Names] 
		    self.Handles  = Handles[:]
		    self.Files    = dict(zip(range(len(Handles)),Handles))
		    self.calendar = Handles[0].variables['time'].calendar
		    if (self.calendar == 'proleptic_gregorian') or (self.calendar == 'standard') or (self.calendar == 'gregorian'):
			self.leap = True
		    else:
			self.leap = False
		    print self.calendar
		    # base year for time computation
		    Year0s = [int(Handle.variables['time'].units.split()[2].split('-')[0]) for Handle in Handles]
		    self.time_units = Handles[np.argmin(Year0s)].variables['time'].units
		    self.Year0  = min(Year0s)
		    print 'Year0 = %s' % (str(self.Year0))
		    tname = '/mnt/climstorage/cian/scripts/times/%s.%s.%s.%s.p' % (Source,Field,ExpType,DataFreq)
		    if os.path.isfile(tname) == False:
		    	Times = [Handle.variables['time'][:] for Handle in Handles]
			toPick(Times,tname)
		    else:
			print tname
			Times = unpick(tname)
		    for i in range(len(Year0s)):
			# Any timesteps = 0 in time axis are replaced with correct time step
			dt              = np.diff(Times[i])
			(values,counts) = np.unique(dt,return_counts=True)
			ind             = np.argmax(counts)
			xs              = np.where(Times[i]==0)[0]
			for jj in xs:
				Times[i][jj] = Times[i][jj-1] + values[ind]
			# Set all times to same base year
			start    = self.getDays(self.Year0,1,1,12)
			Times[i] = Times[i] + self.getDays(Year0s[i],1,1,12) - start
		    self.Times = dict(zip(range(len(Times)),Times))
		    # period covered
#		    self.MinDay  = np.array([Time[0] for Time in Times]).min()
#		    self.MaxDay  = np.array([Time[-1] for Time in Times]).max()
		    self.MinDay  = np.array([Time.min() for Time in Times]).min()
		    self.MaxDay  = np.array([Time.max() for Time in Times]).max()
		    self.MinDate = self.getDate(self.MinDay)
		    self.MaxDate = self.getDate(self.MaxDay)
		    self.Yearmin = self.MinDate[0]
		    print 'Data from ',self.MinDate,' to ',self.MaxDate

    		# Initialize field 
		self.DataFreq  = DataFreq
		self.Field     = Field
		self.FieldName = self.FieldNames[Field]
		self.units     = Handles[0].variables[self.FieldName].units
		self.long_name = Handles[0].variables[self.FieldName].long_name

		# Initialize coord axes
		lat = np.array(Handles[0].variables[self.FieldNames['lat']][:])*1.
		(self.lat, self.j0, self.j1, self.InvertLatAxis) = \
			   self._setupAxis(lat,LatRange)
		self.nlat = len(self.lat)
		lon = np.array(Handles[0].variables[self.FieldNames['lon']][:])*1.
		(self.lon, self.i0, self.i1, self.InvertLonAxis) = \
			   self._setupAxis(lon,LonRange)
		self.nlon = len(self.lon)
		try:
		    lev = self.Files[0].variables[self.FieldNames['plev']][:]/100.
		    (self.lev, self.k0, self.k1, self.InvertLevAxis) = \
			       self._setupAxis(lev,LevRange)
		    self.nlev = len(self.lev)
		except:
		    pass
		if close == True: self.closeFiles(Handles)

    def _setupAxis(self,axis,LevRange):
		if axis[0] > axis[-1]:
		    axis = axis[::-1]
		    invert = True
		else: invert = False
		i0  = np.argmin(abs(axis-LevRange[0]))
		i1  = np.argmin(abs(axis-LevRange[1]))+1
		axis = axis[i0:i1]
		return axis,i0,i1,invert

    def closeFiles(self):
		for File in self.Handles: File.close()
			
    def setPlevs(self,p1=None,p2=None):
    # Only 1 pressure level: Use p1 and p2 = None. If using 2 levels: p1 < p2.
        P1 = None
        P2 = None
        if p1 != None:
                x = np.argmin((self.lev-p1)**2)
                P1 = self.lev[x]
                if p2 == None:
                        v = x + 1
                        if v < len(self.lev):
                                P2 = self.lev[v]
                        else:
                                P2 = self.lev[x]
 #                      y = 0
                else:
                        v = np.argmin((self.lev-p2)**2)
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
    def getDate(self,Days,leap=False,nd=None):
		# Given days elapsed since midnight on 1 January Year0,
		# returns date (year, month, day, hour) 
		# Leap == False (do not include leap years)
		if self.leap == False:
			Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days=Days)
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
#				Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(Days) + datetime.timedelta(m)
#				Date  = Date + datetime.timedelta(days=m-S)
				if m-S == 0:
					break
				else:
					Date  = Date + datetime.timedelta(days=float(m-S))
					S = m
		elif self.leap == True:
			Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(Days)

		return Date.year,Date.month,Date.day,Date.hour
    """

    def getHours(self,Year,Month,Day,Hour):
		hours = self.getDays(Year,Month,Day,Hour,hours=True)
		return hours

    def leapDaysBetweenDates(self,date1):
        if (date1[1]==2 and date1[2]<=28) or (date1[1]==1):
                year1,dn1 = date1[0],0
        elif (date1[1]==2) and (28 < date1[2]):
                year1,dn1 = date1[0],(date1[3]+6)/24.
        else:
                year1,dn1 = date1[0]+1,0
        n = calendar.leapdays(self.Year0,year1) + dn1
        return n

    def getDate(self,Days,leap=False,nd=None):
        # Given hours elapsed since midnight on 1 January Year0,
        # returns date (year, month, day, hour)
        if self.leap == False:
                Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days=Days)
                date = (Date.year,Date.month,Date.day,Date.hour)
                n    = self.leapDaysBetweenDates(date)
                N    = 0
                while n-N > 0:
                        N    = n
                        Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(days= n + Days)
                        date = (Date.year,Date.month,Date.day,Date.hour)
                        n    = self.leapDaysBetweenDates(date)
        elif self.leap == True:
                Date = datetime.datetime(self.Year0,1,1,0) + datetime.timedelta(Days)
        return Date.year,Date.month,Date.day,Date.hour

    def getDays(self,Year,Month,Day,Hour,hours=False,leap=False):
		# Given date (year, month, day, hour),
		# returns Days elapsed since midnight on 1 January Year0
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
			Days  = Hours/24.

		elif self.leap == True:
			Days = datetime.datetime(int(Year),int(Month),int(Day),int(Hour)) \
                       - datetime.datetime(self.Year0,1,1,0)

			Hours = Days.days*24 + Days.seconds/3600
                	Days  = (Hours/24.)

		if hours == False:
	       	 	return Days
		elif hours == True:
			return Hours
    
    def getNdaysInMonth(self,Year,Month):
        if Month==12:
            Ndays = 31
        else:
            h1 = self.getDays(Year,Month,1,0)
            h2 = self.getDays(Year,Month+1,1,0)
            Ndays = int(h2-h1)
        return Ndays

    def getHourList(self,Year=None,Month=None,Day=None,Hour=None,Season=None,step=None):
        datelist = self.getDateList(Year,Month,Day,Hour,Season,step)
        return [self.getHours(*date) for date in datelist]

    def getDateList(self,Year=None, Month=None, Day=None, Hour=None, \
                    Season=None,step=None):

 	if (Year is not None) and (self.DataFreq == 'day'):
            if Month is not None:
                if Day is not None:
                    if Hour is not None:
                        Nsnapshots = 1
                        h = self.getDays(Year,Month,Day,Hour)
                    else: 
                        Nsnapshots = 1
                        h = self.getDays(Year,Month,Day,12)
                else:
                    Nsnapshots = self.getNdaysInMonth(Year,Month)
                    h = self.getDays(Year,Month,1,12)
	    else:
		Nsnapshots = 0
		for i in range(1,13):
			Nsnapshots += self.getNdaysInMonth(Year,i)
			h = self.getDays(Year,1,1,12)
            if Season is not None:
		Nsnapshots = 0
                if Season == 'DJF':
                    Months = [12,1,2]
                    h = self.getDays(Year-1,12,1,12)
                if Season == 'ONDJ':
                    Months = [10,11,12,1]
                    h = self.getDays(Year-1,10,1,12)
                if Season == 'ONDJF':
                    Months = [10,11,12,1,2]
                    h = self.getDays(Year-1,10,1,12)
                if Season == 'NDJ':
                    Months = [11,12,1]
                    h = self.getDays(Year-1,11,1,12)
                if Season == 'NDJFM':
                    Months = [11,12,1,2,3]
                    h = self.getDays(Year-1,11,1,12)
                if Season == 'NDJF':
                    Months = [11,12,1,2]
                    h = self.getDays(Year-1,11,1,12)
                if Season == 'NDJFMi':
                    Months = [11,12,1,2,3]
                    h = self.getDays(Year-1,11,2,12)
                if Season == 'MAM':
                    Months = [3,4,5]
                    h = self.getDays(Year,3,1,12)
                if Season == 'JJA':
                    Months = [6,7,8]
                    h = self.getDays(Year,6,1,12)
                if Season == 'SON':
                    Months = [9,10,11]
                    h = self.getDays(Year,9,1,12) 
                if Season == 'Annual':
                    Months = range(1,12+1,1)
                    h = self.getDays(Year,1,1,12)
                if Season == 'J':
                    Months = [1]
                    h = self.getDays(Year,1,1,12)
                Nsnapshots = 0
                for Month in Months:
                    Nsnapshots += self.getNdaysInMonth(Year,Month)

	    r = np.ones(Nsnapshots)

	elif (Year is not None) and (self.DataFreq == 'mon'):
	    if Month is not None:
		h = self.getDays(Year,Month,1,12)
		r = [1.]
		Nsnapshots = 1
	    elif Season is not None:
		if Season == 'DJF':
                    Months = [12,1,2]
                    h = self.getDays(Year-1,12,1,12)
		    r = [self.getNdaysInMonth(Year-1,12)]
		    r.append(self.getNdaysInMonth(Year,1))
		    r.append(self.getNdaysInMonth(Year,2))
                if Season == 'NDJF':
                    Months = [11,12,1,2]
                    h = self.getDays(Year-1,11,1,12)
                    r = [self.getNdaysInMonth(Year-1,11)]
		    r.append(self.getNdaysInMonth(Year-1,12))
                    r.append(self.getNdaysInMonth(Year,1))
                    r.append(self.getNdaysInMonth(Year,2))
                if Season == 'MAM':
                    Months = [3,4,5]
                    h = self.getDays(Year,3,1,12)
		    r = [self.getNdaysInMonth(Year,m) for m in Months]
                if Season == 'JJA':
                    Months = [6,7,8]
                    h = self.getDays(Year,6,1,12)
		    r = [self.getNdaysInMonth(Year,m) for m in Months]
                if Season == 'SON':
                    Months = [9,10,11]
                    h = self.getDays(Year,9,1,12)
		    r = [self.getNdaysInMonth(Year,m) for m in Months]
		Nsnapshots = 3
            else:
                h = self.getDays(Year,1,1,12)
                Nsnapshots = 12
                r = [self.getNdaysInMonth(Year,m) for m in range(1,13)]
        Dates = []
        for i in range(Nsnapshots):
            Dates.append( self.getDate(h) )
            h += r[i]
        return Dates

    def snapshot(self, Year=0, Month=1, Day=1, Hour=0):
        """
        # Extracts a single snapshot of Field.
        # Note that netCDF4 will automatically
        # - scale/offset values
        # - output a masked np.array if nc file has fill value specify
        """

	if self.DataFreq == 'mon':
#		n = self.getNdaysInMonth(Year,Month)/2.
#		Day,Hour = 1 + n,12
		H    = self.getDays(Year,Month,1,0)
		n    = self.getNdaysInMonth(Year,Month)/2.
		date = self.getDate(H+n)
		Year,Month,Day,Hour = date[0],date[1],date[2],date[3]

        # select file and time index
        now = self.getDays(Year,Month,Day,Hour)
        if now < self.MinDay or now > self.MaxDay:
            raise ValueError('Date ',self.getDate(now),' not in dataset!!')
        for key in self.Files:
	    if now in self.Times[key]:
                File = self.Files[key]
                l    = np.argmin(abs(self.Times[key]-now))
                break
	else:
		raise ValueError('Date ',self.getDate(now),' not in dataset!!')
        # retrieve variable
#       f = File.variables[self.FieldName][l][self.x:self.v]
#	f = np.squeeze(f)
	f = File.variables[self.FieldName][l] 
        # swap axes if necessary and select slab
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

	if type(f) == np.ma.core.MaskedArray:	
		np.ma.set_fill_value(f,0)
	f = np.squeeze(f)
#	if self.Field == 'rlus':
#		f = np.ma.masked_np.array(abs(f))
        return f

    def getDataSnaps(self, Year=None, Month=None, Day=None, Hour=None, Season=None, datelist=None, dailymean=None,step=None):
	if datelist == None:        
		datelist = self.getDateList(Year,Month,Day,Hour,Season,step)
        f = np.ma.array([self.snapshot(*date) for date in datelist])
        f = f.squeeze()	
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
            h    = self.getDays(*date) - TimeZoneOffset
            date = self.getDate(h)
            f.append(self.snapshot(*date))
        f = np.array(f)
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
        return np.array(f)

    def getSeason(self, Year=1959, Season='DJF'):
	"""
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
        #return self.getTimeSlice(DateStart,DateEnd)
	"""
	return self.getDataSnaps(Year=Year,Season=Season).mean(axis=0)

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
        h0 = self.getDays(*DateStart)
        h1 = self.getDays(*DateEnd)        
        N = int((h1-h0)/6+1)
        f = self.snapshot(*self.getDate(h0))
        shape = (N,) + f.shape
        if hasattr(f,'mask'):
            f = np.ma.zeros(shape,dtype=float)
        else:
            f = zeros(shape,dtype=float)
        meter = ProgressMeter(total=N)
        for l in range(N):
            meter.update(1)
            f[l] = self.snapshot(*self.getDate(h0)) 
            h0 += 6
        return f

if __name__ == '__main__':

	from LambertProjector import *
	from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
	import matplotlib.pyplot as pl

	Source = str(sys.argv[1])
	vds = DataServer(Field='hus' ,LevType='plev',DataFreq='mon',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
	print vds.lev
