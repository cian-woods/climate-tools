#!/usr/bin/env python

import sys
from numpy import *
from netCDF4 import Dataset as NetCDFFile
from mpl_toolkits.basemap import Basemap as Basemap
import time,datetime
from toPick import *
from UnPickle import *
from crossNetCDF import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from scipy import stats,interpolate
from LambertProjector import *
import calendar
import glob

class DataServer:
    def __init__(  self,
		   Source =  'ERAInt',
		   Type	  =    'fwrd',
		   blat   =        70,
		   steps  =         8	):

	if Source[0:3] == 'CAM': self.leap = False
	if Source[0:3] == 'ERA': self.leap = True
	self.Type  = Type
        self.x     = steps
        self.FieldNames = {}
        self.FieldNames['time']      = 'time'
        self.FieldNames['Latitude']  = 'Latitude'
        self.FieldNames['Longitude'] = 'Longitude'
	self.FieldNames['Pressure']  = 'Pressure'
	self.FieldNames['pv']        = 'pv'
	self.FieldNames['lon']       = 'Longitude'
	self.FieldNames['steps']     = 'steps'
	# Attributes
	self.Dir       = '/mnt/climstorage/cian/Tracks/%s/%s/%s' % (Source,Type,blat)
	print self.Dir
	self.Files     = {}
	FileName       = glob.glob('%s/%s_*.nc' % (self.Dir,Type))[0]
	self.datestamp = ['%s_%s' % (FileName.split('_')[-2],FileName.split('_')[-1][:-3])]
	self.Files[self.datestamp[0]] = NetCDFFile(FileName,'r')
	self.steps = self.Files[self.datestamp[0]].variables['steps'][:]
	# Initisalize some DataServer parameters
	self.Year0 = self.Files[self.datestamp[0]].variables['time'].units.split()[2]
        if len(self.Year0) != 4:
            self.Year0 = self.Year0.split('-')[0]
        self.Year0 = int(self.Year0)
        # Initialize coord axes
        File     = self.Files.values()[0]
	self.lon = array(File.variables[self.FieldNames['lon']][:])*1.
	# Mercator
	self.m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')


    def getTangents(self,x,y,mag=None):
        # Tracks should be in cartesian coordinates
        ys,xs = [],[]
        # Add tangent of first point
        y0,x0   = y[0],x[0]
        y1,x1   = y[1],x[1]
        dyt,dxt = y1-y0,x1-x0
        R       = np.sqrt(dyt**2 + dxt**2)
        if mag != None: dyt,dxt = dyt/R,dxt/R
        ys.append(dyt)
        xs.append(dxt)
        for i in range(len(y)-2):
                y0,x0 = y[i],x[i]
                y1,x1 = y[i+1],x[i+1]
                dy,dx = y1-y0,x1-x0
                y2,x2 = y[i+2]-dy,x[i+2]-dx
                yt,xt = (y1+y2)/2.,(x1+x2)/2.
                dyt,dxt = yt-y0,xt-x0
                R       = np.sqrt(dyt**2 + dxt**2)
                if mag != None: dyt,dxt = dyt/R,dxt/R
                ys.append(dyt)
                xs.append(dxt)
        # Add tanget for last point
        y0,x0 = y[-2],x[-2]
        y1,x1 = y[-1],x[-1]
        dyt,dxt = y1-y0,x1-x0
        R       = np.sqrt(dyt**2 + dxt**2)
        if mag != None: dyt,dxt = dyt/R,dxt/R
        ys.append(dyt)
        xs.append(dxt)
        return np.array(ys),np.array(xs)

    def trackDensity(self,xs,ys,xx,yy):
        # Centroid trajectories
        xs  = np.array(xs).mean(axis=0) # centroid
        ys  = np.array(ys).mean(axis=0) # centroid
        # Tangents
        yts,xts = self.getTangents(xs,ys,mag=1)
        # Data holders
        N   = np.zeros((len(xx),len(yy)))
        T   = np.zeros((len(xx),len(yy),2))
        ind = []
        for i in range(len(xs)):
                        xi,yi      = np.argmin(abs(xx-xs[i])),np.argmin(abs(yy-ys[i])) 
                        N[yi,xi]   = N[yi,xi] + 1
                        T[yi,xi,0] = T[yi,xi,0] + xts[i]
                        T[yi,xi,1] = T[yi,xi,1] + yts[i]
        Ntile = np.tile(N[:,:,np.newaxis],(1,1,2))
        T     = np.ma.masked_where(Ntile==0,T)/np.ma.masked_where(Ntile==0,Ntile)
        T     = T.data
        return N,T

    def density(self,LON,LAT,dates,proj,nend=None):
	# Lambert Projector
	xx,yy = proj.x[0,:],proj.y[:,0]
        # Intrusions trajectories
        nn,Ttot = [],[]
        for i in range(len(LON)): 
                N,T = [],[]
                for t in range(len(LON[i])):
                        lons,lats = LON[i][t],LAT[i][t]
                        xs,ys     = proj.m(lons,lats)
                        xs,ys     = xs[:,:nend],ys[:,:nend]
                        Nit,Tit   = self.trackDensity(xs,ys,xx,yy)
                        N.append(Nit)
                        Ttot.append(Tit)
                N = np.array(N).sum(axis=0)
                N[np.where(N>=1)] = 1
		#N = N/len(LON[i])
                nn.append(N)
        nn   = np.array(nn)
        Ntot = nn.sum(axis=0)
        Ttot = np.array(Ttot)
        Ttot = np.ma.masked_where(Ttot==0,Ttot).mean(axis=0)
	years,ns = self.seasonalIntrusions(nn,dates)
        return nn,ns,years,Ntot,Ttot,xx,yy

    def seasonalIntrusions(self,nn,dates):
	N = {}
	for ii in range(len(dates)):
		year,month = dates[ii][0],dates[ii][1]
		if month >= 9: year = year + 1
		try:
			N[year].append(nn[ii])
		except:
			N[year] = []
			N[year].append(nn[ii])
	return N.keys(),N.values()

    def getIntrusionTrajectories(self,G,D):
        LON,LAT,P = [],[],[]
        for ii in range(len(D)):
                LON.append([])
                LAT.append([])
                P.append([])
                for jj in range(len(D[ii])):
                        LON[ii].append([])
                        LAT[ii].append([])
                        P[ii].append([])
                        snap      = self.snapshot(*D[ii][jj])
                        lon,lat,p = snap[:,:,1],snap[:,:,0],snap[:,:,2]
                        for kk in range(len(G[ii][jj])): 
                                x = int(G[ii][jj][kk])
                                LON[ii][jj].append(lon[x,:])
                                LAT[ii][jj].append(lat[x,:])
                                P[ii][jj].append(p[x,:])
                        LON[ii][jj] = np.array(LON[ii][jj])
                        LAT[ii][jj] = np.array(LAT[ii][jj])
                        P[ii][jj]   = np.array(P[ii][jj])
        return LON,LAT,P

    def julianDay(self,hours,Year0=None):
        Year,Month,Day,Hour = self.getDate(hours)
        if self.getNdaysInMonth(Year,2) == 28:
		    months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            months = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        jday = 0
        if Month == 1:
        	jday = Day + Hour/24.
        else:
        	for i in range(Month-1):
                 jday = jday + months[i]
        	jday = jday + Day + (Hour/24.)
        m = 0
        if Year0 != None:
        	for i in range(Year0,Year):
        		m = m + self.getNdaysInYear(i)
        jday = jday + m
        return jday
        
    def julianToDate(self,jday,Year):
        hour = self.getHours(Year,1,1,0)
        now = hour+(jday-1)*24
        date = self.getDate(now)
        return date

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
        	Hours = Days.days*24 + Days.seconds/3600
        return Hours
    
    def getNdaysInMonth(self,Year,Month):
        if Month==12:
            Ndays = 31
        else:
            h1 = self.getHours(Year,Month  ,1,0)
            h2 = self.getHours(Year,Month+1,1,0)
            Ndays = int((h2-h1)/24)
        return Ndays
        
    def getNdaysInYear(self,Year):
		if self.getNdaysInMonth(Year,2) == 28:
			Ndays = 365
		else:
			Ndays = 366
		return Ndays

    def getNdaysInSeason(self,Year,Season='DJF'):
		# Year is the year of December in the case of DJF.
		Ndays = 0
		if Season == 'DJF':
			Ndays = self.getNdaysInMonth(Year,2)+62
		elif Season == 'MAM':
			Ndays = 92
		elif Season == 'JJA':
			Ndays = 92
		else:
			Ndays = 91
		return Ndays

    def getDateList(self,Year=None, Month=None, Day=None, Hour=None, \
                    Season=None):
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
            else:
            	Nsnapshots = self.getNdaysInYear(Year)*4
            	h = self.getHours(Year,1,1,0)
            if Season is not None:
                if Season == 'DJF':
                    Months = [12,1,2]
                    h = self.getHours(Year-1,12,1,0)
                if Season == 'ON':
                    Months = [10,11]
                    h = self.getHours(Year-1,10,1,0)
                if Season == 'ONDJ':
                    Months = [10,11,12,1]
                    h = self.getHours(Year-1,10,1,0)
                if Season == 'MAM':
                    Months = [3,4,5]
                    h = self.getHours(Year,3,1,0)
                if Season == 'JJA':
                    Months = [6,7,8]
                    h = self.getHours(Year,6,1,0)
                if Season == 'SON':
                    Months = [9,10,11]
                    h = self.getHours(Year,9,1,0)
                Nsnapshots = 0
                for Month in Months:
                    Nsnapshots += self.getNdaysInMonth(Year,Month)*4
        Dates = []
        for i in range(Nsnapshots):
            Dates.append( self.getDate(h) )
            h += 6.
        return Dates
        
    def getHourList(self,Year=None,Month=None,Day=None,Hour=None,Season=None):
    	dates = self.getDateList(Year,Month,Day,Hour,Season)
    	hours = []
    	for i in dates:
    		hours.append(self.getHours(*i))
    	return hours

    def snapshot(self,Year=0,Month=1,Day=1,Hour=0):
	YearMonth = '%4i_%02i' % (Year,Month)
        FileName  = '%s/%s_%s.nc' % (self.Dir,self.Type,YearMonth)
#	print FileName
	try:
		File      = NetCDFFile(FileName,'r')
	except:
		print 'File %s does not exist!' % (FileName)
		sys.exit()
        # select time 
        time = File.variables[self.FieldNames['time']][:]
        now  = self.getHours(Year,Month,Day,Hour)
        l    = argmin(abs(time-now))
        # retrieve variable
	if now == time[l]:
		f1 = squeeze(File.variables['Latitude' ][l,:,0:self.x])
		f2 = squeeze(File.variables['Longitude'][l,:,0:self.x])
		f3 = squeeze(File.variables['Pressure' ][l,:,0:self.x])
		f  = array([f1,f2,f3])
		f  = rollaxis(f,1,0)
		f  = rollaxis(f,2,1)
	else:
		print '(%s,%s,%s,%s) not in dataset!' % (Year,Month,Day,Hour)
		sys.exit()
	File.close()
        return f

    def snappoints(self,Year,Month,Day,Hour):
	# Get all trajectory points at given date
	date0    = (Year,Month,Day,Hour)
	F = []
	for k in range(len(self.steps)):
		date = self.getDate(self.getHours(*date0) - self.steps[k])
		F.append(self.snapshot(*date)[:,k,:])
	F = array(F).reshape((-1,3))
	return F

    def getDataSnaps(self,Year=None,Month=None,Day=None,Hour=None,Season=None,datelist=None):
	if datelist == None:
		datelist = self.getDateList(Year,Month,Day,Hour,Season)
	f = []	
	for i in datelist:		
		f.append(self.snapshot(*i))
	f = squeeze(array(f))
	return f

    def line(self,p1, p2):
	A = (p1[1] - p2[1])
	B = (p2[0] - p1[0])
	C = (p1[0]*p2[1] - p2[0]*p1[1])
	return A, B, -C

    def intersection(self,L1, L2):
	D  = L1[0] * L2[1] - L1[1] * L2[0]
	Dx = L1[2] * L2[1] - L1[1] * L2[2]
	Dy = L1[0] * L2[2] - L1[2] * L2[0]
	x = Dx / D
	y = Dy / D
	return x,y

    def interpolateTrackToLat(self,lons,lats,blat):
	if lats.min() <= blat:
		x           = np.where(lats<=blat)[0][0]
		lat11,lon11 = lats[x-1],lons[x-1]
		lat12,lon12 = lats[x],lons[x]
		j11,k11     = self.m(lon11,lat11)
		j12,k12     = self.m(lon12,lat12)
		j21,k21     = self.m(-180,blat)
		j22,k22     = self.m(180,blat)
		L1,L2       = self.line([j11,k11],[j12,k12]),self.line([j21,k21],[j22,k22])
		j0,k0       = self.intersection(L1, L2)
		P           = self.m(j0,k0,inverse=True)
		dx          = np.sqrt((j11-j0)**2 + (k11-k0)**2)/np.sqrt((j11-j12)**2 + (k11-k12)**2)
		return P,x,1-dx
	else:
		return -999,-999,-999

    def interpTrack(self,xs,n):
	xold = np.linspace(0,1,len(xs))
	xnew = np.linspace(0,1,len(xs)*(n+1) - n)
	f  = interpolate.interp1d(xold,xs)
	xs = f(xnew)
	return xs

    def crossing(self,Year=1998,Month=1):
	fname = '../Tracks/back1_crossings/%s_%s_%02d.nc' % ('cross',Year,Month)
	print fname
	if not os.path.isfile(fname):
		# Create array of all crossings of 70N, meridional moisture flux and time since initialisation
		proj       = LambertProjector(boundinglat=80,resolution=200.)
		ds         = reDataServer(Field='pw',LevType='surface_analysis')
		vq,Dates   = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1979-2016.30-1000hPa.70N.Annual.p')
		F          = self.getDataSnaps(Year=Year,Month=Month)
		datelist   = self.getDateList(Year=Year,Month=Month)
		hourlist   = [self.getHours(*date) for date in datelist]
		fluxlons   = np.array(range(0,180+1,1) + range(-179,0,1))
		inds       = np.arange(-4,4+1,1)
		lat,lon,p  = F[:,:,:,0],F[:,:,:,1],F[:,:,:,2]
		lon0s      = lon[0,:,0]
		x,y        = proj.m(lon,lat)
		n0,n1,n2   = lat.shape
		flux,steps = np.zeros((n0,n1)),np.zeros((n0,n1))
		loncross   = np.zeros((n0,n1))
		latstart   = np.zeros((n0,n1))
		Pcross     = np.zeros((n0,n1))
		Field      = np.zeros((n0,n1))
		for t in range(len(lat)):
			date = datelist[t]
			snap = proj(ds.snapshot(*date),ds.lon,ds.lat)
			for i in range(len(lat[t])):
				P,xii,dx  = self.interpolateTrackToLat(lon[t,i],lat[t,i],70)	
				ix,iy     = np.argmin((proj.x[0,:]-x[t,i,0])**2),np.argmin((proj.y[:,0]-y[t,i,0])**2)
				if xii != -999:	
					latst     = lat[t,i,0]
					field     = snap[iy,ix]
					lonc      = P[0]
					pc        = p[t,i,xii] + dx*(p[t,i,xii-1] - p[t,i,xii])
					stepc     = xii + 1 - dx
					ixi       = Dates.index(self.getDate(self.getHours(*date) - xii*6))
					lonx      = np.argmin((fluxlons - P[0])**2)
					loninds   = lonx+inds
					for il in range(len(loninds)):
					        if loninds[il]<0:    loninds[il] = loninds[il]+360
					        if loninds[il]>=360: loninds[il] = loninds[il]-360
					vqslice = vq[ixi,loninds] + dx*(vq[ixi+1,loninds]-vq[ixi,loninds])
					vqslice = vqslice.mean()
				else:
					latst   = -9999
					field   = -9999
					lonc    = -9999
					pc      = -9999
					stepc   = -9999
					vqslice = -9999
				latstart[t,i] = latst
				Field[t,i]    = field
				loncross[t,i] = lonc
				Pcross[t,i]   = pc
				steps[t,i]    = stepc
				flux[t,i]     = vqslice
		latstart = np.ma.masked_where(latstart==-9999,latstart)
		Field    = np.ma.masked_where(Field==-9999,Field)
		loncross = np.ma.masked_where(loncross==-9999,loncross)
		Pcross   = np.ma.masked_where(Pcross==-9999,Pcross)
		steps    = np.ma.masked_where(steps==-9999,steps)
		flux     = np.ma.masked_where(flux==-9999,flux)
		crossNetCDF(hourlist,lon0s,latstart,Field,loncross,flux,steps,fname)
	else:
		print 'File %s already exists!' % (fname)
	
    def trackField1(self,Field,LevType,b1,b2,YearRange=(1981,2016),blat=85,plot=True,LonRange=(30,90)):
	fname = 'tracked/%s_vq70_%s-%s_%s-%skgm2.%s.p' % (Field,YearRange[0],YearRange[1],b1,b2,blat)
	print fname
	if not os.path.isfile(fname):
		# Dates
		dates    = unpick('../WarmArctic/dates.warm.1000.p')
		vq,Dates = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1979-2016.30-1000hPa.70N.Annual.p')
		fluxlons = np.array(range(0,180+1,1) + range(-179,0,1))
		inds    = np.arange(-4,4+1,1)
		#inds     = np.arange(0,0+1,1)
		# DataServers
		ds   = reDataServer(Field=Field,LevType=LevType)
		dsv  = reDataServer(Field='V',LevType='plev',LatRange=(70,70))
		dsq  = reDataServer(Field='q',LevType='plev',LatRange=(70,70))
		proj = LambertProjector(boundinglat=80,resolution=200.)
		# Get data
		Years     = range(YearRange[0],YearRange[1]+1,1)
		f70,nn0   = np.zeros(360),np.zeros(360)
		vplev     = np.zeros((16,360))
		nnplev    = np.zeros((16,360))
		pps,llons = [],[]
		datelist  = [ds.getDateList(Year=year,Season='DJF') for year in Years]
		datelist  = [datelist[i][j] for i in range(len(datelist)) for j in range(len(datelist[i]))]
		#datelist  = ds.getDateList(2001,Season='NDJFM')
		Flux,S,PP  = [],[],[]
		for date in datelist:
			print date
			s         = proj(ds.snapshot(*date),ds.lon,ds.lat)
			F         = self.snapshot(*date)
			lat,lon,p = F[:,:,0],F[:,:,1],F[:,:,2]
			x,y       = proj.m(lon,lat)
			for i in range(len(x)):
				ix,iy = np.argmin((proj.x[0,:]-x[i][0])**2),np.argmin((proj.y[:,0]-y[i][0])**2)
				if (b1 < s[iy,ix] < b2) and (blat<lat[i,0]):
					P,xii,dx = self.interpolateTrackToLat(lon[i],lat[i],70)
					if xii != -999 and (LonRange[0] < P[0] <= LonRange[1]):
						ixi     = Dates.index(ds.getDate(ds.getHours(*date) - xii*6))
						lonx    = np.argmin((fluxlons - P[0])**2)
						loninds = lonx+inds
						for il in range(len(loninds)):
							if loninds[il]<0:    loninds[il] = loninds[il]+360
							if loninds[il]>=360: loninds[il] = loninds[il]-360
						vqslice = vq[ixi,loninds] + dx*(vq[ixi+1,loninds]-vq[ixi,loninds])	
						#if vqslice.mean() > 0:	
							#vplev[:,loninds]  += (dsv.snapshot(*Dates[ixi+1])*dsq.snapshot(*Dates[ixi+1]))[:,loninds]
							#nnplev[:,loninds] += 1
							#pps.append(p[i,xii])
							#llon0 = lon[i,xii]
							#if llon0 < 0: llon0 = llon0 + 360
							#llons.append(llon0)
						f70[loninds]  += vqslice
						nn0[loninds]  += 1
						Flux.append(vqslice.mean())
						S.append(s[iy,ix])
						PP.append(p[i,:])
		PM,PS    = np.array(PP).mean(axis=0),np.array(PP).std(axis=0)
		f70,nn0  = np.ma.masked_where(nn0==0,f70),np.ma.masked_where(nn0==0,nn0)
		f70      = f70/nn0
		weighted = (f70*nn0).sum()/(nn0.sum())
		Ntot     = nn0.sum()
		#pl.figure(1)
		#pl.plot(S,Flux,'k.',alpha=0.15)
		#pl.figure(2)
		#days = np.arange(-10,0+0.25,0.25)
		#pl.plot(days,PM[::-1],'k')
		#pl.ylim(1000,700)
		#pl.plot(days,PM[::-1]-PS[::-1],'k--')
		#pl.plot(days,PM[::-1]+PS[::-1],'k--')
		#pl.ylim(1000,700)
		#pl.show()
                #vplev,nnplev = np.ma.masked_where(nnplev==0,vplev),np.ma.masked_where(nnplev==0,nnplev)
                #vplev        = vplev/nnplev
		#cseq         = np.arange(-0.012,0.012+0.002,0.002)
		#pl.contourf(dsv.lon,dsv.lev,vplev,cseq,cmap=pl.cm.RdBu_r,extend='both')
		#pl.plot(llons,pps,'k.',alpha=0.15)
		#pl.colorbar()
		#pl.savefig('figs/neg.test.pdf',format='pdf')
		#pl.ylim(1000,30)
		#pl.show()
		toPick([f70,nn0,weighted,Ntot],fname)
	else:
		f70,nn0,weighted,Ntot = unpick(fname)
	# Plot
	if plot:
		fig,ax1 = pl.subplots(num=2)
		ax2     = ax1.twinx()
		ax1.plot(range(360),f70,'k')
		ax1.plot([0,359],[0,0],'k--')
		ax2.plot(range(360),nn0,'r')
		ax1.set_ylim(-100,600)
		ax2.set_ylim(0,1200)
		ax1.set_ylabel('Mean meridional moisture flux')
		ax2.set_ylabel('Sample number')
		pl.xlabel('Longitude')
		pl.xlim(0,359)
		pl.title('DJF %s to %s: %s' % (YearRange[0],YearRange[1],weighted))
		pl.savefig('figs/characflux/%s_vq70_%s-%s_%s-%skgm2.pdf' % (Field,YearRange[0],YearRange[1],b1,b2), format='pdf')
		pl.close()
	return f70,nn0,weighted,Ntot

    def trackField2(self,Field,LevType,b1,b2,blat,LonRange):
        # Dates
        dates    = unpick('../WarmArctic/dates.warm.1000.p')
	days     = np.arange(-5,0+0.25,0.25)
        vq,Dates = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1979-2016.30-1000hPa.70N.Annual.p')
        # Flux axes
        fluxlons = np.array(range(0,180+1,1) + range(-179,0,1))
        inds    = np.arange(-4,4+1,1)
	#inds     = np.arange(0,0+1,1)
        # DataServers
        ds   = reDataServer(Field=Field,LevType=LevType)
        proj = LambertProjector(boundinglat=45,resolution=200.)
        # Get data 
        datelist = ds.getDateList(Year=2000,Month=1)
	#datelist = [ds.getDate(ds.getHours(*dates[8]) + k*24) for k in days]
        for date in datelist:
                s         = proj(ds.snapshot(*date),ds.lon,ds.lat)
                f         = interpolate.interp2d(proj.x, proj.y, s, kind='linear')
                F         = self.snapshot(*date)
                lat,lon,p = F[:,:,0],F[:,:,1],F[:,:,2]
                x,y       = proj.m(lon,lat)
                cseq      = np.arange(0,20+1,1)
                f70,nn0   = np.zeros(360),np.zeros(360)
                pl.figure(1)
                cf   = pl.contourf(proj.x,proj.y,s,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar = pl.colorbar()
                cbar.set_label('%s [%s]' % (ds.long_name,ds.units))
                pl.contour(proj.x,proj.y,s,levels=[-20],colors='k',linewidths=1.5)
                for i in range(len(x)):
                        if (b1 < f(x[i][0],y[i][0]) < b2) and (blat<lat[i,0]):
                                P,xii,dx = self.interpolateTrackToLat(lon[i],lat[i],70)
                                pl.plot(x[i][0],y[i][0],'k.')
                                pl.plot(x[i],y[i],'k',alpha=0.3)
                                if xii != -999 and (LonRange[0] < P[0] <= LonRange[1]):	
                                        ixi             = Dates.index(ds.getDate(ds.getHours(*date) - xii*6))
                                        lonx            = np.argmin((fluxlons - P[0])**2)
                                        loninds         = lonx+inds
                                        for il in range(len(loninds)):
                                                if loninds[il]<0:    loninds[il] = loninds[il]+360
                                                if loninds[il]>=360: loninds[il] = loninds[il]-360 
					f70[loninds] += vq[ixi,loninds] + dx*(vq[ixi+1,loninds]-vq[ixi,loninds])
                                        nn0[loninds] += 1
                                        xx,yy      = proj.m(P[0],P[1])
                                        pl.plot(xx,yy,'r.',alpha=0.6)
        	f70,nn0  = np.ma.masked_where(nn0==0,f70),np.ma.masked_where(nn0==0,nn0)
        	f70      = f70/nn0
        	weighted = (f70*nn0).sum()/(nn0.sum())
                proj.m.drawcoastlines()
                proj.m.drawparallels([70,80],latmax=90)
                pl.title('%s' % (str(date)))   
        	fig,ax1 = pl.subplots(num=2)
        	ax2     = ax1.twinx()
        	ax1.plot(range(360),f70,'k')
        	ax1.plot([0,359],[0,0],'k--')
        	ax2.plot(range(360),nn0,'r')
        	#ax1.set_ylim(-100,600) 
        	ax1.set_ylabel('Mean meridional moisture flux')
        	ax2.set_ylabel('Sample number')
        	pl.xlabel('Longitude')
        	pl.xlim(0,359)
        	pl.title('%s %s' % (str(date),weighted))  
        	pl.show()

    def interp2d(self,field,x,y,n,kind='cubic'):
	# New lat and lon grid (assumes axes are increasing)
	nx,ny   = int(len(x)*n),int(len(y)*n)
	new_x   = np.linspace(x[0],x[-1],nx)
	new_y   = np.linspace(y[0],y[-1],ny)
	# 2d interp function
	f = interpolate.interp2d(x,y,field,kind=kind,bounds_error=True)
	# Interpolate to new grid
	field = f(new_x,new_y)
	return field,new_x,new_y

    def boxBounds(self,m,LatRange,LonRange,n=35):
	if LonRange[1] < LonRange[0]:
		lons = np.append(np.linspace(LonRange[0],180,np.ceil(n/2.)),np.linspace(-180,LonRange[1],np.floor(n/2.)),axis=0)
	else:
		lons = np.linspace(LonRange[0],LonRange[1],n)
	lats = np.linspace(LatRange[0],LatRange[1],n)
	mer1 = m([LonRange[0] for i in range(n)],lats)
	mer2 = m([LonRange[1] for i in range(n)],lats)
	zon1 = m(lons,[LatRange[0] for i in range(n)])
	zon2 = m(lons,[LatRange[1] for i in range(n)])
	return [mer1,mer2,zon1,zon2]

    def density0(self,year0=1981,N=23):
	# Attributes
	proj        = LambertProjector(boundinglat=40,resolution=390.)
	projt       = LambertProjector(boundinglat=80,resolution=200.)
	xlim_,ylim0 = proj.m(0,50)
	xlim_,ylim1 = proj.m(180,50)
	xlim0,ylim_ = proj.m(-90,50)
	xlim1,ylim_ = proj.m(90,50)
	print proj.res
	# Axes
	xpole,ypole = proj.m(0,90)
	xyR2        = (proj.x-xpole)**2 + (proj.y-ypole)**2
	xyR2        = xyR2/1e13
	x,y         = proj.x[0,:],proj.y[:,0]
	# Masking of trakectories
	LatRange = [80,85]
	LonRange = [-30,90]
#	LonRange = [90,-150]
#	LonRange = [-150,-30]
	bounds   = self.boxBounds(proj.m,LatRange,LonRange,35)
	if LonRange[1] < LonRange[0]:
		inds = np.where((projt.lon.reshape(-1)<LonRange[1])|(projt.lon.reshape(-1)>LonRange[0])\
			       &(projt.lat.reshape(-1)<LatRange[1])&(projt.lat.reshape(-1)>LatRange[0])==True)[0]
	else:
                inds = np.where((projt.lon.reshape(-1)<LonRange[1])&(projt.lon.reshape(-1)>LonRange[0])\
                               &(projt.lat.reshape(-1)<LatRange[1])&(projt.lat.reshape(-1)>LatRange[0])==True)[0]
	# Data holders
	years1 = range(year0,year0+N/2,1)
	years2 = range(year0+N/2,year0+N,1)
	print years1,len(years1)
	print years2,len(years2)
	#years1 = range(1981,1998+1,1)
	#years2 = range(1999,2016+1,1)
	N1,N2  = np.array([np.zeros((proj.nx,proj.ny)) for i in years1]),np.array([np.zeros((proj.nx,proj.ny)) for i in years2])
	for t in range(len(years1)):
		year      = years1[t]
		print year
		s         = self.getDataSnaps(Year=year,Season='DJF')
		lat,lon,p = s[:,inds,:28,0],s[:,inds,:28,1],s[:,inds,:28,2]
		xs,ys     = proj.m(lon,lat)
		xs,ys     = xs.reshape(-1),ys.reshape(-1)
		for i in range(len(xs)):
			xi,yi    = np.argmin(abs(x-xs[i])),np.argmin(abs(y-ys[i]))
			N1[t,yi,xi] = N1[t,yi,xi] + 1
        for t in range(len(years2)):
		year      = years2[t]
                print year
                s         = self.getDataSnaps(Year=year,Season='DJF')
                lat,lon,p = s[:,inds,:28,0],s[:,inds,:28,1],s[:,inds,:28,2]
                xs,ys     = proj.m(lon,lat)
                xs,ys     = xs.reshape(-1),ys.reshape(-1)
                for i in range(len(xs)):
                        xi,yi    = np.argmin(abs(x-xs[i])),np.argmin(abs(y-ys[i]))
                        N2[t,yi,xi] = N2[t,yi,xi] + 1
	# Full time array for trend
	N0 = np.append(N1,N2,axis=0)
	Ntrend,pval = self.getTrend(years1+years2,N0)
	# Time series
	projlat1 = np.tile(proj.lat[np.newaxis,:,:],(len(years1),1,1))
	projlon1 = np.tile(proj.lon[np.newaxis,:,:],(len(years1),1,1))
	projlat2 = np.tile(proj.lat[np.newaxis,:,:],(len(years2),1,1))
	projlon2 = np.tile(proj.lon[np.newaxis,:,:],(len(years2),1,1))
	nn1 = np.ma.masked_where(projlat1>85,N1)
	nn1 = np.ma.masked_where(projlat1<70,nn1)
	nn1 = np.ma.masked_where(projlon1<20,nn1)
	nn1 = np.ma.masked_where(projlon1>90,nn1)
        nn2 = np.ma.masked_where(projlat2>85,N2)
        nn2 = np.ma.masked_where(projlat2<70,nn2)
        nn2 = np.ma.masked_where(projlon2<20,nn2)
        nn2 = np.ma.masked_where(projlon2>90,nn2)
	nn1 = nn1.sum(axis=1).sum(axis=1)/10000.
	nn2 = nn2.sum(axis=1).sum(axis=1)/10000.
	"""
	pl.plot(years1,nn1,'b',linewidth=1.5)
	pl.plot([years1[-1],years2[0]],[nn1[-1],nn2[0]],'r',linewidth=1.5)
	pl.plot(years2,nn2,'r',linewidth=1.5)
	pl.ylabel('Number density [10$^{4}$]')
	pl.xlabel('Year')
	pl.grid()
	pl.savefig('figs/dshed_timeseries.pdf',format='pdf')
	pl.show()
	"""
	# Format data
	N1,N2        = N1.sum(axis=0),N2.sum(axis=0)
	N1,N2        = N1[1:-1,1:-1]/(120.*len(years1)),N2[1:-1,1:-1]/(120.*len(years2))
	Ntrend       = Ntrend[1:-1,1:-1]/120.
	pval         = pval[1:-1,1:-1]
	x,y          = x[1:-1],y[1:-1]
	N1,xx,yy     = self.interp2d(N1,x,y,6,kind='cubic')
	N2,xx,yy     = self.interp2d(N2,x,y,6,kind='cubic')
	#Ntrend,xx,yy = self.interp2d(Ntrend,x,y,6,kind='cubic')
	#pval,xx,yy   = self.interp2d(pval,x,y,6,kind='cubic')
	# Plot
	cseqa,cseqf = np.arange(-3,3+0.5,0.5),np.arange(0,80+5,5)
	pl.figure(1)
	cf   = pl.contourf(xx,yy,N2-N1,cseqa,cmap=pl.cm.RdBu_r,extend='both')
	cbar = pl.colorbar(cf)
	cbar.set_label(r'Number density {400$\times$400 km$^{2}$ day}$^{-1}$')
	drawCoastlinesNoRivers(proj.m)
	proj.m.drawparallels([70,80],latmax=90)
	pl.xlim(xlim0,xlim1)
	pl.ylim(ylim0,ylim1)
	#pl.xlim(x[1],x[-2])
	#pl.ylim(y[1],y[-2])
	pl.title('%s to %s minus %s to %s' % (years2[0],years2[-1],years1[0],years1[-1]))
	pl.savefig('figs/dshed.pdf',format='pdf')
	pl.clf()
	pl.figure(2)
        cf   = pl.contourf(xx,yy,N1,cseqf,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Number density {400$\times$400 km$^{2}$ day}$^{-1}$')
        drawCoastlinesNoRivers(proj.m)
        proj.m.drawparallels([70,80],latmax=90)
        pl.xlim(x[1],x[-2])
        pl.ylim(y[1],y[-2])
	pl.title('%s to %s' % (years1[0],years1[-1]))
	pl.savefig('figs/shed.%s-%s.pdf' % (years1[0],years1[-1]),format='pdf')
	pl.clf()
        pl.figure(3)
        cf   = pl.contourf(xx,yy,N2,cseqf,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Number density {400$\times$400 km$^{2}$ day}$^{-1}$')
        drawCoastlinesNoRivers(proj.m)
        proj.m.drawparallels([70,80],latmax=90)
        pl.xlim(x[1],x[-2])
        pl.ylim(y[1],y[-2])
	pl.title('%s to %s' % (years2[0],years2[-1]))
	pl.savefig('figs/shed.%s-%s.pdf' % (years2[0],years2[-1]),format='pdf')
	pl.clf()
        pl.figure(4)
	mx,my = np.ma.masked_where(pval>0.1,proj.x[1:-1,1:-1]),np.ma.masked_where(pval>0.1,proj.y[1:-1,1:-1])
        cf    = pl.contourf(x,y,Ntrend,cseqa,cmap=pl.cm.RdBu_r,extend='both')
        cbar  = pl.colorbar(cf)
        cbar.set_label(r'Number density {400$\times$400 km$^{2}$ day}$^{-1}$ decade$^{-1}$')
	pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
        drawCoastlinesNoRivers(proj.m)
        proj.m.drawparallels([70,80],latmax=90)
        for xi,yi in bounds:
        	pl.plot(xi,yi,'b',linewidth=1.5,alpha=0.6)
        pl.xlim(x[1],x[-2])
        pl.ylim(y[1],y[-2])
        pl.title('Trend %s to %s' % (years1[0],years2[-1]))
        pl.savefig('figs/shed_trend.%s.%syears.%s-%sN.%s-%sE.pdf' % (year0,N,LatRange[0],LatRange[1],LonRange[0],LonRange[1]),format='pdf')
	pl.clf()
	#pl.show()

    def tracksInBox(self,YearRange=(1981,2016)):
	maskc     = 9.96921e36/2.
	Dir       = '../Tracks/back1_crossings'
	LATS,LONS = [],[]
	PWS,FS    = [],[]
	Years     = range(YearRange[0],YearRange[1]+1,1)
	for Year in Years:
		s     = self.getDataSnaps(Year=Year,Season='DJF')
		fname = '%s/cross_%s_%02d.nc' % (Dir,Year-1,12)
		print fname
		File  = Dataset(fname,'r')
		Flux  = File.variables['Flux'][:]
		Step  = File.variables['Step'][:]
		PW    = File.variables['PW'][:]
		for Month in [1,2,12]:	
			fname = '%s/cross_%s_%02d.nc' % (Dir,Year,Month)
			print fname
			File  = Dataset(fname,'r')	
			Flux  = np.append(Flux,File.variables['Flux'][:],axis=0)
			Step  = np.append(Step,File.variables['Step'][:],axis=0)
			PW    = np.append(PW,File.variables['PW'][:],axis=0)
		lat,lon,p = s[:,:,:,0],s[:,:,:,1],s[:,:,:,2]
		for t in range(len(lat)):
			for i in range(len(lat[t])):
				lats,lons = lat[t,i,:],lon[t,i,:]
				if ((lats<85)&(lats>70)==True).any() and ((lons<80)&(lons>45)==True).any() and (0<=Step[t,i]<=41):
					LATS.append(lats)
					LONS.append(lons)
					if Flux[t,i]<=maskc:
						FS.append(Flux[t,i])
						PWS.append(PW[t,i])
	LATS,LONS,FS = np.array(LATS),np.array(LONS),np.array(FS)
	return FS,PWS

    def histo(self):
	FS1,PWS1 = self.tracksInBox(YearRange=(1981,1998))
	FS2,PWS2 = self.tracksInBox(YearRange=(1999,2016))
	print len(FS1),len(FS2)
	slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(FS1,PWS1)
	slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(FS2,PWS2)
	xs    = np.linspace(-1000,2000,50)
	line1 = [x*slope1 + intercept1 for x in xs]
	line2 = [x*slope2 + intercept2 for x in xs]
	pl.plot(FS1,PWS1,'b.',alpha=0.1)
	pl.plot(FS2,PWS2,'r+',alpha=0.1)
	pl.plot(xs,line1,'k--',linewidth=1.2)
	pl.plot(xs,line2,'g--',linewidth=1.2)
	pl.grid()
	pl.show()
	"""
	h1,e1    = np.histogram(FS1,40,normed=True)
	h2,e2    = np.histogram(FS2,40,normed=True)
	pl.plot(e1[0:-1],h1,'b',linewidth=1.5,alpha=0.5)
	pl.plot(e2[0:-1],h2,'r',linewidth=1.5,alpha=0.5)
	pl.xlabel('Meridional moisture flux [Tg day$^{-}$ deg$^{-1}$]')
	pl.ylabel('Frequency')
	pl.grid()
	pl.show()
	"""

    def getTrend(self,years,field):
        # Field of shape (time)x(X)x(Y)
        n0,n1,n2 = field.shape
        trend    = np.zeros((n1,n2))
        pval     = np.zeros((n1,n2))
        for i in range(n1):
                for j in range(n2):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(years,field[:,i,j])
                        trend[i,j] = slope
                        pval[i,j]  = p_value
        return 10.*trend,pval

if __name__ == "__main__":

	import matplotlib.pyplot as pl
	from ReanalysisDataServer import DataServer as reDataServer
	from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
	from UnPickle import *
	from scipy import interpolate

	d = DataServer(Source='CAM4xCO2',Type='fwrd',steps=21)
	s = d.getDataSnaps(Year=1900,Month=1)
	print s.shape
	print s[0,:,0,0]
	print s[0,:,0,1]
	print s[0,:,0,2]
	
	"""
	Field   = 'fls'
	LevType = 'surface_forecast'
	blat    = 85
	plot    = True
	f701,nn01,weighted1,Ntot1,Flux1 = d2.trackField1(Field,LevType,0,6,YearRange=(1981,2016),blat=blat,plot=plot)
	f702,nn02,weighted2,Ntot2,Flux2 = d2.trackField1(Field,LevType,2,8,YearRange=(1981,2016),blat=blat,plot=plot)
	f703,nn03,weighted3,Ntot3,Flux3 = d2.trackField1(Field,LevType,3,6,YearRange=(1981,2016),blat=blat,plot=plot)

	bins  = np.linspace(-500,2000,125)
	h1,e1 = np.histogram(Flux1,bins)
	h2,e2 = np.histogram(Flux2,bins)
	h3,e3 = np.histogram(Flux3,bins)
	pl.plot(e1[0:-1],h1,'k',label='0-6 kg m$^{-2}$')
	pl.plot(e2[0:-1],h2,'b',label='2-6 kg m$^{-2}$')
	pl.plot(e3[0:-1],h3,'r',label='3-6 kg m$^{-2}$')
	lg = pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
	pl.show()
	"""

	"""
	Field       = 'pw'
	LevType     = 'surface_analysis'
	blat,plot   = 85,False
	#bounds,bins = zip(range(-90,10+10,10),range(-80,20+10,10)),np.arange(-85,15+10,10)
	bounds,bins = zip(range(0,8+1,1),range(1,9+1,1)),np.arange(0.5,8.5+1,1)
	W1,W2 = [],[]
	N1,N2 = [],[]
	for b1,b2 in bounds:
		f701,nn01,weighted1,Ntot1 = d2.trackField1(Field,LevType,b1,b2,YearRange=(1981,1998),blat=blat,plot=plot,LonRange=LonRange)
		f702,nn02,weighted2,Ntot2 = d2.trackField1(Field,LevType,b1,b2,YearRange=(1999,2016),blat=blat,plot=plot,LonRange=LonRange)		
		W1.append(weighted1)
		W2.append(weighted2)
		N1.append(Ntot1/1000.)
		N2.append(Ntot2/1000.)
	fig,ax1 = pl.subplots(num=1)
	ax2     = ax1.twinx()
	ax1.plot(bins,W1,'b',label='1981-1998')
	ax1.plot(bins,W2,'r',label='1999-2016')
	ax1.set_ylabel('Characteristic flux [Tg day$^{-1}$ deg$^{-1}$]')
        ax2.plot(bins,N1,'b--',label='1981-1998')
        ax2.plot(bins,N2,'r--',label='1999-2016')
	ax2.set_ylabel('Frequency [10$^{3}$]')
	ax1.set_xlabel('Precipitable water [kg m$^{-2}$]')
	#ax1.set_xlabel('Net longwave radiation [W m$^{-2}$]')
	ax1.grid()
	lg = ax1.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
	pl.savefig('figs/%s_characteristic_vq70.pdf' % (Field),format='pdf')
	pl.show()
	"""
