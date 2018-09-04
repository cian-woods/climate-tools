import sys
sys.path.insert(0, '/home/cian/scripts')
from LambertProjector import *
from ReanalysisDataServer import *
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from mpl_toolkits.basemap import Basemap as Basemap
from mpl_toolkits import basemap
from scipy import stats,interpolate
from UnPickle import *
from toPick import *
from stipling import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers

import matplotlib.pyplot as pl
import numpy as np
import sys,glob

SMs = {}
SMs['DJF']       = [12,1,2]
SMs['NDJF']      = [11,12,1,2]
SMs['DJ']        = [12,1]
SMs['JJA']       = [6,7,8]
SMs['SON']       = [9,10,11]
SMs['MAM']       = [3,4,5]
SMs['Annual']    = range(1,12+1,1)
SMs['Annual_wc'] = range(7,12+1,1) + range(1,6+1,1)

class Cyclone:
	def __init__(	self,
			lats,
			lons,
			times,
			slp	):

		self.lats  = lats
		self.lons  = lons
		self.times = times
		self.slp   = slp

class CycloneServer:

	def __init__(	self,
			method     = 'M13',
			hemisphere =  'NH',
			blat       =   55 ,
			resolution =  350.	):

		# File directory and pathname
		if method != 'MSS':
			self.Dir = '/mnt/climstorage/obs/IMILAST/ERAInterim_1.5_%s_ST_MTEX' % (hemisphere)
		else:
			self.Dir = '/mnt/climstorage/cian/cyclones/'
		self.File       = glob.glob(self.Dir+'/*%s*.txt' % (method))[0].split('/')[-1]
                self.method     = method
                self.hemisphere = hemisphere
		# Date start and end
		self.datestart = (int(self.File.split('_')[4][0:4]),int(self.File.split('_')[4][4:6]),int(self.File.split('_')[4][6:]),0)
		self.dateend   = (int(self.File.split('_')[5][0:4]),int(self.File.split('_')[5][4:6]),int(self.File.split('_')[5][6:]),18)
		print 'Cyclones from %s to %s' % (str(self.datestart),str(self.dateend))
		# Make Cyclones
		self.Cyclones = self.readCyclones()
		# Basemap and DataServer
		self.m        = Basemap(projection='nplaea',boundinglat=blat,lon_0=0)
		self.proj     = LambertProjector(boundinglat=blat,resolution=resolution)
		self.reds     = DataServer(Field='slp',LevType='surface_analysis',Source='ERAInt')
		self.x,self.y = self.proj.x[0,:],self.proj.y[:,0]
		# Datelist of data
		self.datelist = [self.reds.getDate(hour) for hour in np.arange(self.reds.getHours(*self.datestart),self.reds.getHours(*self.dateend)+6,6)]
		self.hourlist = [self.reds.getHours(*date) for date in self.datelist]
		#self.datelist = self.datelist[:100]
		#self.hourlist = self.hourlist[:100]

        def detrend2d(self,years,field):
                n0,n1,n2 = field.shape
                m,c,p    = np.zeros((n1,n2)),np.zeros((n1,n2)),np.zeros((n1,n2))
                for i in range(n1):
                        for j in range(n2):
                                slope, intercept, r_value, p_value, std_err = stats.linregress(years,field[:,i,j])
                                m[i,j] = slope
                                c[i,j] = intercept
                                p[i,j] = p_value
                line  = np.array([m*year + c for year in years])
                field = field - line
                return field,m,p

	def ice(self):
		mds = MMDataServer(Field='ci',LevType='surface_analysis')
		s   = np.array([mds.getSeason(Year=year,Season='NDJFM') for year in range(1989,2010+1,1)]).mean(axis=0)
		s   = self.proj(s,mds.lon,mds.lat)
		pl.contour(self.x,self.y,s,levels=[0.15])
		#self.proj.m.drawcoastlines()
		drawCoastlinesNoRivers(self.proj.m)
		pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/sea-ice.pdf',format='pdf')
		pl.show()

	def readCyclones(self):
		print 'Reading cyclones from %s' % self.File
    		# read all text in file
    		lines = open('%s/%s' % (self.Dir,self.File),'r').readlines()
    		lines = [ line.split() for line in lines ]
    		# find indices corresponding to breaks between cyclones
    		N = len(lines)
    		ix = []
    		for i in range(N):
      			if lines[i][0] == '99': continue
      			if lines[i][0] == '90': ix.append(i)
    		ix.append(N)
    		# iterate over cyclones
    		Cyclones = []
    		for j in range(len(ix)-1):
      			i0 = ix[j]+1
      			i1 = ix[j+1]
      			lats = [ float(line[9]) for line in lines[i0:i1] ]
      			lons = [ float(line[8]) for line in lines[i0:i1] ]
      			times = [ (int(line[4]),int(line[5]),int(line[6]),int(line[7])) for line in lines[i0:i1] ]
      			slp = [ float(line[10]) for line in lines[i0:i1] ]
      			Cyclones.append(Cyclone(lats,lons,times,slp))
    		print 'Read %s cyclones' % len(Cyclones)
    		return Cyclones

	def cycloneSnapshot(self,date):
		cs,xs = [],[]
		for cyclone in self.Cyclones:
			if date in cyclone.times:	
				x = cyclone.times.index(date)
				cs.append(Cyclone([cyclone.lats[x]],[cyclone.lons[x]],[cyclone.times[x]],[cyclone.slp[x]]))
				xs.append(1.*x/(len(cyclone.times)-1))
		return cs,xs

	def cyclonesInDateList(self,Cyclones,datelist,lags=[-8,-2]):
		# Get all cyclones that existed within datelist and sector
		# Transform datelist into lag window
		dates,days  = [],np.arange(lags[0],lags[1]+0.25,0.25)
		for date0 in datelist:
			dates = dates + [self.reds.getDate(self.reds.getHours(*date0) + k*24) for k in days]
		# Only unique dates
		dates_ = []
		for date in dates:
			if dates_.count(date)==0: dates_.append(date)
		dates = dates_
		# Take only dates in Cyclone dataset range
		dates = [date for date in dates if (self.hourlist[0]<=self.reds.getHours(*date)<=self.hourlist[-1])]	
		C     = []
		for ii in range(len(Cyclones)):
			Cyclone   = Cyclones[ii]
			dateunion = list(set(Cyclone.times)&set(dates))
			if (dateunion != []):
				C.append(Cyclone)
		return len(dates),C

	def cyclonesInSector(self,Cyclones,sector=[[0,90],[60,70]]):
		# Find the cyclones that enetered the area
		C = []
		for Cyclone in Cyclones:
			lons,lats = Cyclone.lons,Cyclone.lats
			for i in range(len(lons)):
				if (sector[0][0] <= lons[i] <= sector[0][1]) and (sector[1][0] <= lats[i] <= sector[1][1]):
					C.append(Cyclone)
					break
		return C

	def cyclonesInDateListAndSector(self,Cyclones,datelist,lags=[-8,-2],sector=[[0,90],[60,70]],when='all'):
                # Get all cyclones that existed within datelist and sector simultaneously
                # Transform datelist into lag window
                dates,days  = [],np.arange(lags[0],lags[1]+0.25,0.25)
                for date0 in datelist:
                        dates = dates + [self.reds.getDate(self.reds.getHours(*date0) + k*24) for k in days]
                # Only unique dates
                dates_ = []
                for date in dates:
                        if dates_.count(date)==0: dates_.append(date)
                dates = dates_
                # Take only dates in Cyclone dataset range
                dates = [date for date in dates if (self.hourlist[0]<=self.reds.getHours(*date)<=self.hourlist[-1])]
                C     = []
                for ii in range(len(Cyclones)):
                        Cyclone   = Cyclones[ii]
			# Select cyclones within in daterange over any, genesis or lysis timesteps
			if when == 'all': times = Cyclone.times
			if when == 'gen': times = [Cyclone.times[0]]
			if when == 'lys': times = [Cyclone.times[-1]]
                        dateunion = list(set(times)&set(dates))
                        if (dateunion != []):
				for date in dateunion:
					x       = times.index(date)	
					lon,lat = Cyclone.lons[x],Cyclone.lats[x]
					if (sector[0][0] <= lon <= sector[0][1]) and (sector[1][0] <= lat <= sector[1][1]):
                                		C.append(Cyclone)
						break
		return C

	def getMonthsCyclones(self,Months=[12,1,2]):
		# returns list of N cyclones

		Cyclones = []
		print 'Extracting cyclones for the months %s...' % (Months)
		for Cyclone in self.Cyclones:
			time0,time1 = Cyclone.times[0][1],Cyclone.times[-1][1]
			if (time0 in Months):# or (time1 in Months):	
				Cyclones.append(Cyclone)
		print 'Extracted %s cyclones' % (len(Cyclones))
		return Cyclones

	def getMonthListCyclones(self,hourlist,FC,Months,YearRange):
		years = range(YearRange[0],YearRange[1]+1,1)
		cross = False
		if Months[0]>Months[-1]: cross = True
		fc = {}
		for i in range(len(hourlist)):
			date  = self.reds.getDate(hourlist[i])
			year  = date[0]
			month = date[1]
			if (cross == True) and (month>=7): year = year + 1
			if (month in Months) and (year in years):
				try:
					fc[year].append(FC[i])
				except:
					fc[year] = []
					fc[year].append(FC[i])
		return fc.keys(),[np.array(fc[i]) for i in fc]

	def getSeasonCyclones(self,YearRange=(1990,2009),Season='DJF'):
		# returns list of (season)x(N) cyclones

		if 'DJ' in Season:
			i0,i1 = Season.index('D'),len(Season)-1
			m0,m1 = 12-i0,i1-i0
			cross = True
		else:
			# cc is for decdiding if the Season goes in same year or 
			# next year eg. Oct 1990 is in OrRd of 1991
			if Season ==  'ON': ms,cc =  [10,11],1
			if Season == 'JJA': ms,cc =  [6,7,8],0
			if Season == 'JFM': ms,cc =  [1,2,3],0
			if Season == 'OND': ms,cc =  [10,11,12],1
			cross = False

		# years in YearRange and Cyclone holder setup
		years    = np.arange(YearRange[0],YearRange[1]+1,1)
		Cyclones = {}
		for year in years:
			Cyclones[year] = []

		# Seasons where year is crosses i.e. containing DJ
		if cross == True:
			for Cyclone in self.Cyclones:
				time0 = Cyclone.times[0]
				if time0[1] <= m1:
					try:
						Cyclones[time0[0]].append(Cyclone)
					except:
						pass
				elif time0[1] >= m0:
					try:
						Cyclones[time0[0]+1].append(Cyclone)
					except:
						pass
		# Seasons where all months fall in same year
		elif cross == False:
			for Cyclone in self.Cyclones:
				time0 = Cyclone.times[0]
				if time0[1] in ms:
					try:
						Cyclones[time0[0]+cc].append(Cyclone)
					except:
						pass

		N = [len(Cyclones[year]) for year in Cyclones.keys()]	
		return Cyclones

	def bearing(self,lon1,lat1,lon2,lat2):

		lon1,lat1,lon2,lat2 = lon1*np.pi/180.,lat1*np.pi/180.,lon2*np.pi/180.,lat2*np.pi/180.
		b = np.arctan2( np.sin(lon2-lon1) * np.cos(lat2), np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
		return b*180/np.pi

	def haversine(self,lon1,lat1,lon2,lat2):
    		"""
    		Calculate the great circle distance (in km) between two points 
    		on the earth
    		"""
    		lon1=lon1*np.pi/180.
    		lat1=lat1*np.pi/180.
    		lon2=lon2*np.pi/180.
    		lat2=lat2*np.pi/180.
    		# haversine formula 
    		dlon = lon2 - lon1 
    		dlat = lat2 - lat1 
    		a = np.sin(dlat/2.)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.)**2
    		c = 2 * np.arcsin(np.sqrt(a)) 
    		d = 6367. * c
    		return d

	def getDensity(self,Cyclones,res=3,r=564):
		'''
  		Get counts of features and cyclone 
  		r is search radius (km)
  		res is resolution of grid (degrees)
  		'''
  		# set up lat/lon grid
  		nlon    = 360/res
  		nlat    = 180/res
  		lon     = np.arange(nlon)*res - 180.
  		lat     = np.arange(nlat+1)*res - 90.
  		lon,lat = np.meshgrid(lon,lat)
  		# build up counts of features and tracks
  		FeatureCounts = np.zeros(lon.shape)*0.
  		TrackCounts   = np.zeros(lon.shape)*0.
		GenCounts     = np.zeros(lon.shape)*0.
		LysCounts     = np.zeros(lon.shape)*0.
  		cc            = 0
  		Ncyclones     = len(Cyclones)
  		for Cyclone in Cyclones:
    			#if cc%1000 == 0: print 'At cyclone %s of %s' % (cc,Ncyclones)
    			count    = np.zeros(lon.shape)*0.
			gencount = np.where(self.haversine(Cyclone.lons[0],Cyclone.lats[0],lon,lat) < r, 1, 0)
			lyscount = np.where(self.haversine(Cyclone.lons[-1],Cyclone.lats[-1],lon,lat) < r, 1, 0)
    			for i in range(len(Cyclone.lons)):
      				count += np.where(self.haversine(Cyclone.lons[i],Cyclone.lats[i],lon,lat) < r, 1, 0)
    			FeatureCounts += count
    			TrackCounts   += np.where(count > 0, 1, 0)
			GenCounts     += gencount
			LysCounts     += lyscount
    			cc = cc + 1
  		return lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts

        def getDensity0(self,lon,lat,Cyclones,r=564):
                '''
                Get counts of features and cyclone 
                r is search radius (km)
                res is resolution of grid (degrees)
                '''
                # build up counts of features and tracks
                Counts = np.zeros(lon.shape)*0. 
                for Cyclone in Cyclones: 
                        Counts += np.where(self.haversine(Cyclone.lons[0],Cyclone.lats[0],lon,lat) < r, 1, 0)
                return Counts

	def getClim(self,Season='NDJFM',Months=[11,12,1,2,3]):
		fname         = 'cycloneclims/seasons/%s.%s.p' % (self.method,Season)
		if not os.path.isfile(fname):
			Cyclones = self.getMonthsCyclones(Months)
			N        = len([date for date in self.datelist if date[1] in Months])
			lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts = self.getDensity(Cyclones)
			toPick([N,lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts],fname)
		else:
			N,lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts = unpick(fname)
		return 1.*N,lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts

	def getDensityN(self,Cyclones,res=3,r=564):
		# Same as getDensity but for Cyclones organised
		# by (year)x(N) from getSeasonCyclones

		FC,TC,GC,LC = [],[],[],[]
		for year in Cyclones.keys():
			lon,lat,fc,tc,gc,lc = self.getDensity(Cyclones[year],res=res,r=r)
			FC.append(fc)
			TC.append(tc)
			GC.append(gc)
			LC.append(lc)
		return lon,lat,np.array(FC),np.array(TC),np.array(GC),np.array(LC),Cyclones.keys()

	def getCycloneTracks(self,Cyclones):
		# Extract the track coordinates
		T = []
		for Cyclone in Cyclones:
			lats,lons,times,slp = Cyclone.lats,Cyclone.lons,Cyclone.times,Cyclone.slp
			T.append([lats,lons,times,slp])
		return T

	def getTrend(self,FC,TC,GC,LC,years):
		# Takes getDensityN counts as input

		FCtrend,TCtrend,GCtrend,LCtrend = np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:])
		FCp,TCp,GCp,LCp                 = np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:]),np.zeros(FC.shape[1:])
		for i in range(FC.shape[1]):
			# Trends
			for j in range(FC.shape[2]):
				mF,cF,rF,pF,errF = stats.linregress(years,FC[:,i,j])
				mT,cT,rT,pT,errT = stats.linregress(years,TC[:,i,j])
				mG,cG,rG,pG,errG = stats.linregress(years,GC[:,i,j])
				mL,cL,rL,pL,errL = stats.linregress(years,LC[:,i,j])
				FCtrend[i,j]     = mF
				TCtrend[i,j]     = mT
				GCtrend[i,j]     = mG
				LCtrend[i,j]     = mL
                                FCp[i,j]         = pF/2.
                                TCp[i,j]         = pT/2.
                                GCp[i,j]         = pG/2.
                                LCp[i,j]         = pL/2.
		return FCtrend,TCtrend,GCtrend,LCtrend,FCp,TCp,GCp,LCp

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

        def plotCounts1(self,TC1,TC2,savename='',num=0,stip1=None,stip2=None,cseq=20,clabel='',cmap=pl.cm.RdBu_r):
		TC1,x,y  = self.interp2d(TC1,self.x,self.y,6,kind='linear')
		TC2,x,y  = self.interp2d(TC2,self.x,self.y,6,kind='linear')
                # Plot
                fig  = pl.figure(figsize=(7,10))
                # Subplot 1
                pl.subplot(211)
                pl.figure(1)
                cf   = pl.contourf(x,y,TC1,cseq,cmap=cmap,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('%s' % (clabel))
                #self.proj.m.drawcoastlines()
		drawCoastlinesNoRivers(self.proj.m)
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title('%s' % (num))
                if stip1 != None:
                        pl.plot(stip1[0][::2,::2],stip1[1][::2,::2],'k.',alpha=0.5)
                # Subplot 2
                pl.subplot(212)
                cf   = pl.contourf(x,y,TC2,cseq,cmap=cmap,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('%s' % (clabel))
                #self.proj.m.drawcoastlines()
		drawCoastlinesNoRivers(self.proj.m)
                self.proj.m.drawparallels([70,80],latmax=90)
                if stip2 != None:
                       pl.plot(stip2[0][::2,::2],stip2[1][::2,::2],'k.',alpha=0.5)
                pl.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
                if savename != '':
                        pl.savefig(savename,format='pdf')
                        pl.close()
                else:
                        pl.show()

	def plotCounts2(self,TC1,TC2,savename='',num=0,stip1=None,stip2=None,cseq=20,clabel='',cmap=pl.cm.RdBu_r):
                TC1,x,y  = self.interp2d(TC1,self.x,self.y,6,kind='linear')
                TC2,x,y  = self.interp2d(TC2,self.x,self.y,6,kind='linear')
		# Plot
		pl.figure(1)
        	cf   = pl.contourf(x,y,TC1,cseq,cmap=cmap,extend='both')
		#cl1 = pl.contour(x,y,TC2,[0.50],colors='k',linewidths=2,alpha=0.75)
		#cl1 = pl.contour(x,y,TC2,[0.65,0.80,0.95,1.1],colors='k',linewidths=1.5,alpha=0.5)
		cl1  = pl.contour(x,y,TC2,[0.40],colors='k',linewidths=2,alpha=0.75)
		cl1  = pl.contour(x,y,TC2,[0.55,0.70,0.85,1.0,1.15],colors='k',linewidths=1.5,alpha=0.5)
		#pl.clabel(cl1,colors='k',fontsize=10)
        	cbar = pl.colorbar(cf)
		cbar.set_label('%s' % (clabel))
		drawCoastlinesNoRivers(self.proj.m,color='0.5',linewidth=0.7)
        	#self.proj.m.drawcoastlines(color='0.5',linewidth=0.7)
        	self.proj.m.drawparallels([70,80],latmax=90)
		pl.title('%s' % (num))
                if stip1 != None:
                        pl.plot(stip1[0][::2,::2],stip1[1][::2,::2],'k.',alpha=0.5)
		if savename != '':
			pl.savefig(savename,format='pdf')
			pl.close()
		else:
			pl.show()

	def plotTracks(self,Cyclones):

		for lats,lons,times,slp in Cyclones:
			j,k = self.m(lons,lats)
			pl.plot(j,k,'0.5',linewidth=0.5)
		self.m.drawcoastlines()
		self.m.drawparallels([70,80],latmax=90)
		pl.show()

	def makeClim(self):
                # set up lat/lon grid
		res     = 3.
                nlon    = 360/res
                nlat    = 180/res
                lon     = np.arange(nlon)*res - 180.
                lat     = np.arange(nlat+1)*res - 90.
                lon,lat = np.meshgrid(lon,lat)
		# FileName
		fname    = '/mnt/climstorage/cian/WarmArctic/cycloneclims/%s.nc' % (self.method)
		if not os.path.isfile(fname):
			FC,GC,LC = [],[],[]
			hourlist = self.hourlist
			for date in self.datelist:
				print date
				cs,xs = self.cycloneSnapshot(date)
				gcs   = [cs[i] for i in range(len(cs)) if xs[i]==0]
				lcs   = [cs[i] for i in range(len(cs)) if xs[i]==1]
				#lon,lat,FeatureCounts,TrackCounts,GenCounts,LysCounts = self.getDensity(cs)
				FeatureCounts = self.getDensity0(lon,lat,cs)
				GenesisCounts = self.getDensity0(lon,lat,gcs)
				LysisCounts   = self.getDensity0(lon,lat,lcs)
				FC.append(FeatureCounts)
				GC.append(GenesisCounts)
				LC.append(LysisCounts)
			FC,GC,LC = np.array(FC),np.array(GC),np.array(LC)
			self.saveFile(hourlist,lon,lat,FC,GC,LC,fname)
		else:
			File     = Dataset(fname,'r')
			hourlist = File.variables['time'][:]
			lon      = File.variables['lon'][:]
			lat      = File.variables['lat'][:]
			FC       = File.variables['FeatureCount'][:]
			GC       = File.variables['GenesisCount'][:]
			#LC       = File.variables['LysisCount'][:]	
			LC       = File.variables['GenesisCount'][:] 
		return hourlist,lon,lat,FC,GC,LC

	def saveFile(self,times,lon,lat,FC,GC,LC,FileName):
		if not os.path.isfile(FileName):
			# Create file
        		File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        		# Define some global attribs
			File.Conventions='COARDS'
			# Time is record dimension
			File.createDimension('time',len(times))
			var           = File.createVariable('time','d',('time',))
			var.long_name = 'Time'
			var.units     = 'hours since 1900-01-01 00:00:0.0'
			var[:]        = times
			# Horizontal axes
			File.createDimension('Y',lat.shape[0])
			File.createDimension('X',lat.shape[1])
			var           = File.createVariable('lat','f',('Y','X',))
			var.long_name = 'Latitude'
			var.units     = 'degrees_north'
			var[:]        = lat.astype('f')
			var           = File.createVariable('lon','f',('Y','X',))
			var.long_name = 'Longitude'
			var.units     = 'degrees_east'
			var[:]        = lon.astype('f')
			# Create Variables
			var           = File.createVariable('FeatureCount','f',('time','Y','X',))
			var.long_name = 'Cyclone feature counts within 564 km'
			var.units     = '0-1'
			var[:]        = FC
			var           = File.createVariable('GenesisCount','f',('time','Y','X',))
                        var.long_name = 'Cyclogenesis counts within 564 km'
                        var.units     = '0-1'
                        var[:]        = GC
			var           = File.createVariable('LysisCount','f',('time','Y','X',))
                        var.long_name = 'Cyclolysis counts within 564 km'
                        var.units     = '0-1'
                        var[:]        = LC
			# Close file
			File.close()

	def bivarPDF(self,x,y,xedges,yedges,norm=True):
		H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=None,normed=norm,weights=None)
		xedges,yedges   = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
		return xedges,yedges,H

	def cyclonesTrend(self,YearRange,Season,Months,method):
		# Get cyclone data
		hourlist,lon,lat,FC,GC,LC = self.makeClim()
		# Compile into months
		years,fc = self.getMonthListCyclones(hourlist,FC,Months,YearRange=YearRange)
		fc       = np.array([i.mean(axis=0) for i in fc])
		fc       = self.proj(fc,lon[0,:],lat[:,0])
		# Take trend
		fc,m,p   = self.detrend2d(years,fc)
		fc_sigma = fc.std(axis=0)
		# Plot
		pl.figure(1)
		cseq  = 13#np.arange(-24,24+4,4)
		sx,sy = np.ma.masked_where(p/2.>0.05,self.proj.x),np.ma.masked_where(p/2.>0.05,self.proj.y)
		cf    = pl.contourf(self.proj.x,self.proj.y,len(Months)*121.7*m*10,cseq,cmap=pl.cm.coolwarm,extend='both') # 121.7 = 6xhourly timesteps in an average Month i.e. 365*4/12 = 121.7
		cbar  = pl.colorbar(cf)
		cbar.set_label('Number denisty [%s$^{-1}$ decade$^{-1}$]' % (Season))
		pl.plot(sx,sy,'k.',alpha=0.5)
		self.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(self.proj.m,color='0.5',linewidth=0.7)
		pl.title(method)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/cyclones/N.trend.%s-%s.%s.%s.pdf' % (YearRange[0],YearRange[1],Season,method))
		pl.show()


if __name__ == '__main__':

        #methods  = ['M02','M03','M06','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        #for method in methods:
        #       CS = CycloneServer(method=method)
        #       d,lon,lat,FC,GC,LC = CS.makeClim()

        #CS = CycloneServer(method='M13',resolution=80.)
        #CS.ice()
        #h,lon,lat,FC,GC,LC = CS.makeClim()
        #print FC.shape,GC.shape,LC.shape
        #FC,GC,LC = CS.proj(FC.mean(axis=0),lon[0,:],lat[:,0]),CS.proj(GC.mean(axis=0),lon[0,:],lat[:,0]),CS.proj(LC.mean(axis=0),lon[0,:],lat[:,0])
        #CS.plotCounts(FC,LC,savename='',num=1,cseq=20)

	"""
        # Cyclone densities in sector and time range
        case        = 'warm'
        when        = 'all'
        lags,sector = [2.5,7.5],[[-180,180],[80,90]]
        datelist    = unpick('dates.%s.50.ERAInt.p' % (case))[:50]
        #methods    = ['M02','M03','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        methods     = ['M13']
        LAT0,LON0   = [],[]
        LAT1,LON1   = [],[]
        LAT2,LON2   = [],[]
        weights     = []
        FC,GC,LC    = [],[],[]
        Months      = [11,12,1,2,3]
        for method in methods:
                CS                             = CycloneServer(method=method)
                h,lon,lat,FCclim,GCclim,LCclim = CS.makeClim()
                xs                             = [i for i in range(len(CS.datelist)) if CS.datelist[i][1] in Months]
                Fclim,Gclim,Lclim              = FCclim[xs,:,:].mean(axis=0),GCclim[xs,:,:].mean(axis=0),LCclim[xs,:,:].mean(axis=0)
                Cyclones0                      = CS.cyclonesInDateListAndSector(CS.Cyclones,datelist,lags=lags,sector=sector,when=when)
                weights.append(1)
                LAT0.append([cyclone.lats[0]  for cyclone in Cyclones0])
                LON0.append([cyclone.lons[0]  for cyclone in Cyclones0])
                LAT1.append([cyclone.lats[-1] for cyclone in Cyclones0])
                LON1.append([cyclone.lons[-1] for cyclone in Cyclones0])
                lat_,lon_ = [],[]
                for cyclone in Cyclones0:
                        for i in range(len(cyclone.lats)):
                                lat_.append(cyclone.lats[i])
                                lon_.append(cyclone.lons[i])
                LAT2.append(lat_)
                LON2.append(lon_)
                FC.append(Fclim)
                GC.append(Gclim)
                LC.append(Lclim)
        FC,GC,LC = np.array(FC).mean(axis=0),np.array(GC).mean(axis=0),np.array(LC).mean(axis=0)
        weights = 1./np.array(weights)
        # set up lat/lon grid
        res           = 3.
        nlon          = 360/res
        nlat          = 180/res
        lon           = np.arange(nlon)*res - 180.
        lat           = np.arange(nlat+1)*res - 90.
        lon,lat       = np.meshgrid(lon,lat)
        C0,C1,C2      = [],[],[]
        for i in range(len(LAT0)):
                Count0 = np.zeros(lon.shape)
                for j in range(len(LAT0[i])):
                        Count0 += np.where(CS.haversine(LON0[i][j],LAT0[i][j],lon,lat) < 564, 1, 0)
                C0.append(Count0)
        for i in range(len(LAT1)):
                Count1 = np.zeros(lon.shape)
                for j in range(len(LAT1[i])):
                        Count1 += np.where(CS.haversine(LON1[i][j],LAT1[i][j],lon,lat) < 564, 1, 0)
                C1.append(Count1)
        for i in range(len(LAT2)):
                Count2 = np.zeros(lon.shape)
                for j in range(len(LAT2[i])):
                        Count2 += np.where(CS.haversine(LON2[i][j],LAT2[i][j],lon,lat) < 564, 1, 0)
                C2.append(Count2)
        Count0        = (np.array(C0)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count1        = (np.array(C1)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count2        = (np.array(C2)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count0,Count1 = CS.proj(Count0,lon[0,:],lat[:,0]),CS.proj(Count1,lon[0,:],lat[:,0])
        Count2        = CS.proj(Count2,lon[0,:],lat[:,0])
        FC,GC,LC      = CS.proj(FC,lon[0,:],lat[:,0]),CS.proj(GC,lon[0,:],lat[:,0]),CS.proj(LC,lon[0,:],lat[:,0])
        #cseqf,cseql,cmap = 14,14,pl.cm.OrRd
        cseqf,cseql,cmap = np.arange(0,8+1,1),np.arange(0,0.8+0.1,0.1),pl.cm.OrRd
        #cseqf,cseql,cmap = np.arange(0,0.35+0.025,0.025),np.arange(0,0.03+0.0025,0.0025),pl.cm.OrRd
        # Plot
        Count0,x,y = CS.interp2d(Count0,CS.x,CS.y,6,kind='linear')
        Count1,x,y = CS.interp2d(Count1,CS.x,CS.y,6,kind='linear')
        Count2,x,y = CS.interp2d(Count2,CS.x,CS.y,6,kind='linear')
        #Count0,Count1   = np.ma.masked_where(Count0<1,Count0),np.ma.masked_where(Count1<1,Count1)
        #Count2          = np.ma.masked_where(Count2<1,Count2)
        pl.figure(1)
        cf   = pl.contourf(x,y,Count0/30.,cseql,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count0.sum()))
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/cyclogensis_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.figure(2) 
        cf   = pl.contourf(x,y,Count1/30.,cseql,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count1.sum()))
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/cyclolysis_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.figure(3)
        cf   = pl.contourf(x,y,Count2/30.,cseqf,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count2.sum()))
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/feature_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.show()
	"""

	"""
        from Intrusions import *
        YearRange         = (1989,2009)
        years             = range(YearRange[0],YearRange[1]+1,1)
        IN                = Intrusions(OpenTraj=True)
        Nl                = 4*IN.intrusionLaggedDates(lags=[-10,0],case='warm',rankN=50,plot=True,YearRange=YearRange)
        print Nl.shape
        #methods          = ['M02','M03','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        methods           = ['M13']
        datelist1         = unpick('/home/cian/WarmArctic/dates.warm.50.ERAInt.p')[:50]
        datelist2         = unpick('/home/cian/WarmArctic/dates.cold.50.ERAInt.p')[:50]
        datelist1         = [date for date in datelist1 if date[0] in years]
        datelist2         = [date for date in datelist2 if date[0] in years]
        Months            = [11,12,1,2,3]
        FFC1,FFC2         = [],[]
        GGC1,GGC2         = [],[]
        LLC1,LLC2         = [],[]
        FCLIM,GCLIM,LCLIM = [],[],[]
        for method in methods:
                # CycloneServer and climatologies
                CS                             = CycloneServer(method=method)
                h,lon,lat,FCclim,GCclim,LCclim = CS.makeClim()
                xs                             = [i for i in range(len(CS.datelist)) if (CS.datelist[i][1] in Months) and (CS.datelist[i][0] in years)]
                Fclim,Gclim,Lclim              = FCclim[xs,:,:].mean(axis=0),GCclim[xs,:,:].mean(axis=0),LCclim[xs,:,:].mean(axis=0)
                # Use only date0s falling in cyclone range
                datelist10 = [date for date in datelist1 if (CS.hourlist[0]<=CS.reds.getHours(*date)<=CS.hourlist[-1])]
                datelist20 = [date for date in datelist2 if (CS.hourlist[0]<=CS.reds.getHours(*date)<=CS.hourlist[-1])]
                print len(datelist10),len(datelist20)
                FC1,FC2    = [],[]
                GC1,GC2    = [],[]
                LC1,LC2    = [],[]
                for date0 in datelist10:
                        x     = CS.datelist.index(date0)
                        fc    = FCclim[x-41+2:x+40+1,:,:]       # -9.75 to 10.0 days; len = 80/4 = 20 days
                        gc    = GCclim[x-41+2:x+40+1,:,:]
                        lc    = LCclim[x-41+2:x+40+1,:,:]
                        FC1.append(fc)
                        GC1.append(gc)
                        LC1.append(lc)
                for date0 in datelist20:
                        x     = CS.datelist.index(date0)
                        fc    = FCclim[x-41+2:x+40+1,:,:]
                        gc    = GCclim[x-41+2:x+40+1,:,:]
                        lc    = LCclim[x-41+2:x+40+1,:,:]
                        FC2.append(fc)
                        GC2.append(gc)
                        LC2.append(lc)
                FC1,FC2 = np.array(FC1).mean(axis=0)-Fclim[np.newaxis,:,:],np.array(FC2).mean(axis=0)-Fclim[np.newaxis,:,:]
                GC1,GC2 = np.array(GC1).mean(axis=0)-Gclim[np.newaxis,:,:],np.array(GC2).mean(axis=0)-Gclim[np.newaxis,:,:]
                LC1,LC2 = np.array(LC1).mean(axis=0)-Lclim[np.newaxis,:,:],np.array(LC2).mean(axis=0)-Lclim[np.newaxis,:,:]
                #FC1,FC2 = np.array(FC1).mean(axis=0),np.array(FC2).mean(axis=0)
                #GC1,GC2 = np.array(GC1).mean(axis=0),np.array(GC2).mean(axis=0)
                #LC1,LC2 = np.array(LC1).mean(axis=0),np.array(LC2).mean(axis=0)
                FFC1.append(FC1)
                FFC2.append(FC2)
                GGC1.append(GC1)
                GGC2.append(GC2)
                LLC1.append(LC1)
                LLC2.append(LC2)
                FCLIM.append(Fclim)
                GCLIM.append(Gclim)
                LCLIM.append(Lclim)
        #Nl               = CS.proj(Nl,lonNl[0,:],latNl[:,0])
        FFC1,FFC2         = np.array(FFC1),np.array(FFC2)
        GGC1,GGC2         = np.array(GGC1),np.array(GGC2)
        LLC1,LLC2         = np.array(LLC1),np.array(LLC2)
        FCLIM,GCLIM,LCLIM = np.array(FCLIM).mean(axis=0),np.array(GCLIM).mean(axis=0),np.array(LCLIM).mean(axis=0)
        FFC1,FFC2         = CS.proj(FFC1,lon[0,:],lat[:,0]),CS.proj(FFC2,lon[0,:],lat[:,0])
        GGC1,GGC2         = CS.proj(GGC1,lon[0,:],lat[:,0]),CS.proj(GGC2,lon[0,:],lat[:,0])
        LLC1,LLC2         = CS.proj(LLC1,lon[0,:],lat[:,0]),CS.proj(LLC2,lon[0,:],lat[:,0])
        FFC1m,FFC2m       = 605*FFC1.mean(axis=0),605*FFC2.mean(axis=0)
        GGC1m,GGC2m       = 605*GGC1.mean(axis=0),605*GGC2.mean(axis=0)
        LLC1m,LLC2m       = 605*LLC1.mean(axis=0),605*LLC2.mean(axis=0)
        #cseqf,cseql      = np.arange(-60,60+10,10),np.arange(-9,9+1.5,1.5)
        cseqf,cseql      = np.arange(-90,90+15,15),np.arange(-12,12+1.5,1.5)
        cmap              = pl.cm.RdBu_r


        segs  = [[-60, -40], [-72, -52], [ -51, -30], [-29, -9]]
        sdays = [[ -5,   0], [ -7,  -3], [-2.5, 2.5], [  3,  7]]
        for ii in range(len(segs)):
                # Mean composite between day -5 to 0
                CS.plotCounts1(FFC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),FFC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Feature/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseqf,clabel=r'Cyclone freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
                CS.plotCounts1(GGC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),GGC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Genesis/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseql,clabel=r'Genesis freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
                CS.plotCounts1(LLC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),LLC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Lysis/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseql,clabel=r'Lysis freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)


        # Lagged composites between day -10 to 0 over 2-day windows and mean composite between day -5 to 0
        for i in range(len(FFC1m)-8):
                CS.plotCounts2(FFC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Feature/%s.pdf' % (i+1),\
                                num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseqf,\
                               clabel=r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        CS.plotCounts2(FFC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Feature/mean_int.pdf',\
                        num='day -5 to 0',stip1=None,stip2=None,cseq=cseqf,clabel=r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        for i in range(len(GGC1m)-8):
                CS.plotCounts2(GGC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Genesis/%s.pdf' % (i+1),\
                                num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseql,
                               clabel=r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        CS.plotCounts2(GGC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Genesis/mean_int.pdf',\
                        num='day -5 to 0',stip1=None,stip2=None,cseq=cseql,clabel=r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        for i in range(len(LLC1m)-8):
                CS.plotCounts2(LLC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Lysis/%s.pdf' % (i+1),\
                                num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseql,\
                               clabel=r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        CS.plotCounts2(LLC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='/mnt/climstorage/cian/WarmArctic/figs/lagged/Lysis/mean_int.pdf',\
                        num='day -5 to 0',stip1=None,stip2=None,cseq=cseql,clabel=r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)

        # Plot CLIMs
        cseqf,cseql = np.arange(0,210+15,15),np.arange(0,18+1.5,1.5)
        #cseqf,cseql = np.arange(0,0.35+0.025,0.025),np.arange(0,0.03+0.0025,0.0025)
        FCLIM     = CS.proj(FCLIM,lon[0,:],lat[:,0])
        GCLIM     = CS.proj(GCLIM,lon[0,:],lat[:,0])
        LCLIM     = CS.proj(LCLIM,lon[0,:],lat[:,0])
        FCLIM,x,y = CS.interp2d(FCLIM,CS.x,CS.y,6,kind='linear')
        GCLIM,x,y = CS.interp2d(GCLIM,CS.x,CS.y,6,kind='linear')
        LCLIM,x,y = CS.interp2d(LCLIM,CS.x,CS.y,6,kind='linear')
        cf    = pl.contourf(x,y,605*FCLIM,cseqf,cmap=pl.cm.OrRd,extend='max')
        cbar  = pl.colorbar()
        cbar.set_label(r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80,85],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/IMILAST.fclim.pdf',format='pdf')
        pl.close()
        cf   = pl.contourf(x,y,605*GCLIM,cseql,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar()
        cbar.set_label(r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/IMILAST.gclim.pdf',format='pdf')
        pl.close()
        cf   = pl.contourf(x,y,605*LCLIM,cseql,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar()
        cbar.set_label(r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/IMILAST.lclim.pdf',format='pdf')
        pl.close()
	"""

	"""
	CS = CycloneServer('MSS','NH',40,200)
	hourlist,lon,lat,FC,GC,LC = CS.makeClim()
	"""

        #YearRange     = (int(sys.argv[1]),int(sys.argv[2]))
        #Season        = str(sys.argv[3])
	#Months        = SMs[Season]

	"""
	CS = CycloneServer('MSS','NH',40,200)
	CS.cyclonesTrend(YearRange,Season,Months,'MSS')
	"""

	"""
	CS = CycloneServer(     method     =  'M13',
                        	hemisphere =   'NH',
                        	blat       =     55,
                        	resolution =    200     )

	years,htmoor  = unpick('/mnt/climstorage/cian/BSO/hfTOT30d.bso.dtrend.%s-%s.%s.p' % (YearRange[0],YearRange[-1],Season))
	htmoor_sigma  = htmoor.std()/1e12

        # Wind data
        lag       = 0
        usds,vsds = MMDataServer(Field='U',LevRange=(1000,1000)),MMDataServer(Field='V',LevRange=(1000,1000))
        VS = CS.proj(np.array([vsds.getSeason(Year=year+lag,Season=Season).squeeze() for year in years]),vsds.lon,vsds.lat)
        US = CS.proj(np.array([usds.getSeason(Year=year+lag,Season=Season).squeeze() for year in years]),usds.lon,usds.lat)
        # Detrend
        us ,m_us ,p_us  = CS.detrend2d(years, US)
        vs ,m_vs ,p_vs  = CS.detrend2d(years, VS)
        # Regress
        mus,pus = np.zeros((CS.proj.nx,CS.proj.ny)),np.zeros((CS.proj.nx,CS.proj.ny))
	mvs,pvs = np.zeros((CS.proj.nx,CS.proj.ny)),np.zeros((CS.proj.nx,CS.proj.ny))
        for i in range(CS.proj.nx):
                for j in range(CS.proj.ny):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(htmoor,us[:,i,j])
                        mus[i,j],pus[i,j] = slope,p_value/2.
                        slope, intercept, r_value, p_value, std_err = stats.linregress(htmoor,vs[:,i,j])
                        mvs[i,j],pvs[i,j] = slope,p_value/2.

	methods     = ['MSS']#['M02','M03','M06','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
	mm,pp,MM,PP = [],[],[],[]
	MR,PR       = [],[]
	for method in methods:
		CS = CycloneServer(	method     = method,
					hemisphere =   'NH',
					blat       =     55,
					resolution =    200	)

		hourlist,lon,lat,FC,GC,LC = CS.makeClim()
		years,fc = CS.getMonthListCyclones(hourlist,FC,Months,YearRange=YearRange)
		fc       = np.array([i.mean(axis=0) for i in fc])
		fc       = CS.proj(fc,lon[0,:],lat[:,0])
		fc,m,p   = CS.detrend2d(years,fc)
		fc_sigma = fc.std(axis=0)
		mm.append(m)
		pp.append(p)

		Mm,Pm = np.zeros((CS.proj.nx,CS.proj.ny)),np.zeros((CS.proj.nx,CS.proj.ny))
		for i in range(CS.proj.nx):
			for j in range(CS.proj.ny):
				slope,intercept,r_value,p_value,std_err = stats.linregress(htmoor,fc[:,i,j])
				Mm[i,j] = slope
				Pm[i,j] = p_value/2.
		MM.append(Mm)
		PP.append(Pm)
	
		pl.figure(1)
		cseq  = np.arange(-60,60+10,10)
		sx,sy = np.ma.masked_where(p/2.>0.05,CS.proj.x),np.ma.masked_where(p/2.>0.05,CS.proj.y)
		cf    = pl.contourf(CS.proj.x,CS.proj.y,360*m*10,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Number denisty [%s$^{-1}$ decade$^{-1}$]' % (Season))
		pl.plot(sx,sy,'k.',alpha=0.5)
		CS.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(CS.proj.m,color='0.5',linewidth=0.7)
		pl.title(method)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/cyclones/N.trend.%s-%s.%s.%s.pdf' % (YearRange[0],YearRange[1],Season,method))
		#pl.close()

		pl.figure(2)
		cseq  = np.arange(-30,30+5,5)
		sx,sy = np.ma.masked_where(Pm>0.05,CS.proj.x),np.ma.masked_where(Pm>0.05,CS.proj.y)
		cf    = pl.contourf(CS.proj.x,CS.proj.y,htmoor_sigma*4*120*1e12*Mm,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Reg. coeff. [%s$^{-1}$ %sTW$^{-1}$]' % (Season,round(htmoor_sigma,2)))
		pl.plot(sx,sy,'k.',alpha=0.5)
                urot,vrot = CS.proj.m.rotate_vector(mus,mvs,CS.proj.lon,CS.proj.lat,returnxy=False)
                rotmask   = (pus>0.05)&(pvs>0.05)
                urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
                Q  = pl.quiver(CS.proj.x[::1,::1],CS.proj.y[::1,::1],htmoor_sigma*urot[::1,::1],htmoor_sigma*vrot[::1,::1],pivot='tail',alpha=0.95)
		#qk = pl.quiverkey(Q, 0.2, 1.02, 0.00005, '%s%s' % (0.05,'m s$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
       	 	CS.proj.m.drawparallels([70,80],latmax=90)
        	drawCoastlinesNoRivers(CS.proj.m,color='0.5',linewidth=0.7)
		pl.title(method)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/cyclones/N_ht.moor.regress.%s-%s.%s.%s.pdf' % (YearRange[0],YearRange[1],Season,method))
		#pl.close()

		pl.show()
	"""

	"""
	mm,pp     = np.array(mm),np.array(pp)
        MM,PP     = np.array(MM),np.array(PP)
	sxmm,symm = stipling(mm,xx=CS.proj.x,yy=CS.proj.y,thresh=0.8)
	sxMM,syMM = np.ma.masked_where(PP.mean(axis=0)>0.2,CS.proj.x),np.ma.masked_where(PP.mean(axis=0)>0.2,CS.proj.y)#stipling(MM,xx=CS.proj.x,yy=CS.proj.y,thresh=0.8)
	mm,pp     = mm.mean(axis=0),pp.mean(axis=0)
	MM,PP     = MM.mean(axis=0),PP.mean(axis=0)

	pl.figure(1)
	cseq  = np.arange(-60,60+10,10)
	cf    = pl.contourf(CS.proj.x,CS.proj.y,360*mm*10,cseq,cmap=pl.cm.coolwarm,extend='both')
	cbar  = pl.colorbar(cf)
	cbar.set_label('Number denisty [%s$^{-1}$ decade$^{-1}$]' % (Season))
	pl.plot(sxmm,symm,'k.',alpha=0.5)
	CS.proj.m.drawparallels([70,80],latmax=90)
	drawCoastlinesNoRivers(CS.proj.m,color='0.5',linewidth=0.7)
	pl.title('all')
	pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/cyclones/N.trend.%s-%s.%s.all.pdf' % (YearRange[0],YearRange[1],Season))

	pl.figure(2)
	cf    = pl.contourf(CS.proj.x,CS.proj.y,htmoor_sigma*360*1e12*MM,13,cmap=pl.cm.coolwarm,extend='both')
	cbar  = pl.colorbar(cf)
	cbar.set_label('Reg. coeff. [%s$^{-1}$ %sTW$^{-1}$]' % (Season,round(htmoor_sigma,2)))
	pl.plot(sxMM,syMM,'k.',alpha=0.5)
	CS.proj.m.drawparallels([70,80],latmax=90)
	drawCoastlinesNoRivers(CS.proj.m,color='0.5',linewidth=0.7)
	pl.title('all')
	pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/cyclones/N_ht.moor.regress.%s-%s.%s.all.pdf' % (YearRange[0],YearRange[1],Season))

	pl.show()
	"""


