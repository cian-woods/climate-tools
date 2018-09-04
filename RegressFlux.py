from TrackDataServer import DataServer as TrackDataServer
from ReanalysisDataServer import DataServer as reDataServer
from LambertProjector import *
from UnPickle import *
from scipy import interpolate

import numpy as np
import matplotlib.pyplot as pl

class RegressFlux:

	def __init__(	self,
			LatRange = (80,90),
			LonRange = (-180,180),
			Type     = 'fwrd'	):

		# TrackDataServer
		self.td = TrackDataServer(Type=Type,steps=25)
                # Attributes
                self.LatRange = LatRange
                self.LonRange = LonRange
		self.blat     = LatRange[0]
                # LambertProjector
                self.proj      = LambertProjector(boundinglat=70,resolution=200.)
                self.x,self.y  = self.proj.x[0,:],self.proj.y[:,0]	
                # Mask
                box = [LatRange,LonRange]
                if box[1][0] >= box[1][1]: lonA,latA = (self.proj.lon>box[1][1]) & (self.proj.lon<box[1][0]),(self.proj.lat>box[0][1]) | (self.proj.lat<box[0][0])
                if box[1][0] <  box[1][1]: lonA,latA = (self.proj.lon>box[1][1]) | (self.proj.lon<box[1][0]),(self.proj.lat>box[0][1]) | (self.proj.lat<box[0][0])
                self.mask    = (lonA | latA)
		self.Nmasked =  len(np.where(self.mask==False)[0])
		#x,y = np.ma.masked_array(self.proj.x,mask=self.mask),np.ma.masked_array(self.proj.y,mask=self.mask)
		#pl.plot(x,y,'k.')
		#self.proj.m.drawcoastlines()
		#self.proj.m.drawparallels([70,80],latmax=90)
		#pl.title('%s points; %s km resolution' % (self.Nmasked,self.proj.res))
		#pl.show()

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

	def filterIntrusions(self,G,Q,D,LON,LAT,P):
		gg,qq,dd  = [],[],[]
		lon,lat,p = [],[],[]
		for ii in range(len(G)):
		        Ntot,Nlat = 0,0
		        for t in range(len(G[ii])):
		                Ntot = Ntot + len(G[ii][t])
		                Nlat = Nlat + len(np.where((LAT[ii][t]>=80).mean(axis=1))[0])
		        if 1.*Nlat/Ntot > 0.4:
		                gg.append(G[ii])
		                qq.append(Q[ii])
		                dd.append(D[ii])
		                lon.append(LON[ii])
		                lat.append(LAT[ii])
		                p.append(P[ii])
		return gg,qq,dd,lon,lat,p

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
				snap      = self.td.snapshot(*D[ii][jj])
				lon,lat,p = snap[:,:,1],snap[:,:,0],snap[:,:,2]
				for kk in range(len(G[ii][jj])):
#				x = [int(G[ii][jj][kk]) for kk in range(len(G[ii][jj]))]
					x = int(G[ii][jj][kk])
					LON[ii][jj].append(lon[x,:])
					LAT[ii][jj].append(lat[x,:])
					P[ii][jj].append(p[x,:])
				LON[ii][jj] = np.array(LON[ii][jj])
				LAT[ii][jj] = np.array(LAT[ii][jj])
				P[ii][jj]   = np.array(P[ii][jj])
		return LON,LAT,P

	def tracksNotMasked(self,Season='DJF',flux=200,dur=12,rad=300,steps=np.arange(0,0+1,1)):
		# Returns (time)x(coordinate) array of all trajectory points in unmasked region
		# If closest gridpoint in unmasked region then use trajectory point

		# Open trajectory and flux files
		years = range(1980,2015+1,1)
       		Tracks     = [self.td.getDataSnaps(year,Season=Season) for year in years]
       		Tracks     = np.array([Tracks[i][j] for i in range(len(Tracks)) for j in range(len(Tracks[i]))])
		vq,Dates   = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1980-2015.30-1000hPa.70N.%s.p' % (Season))
		# Extract sector
		inds       = np.arange(0,359+1,1)
		Tracks     = Tracks[:,inds,:,:]
		vq         = vq[:,inds]
		# Get trajectory components
		LON,LAT,P  = Tracks[:,:,:dur,1],Tracks[:,:,:dur,0],Tracks[:,:,:dur,2]
		if flux != 'all':
			lon,lat,p  = [],[],[]
			for i in range(len(LON)):
				lon.append([])
				lat.append([])
				p.append([])
				for j in range(len(LON[i])):
					if vq[i,j] > flux:
						lon[i].append(LON[i][j])
						lat[i].append(LAT[i][j])
						p[i].append(P[i][j])
			LON,LAT,P = lon,lat,p
		print len(vq),len(LON)	
		# Find trajectory points inside unmasked region
		time,points = [],[]
		for i in range(len(LON)):
			print '%s of %s ...' % (i+1,len(LON))
			xs,ys = self.proj.m(np.array(LON[i]),np.array(LAT[i]))	
			for l in range(len(xs)):
		               	for k in range(len(xs[l])):
		               	        x,y    = xs[l,k],ys[l,k]
					ixiy   = zip(*np.where(np.sqrt((self.proj.x-x)**2+(self.proj.y-y)**2)/1000.<rad))
		               	        #ix,iy = np.argmin(np.abs(self.x-x)),np.argmin(np.abs(self.y-y))
					for iy,ix in ixiy:
		               	        	if not self.mask[iy,ix]:
							for m in steps:	
								#m    = 0
								hour = self.td.getHours(*Dates[i]) + (k+m)*6
								time.append(hour)
								points.append([iy,ix])	
		return time,points

        def intrusionsNotMasked(self,Season='DJF',flux=200,dur=12,rad=300,steps=np.arange(0,0+1,1)):
                # Returns (time)x(coordinate) array of all trajectory points in unmasked region
                # If closest gridpoint in unmasked region then use trajectory point 
                G,Q,Dates           = unpick('../intrusions/ERAInt_intrusions.%s.6x6hr.9deg.%s.6dt.20.5.p' % (Season,flux))
                LON,LAT,P           = self.getIntrusionTrajectories(G,Dates)
                G,Q,Dates,LON,LAT,P = self.filterIntrusions(G,Q,Dates,LON,LAT,P)
                # Find trajectory points inside unmasked region
                time,points = [],[]
                for i in range(len(LON)):
                        print '%s of %s ...' % (i+1,len(LON))
                        for j in range(len(LON[i])):
                        	xs,ys = self.proj.m(LON[i][j],LAT[i][j])
				xs,ys = xs[:,:dur],ys[:,:dur]
                        	for l in range(len(xs)):
                        	        for k in range(len(xs[l])):
                        	                x,y    = xs[l,k],ys[l,k]
						ixiy   = zip(*np.where(np.sqrt((self.proj.x-x)**2+(self.proj.y-y)**2)/1000.<rad))
                        	                #ix,iy = np.argmin(np.abs(self.x-x)),np.argmin(np.abs(self.y-y))
						for iy,ix in ixiy:
                        	                	if not self.mask[iy,ix]:
                        	                        	for m in steps:
                        	                                	hour = self.td.getHours(*Dates[i][j]) + (k+m)*6
                        	                                	time.append(hour)
                        	                                	points.append([iy,ix]) 
                return time,points

	def sortTimePoint(self,time,points):
	        P = {}
	        for i in range(len(time)):
	                hour  = time[i]
	                point = points[i]
	                try:
	                        P[hour].append(point)
	                except:
	                        P[hour] = []
	                        P[hour].append(point)
	        return P.keys(),[P[i] for i in P]

	def getTrackPoints(self,Season='DJF',flux='all',dur=12,rad=300,steps=np.arange(0,0+1,1)):
		fname = 'points/track.points.80N.%s.%sTg.%sx6hr.%skm.%s-%ssteps.p' % (Season,flux,dur,rad,steps[0],steps[-1])
		if not os.path.isfile(fname):
			t,p = self.tracksNotMasked(Season,flux,dur,rad,steps)
			t,p = self.sortTimePoint(t,p)
                	# Unique time and points (no duplicates)
                	pp = []
                	for i in range(len(p)):
                	        pp.append([])
                	        for j in range(len(p[i])):
                	                if pp[i].count(p[i][j])==0:
                	                        pp[i].append(p[i][j])
                	p = pp
                	toPick([t,p],fname)
		else:
			t,p = unpick(fname)
		return t,p

	def getIntrusionPoints(self,Season='DJF',flux=200,dur=12,rad=300,steps=np.arange(0,0+1,1)):
		fname = 'points/intrusion.points.80N.%s.%sTg.%sx6hr.%skm.%s-%ssteps.p' % (Season,flux,dur,rad,steps[0],steps[-1])
                if not os.path.isfile(fname):
                        t,p = self.intrusionsNotMasked(Season,flux,dur,rad,steps)
                        t,p = self.sortTimePoint(t,p)
			# Unique time and points (no duplicates)
                        pp = []
                        for i in range(len(p)):
                                pp.append([])
                                for j in range(len(p[i])):
                                        if pp[i].count(p[i][j])==0:
                                                pp[i].append(p[i][j])
                        p = pp
                        toPick([t,p],fname)
                else:
                        t,p = unpick(fname)
                return t,p

        def getTrackFluxPDF(self,Field='fls',LevType='surface_forecast',Season='DJF',flux=200,dur=12,rad=300,steps=np.arange(0,0+1,1)):
                fname = 'points/pdf.tracks.%s.80N.%s.%sTg.%sx6hr.%skm.%s-%ssteps.p' % (Field,Season,flux,dur,rad,steps[0],steps[-1])
                if not os.path.isfile(fname):
                        # FluxFile
                        vq,Dates = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1980-2015.30-1000hPa.70N.DJF.p')
                        # DataServer
                        ds    = reDataServer(Field=Field,LevType=LevType,LatRange=(70,90))
                        Times = [ds.getHours(*date) for date in Dates]
                        # Time and points
                        time,points = self.getTrackPoints(Season,flux,dur=dur,rad=rad,steps=steps)
                        # Data for points
                        S = []
                        for t in range(len(time)):
                                print '%s of %s ...' % (t+1,len(time))
                                if time[t] in Times:
                                        date = ds.getDate(time[t])
                                        snap = self.proj(ds.snapshot(*date),ds.lon,ds.lat)
                                        for l in range(len(points[t])):
                                                iy,ix = points[t][l]
                                                S.append(snap[iy,ix])
                        h,e = np.histogram(S,50,normed=False)
                        toPick([h,e],fname)
                else:
                        h,e = unpick(fname)
                return h,e

	def getIntrusionFluxPDF(self,Field='fls',LevType='surface_analysis',Season='DJF',flux=200,dur=12,rad=300,steps=np.arange(0,0+1,1)):
		fname = 'points/pdf.intrusions.%s.80N.%s.%s.p' % (Field,Season,flux)
		if not os.path.isfile(fname):
			# FluxFile
			vq,Dates = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1980-2015.30-1000hPa.70N.DJF.p')
			# DataServer
			ds    = reDataServer(Field=Field,LevType=LevType,LatRange=(70,90))
			Times = [ds.getHours(*date) for date in Dates]
			# Time and points
			time,points = self.getIntrusionPoints(Season,flux,dur=dur,rad=rad,steps=steps)
			# Data for points
                	S = []
                	for t in range(len(time)):
                	        print '%s of %s ...' % (t+1,len(time))
                	        if time[t] in Times:
                	                date = ds.getDate(time[t])
                	                snap = self.proj(ds.snapshot(*date),ds.lon,ds.lat)
                	                for l in range(len(points[t])):
                	                        iy,ix = points[t][l]
                	                        S.append(snap[iy,ix])
			h,e = np.histogram(S,50,normed=False)
			toPick([h,e],fname)
		else:
			h,e = unpick(fname)
		return h,e

	def allPointPDF(self,Field,LevType,bins,YearRange=(1981,2016)):
		fname = '%s.%sN.%s-%s.DJF.p' % (Field,self.blat,YearRange[0],YearRange[1])
		if not os.path.isfile(fname):
			# DataServer
			ds    = reDataServer(Field=Field,LevType=LevType,LatRange=(70,90))
			data  = []
			years = range(YearRange[0],YearRange[1]+1,1)
			for year in years:
			        print year
        			data0 = self.proj(ds.getDataSnaps(Year=year,Season='DJF'),ds.lon,ds.lat)
        			data0 = np.ma.masked_array(data0,mask=np.tile(self.mask,(len(data0),1,1))).reshape(-1)
        			data0 = list(data0.data[np.where(data0.mask==False)])
        			data  = data + data0
			h0,e0   = np.histogram(data,bins,normed=False)
			toPick([h0,e0],fname)
		else:
			h0,e0 = unpick(fname)
		return h0,e0

	def regression1(self,Field='pw',LevType='surface_analysis'):
		# Trajectory parameters
		dur,rad,steps = 15,300,np.arange(0,0+1,1)
		# Time and points for intrusion trajectories (Lambert 70N, 200 res)
		h1,e1 = self.getTrackFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=50,dur=dur,rad=rad,steps=steps)
		h2,e2 = self.getTrackFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=100,dur=dur,rad=rad,steps=steps)
		h3,e3 = self.getTrackFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=200,dur=dur,rad=rad,steps=steps)
		h4,e4 = self.getTrackFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=300,dur=dur,rad=rad,steps=steps)
		h5,e5 = self.getTrackFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=400,dur=dur,rad=rad,steps=steps)
		h6,e6 = self.getIntrusionFluxPDF(Field=Field,LevType=LevType,Season='DJF',flux=200,dur=dur,rad=rad,steps=steps)
		# DataServer
                ds = reDataServer(Field=Field,LevType=LevType,LatRange=(70,90))	
		# Data for all points
		#bins    = np.linspace(-100,50,50)
		h0,e0   = self.allPointPDF(Field,LevType,50,YearRange=(1980,2015))
		#h00,e00 = self.allPointPDF(Field,LevType,bins,YearRange=(1999,2016))
		#pl.plot(e00[0:-1],h00,'r',linewidth=1.5,label='all points 1999-2016',alpha=0.5)
		pl.plot(e0[0:-1],h0,'k',linewidth=1.5,label='all points 1980-2015')
		pl.plot(e1[0:-1],h1,'b',label='$f>50$')
		pl.plot(e2[0:-1],h2,'g',label='$f>100$')
		pl.plot(e3[0:-1],h3,'y',label='$f>200$')
		pl.plot(e4[0:-1],h4,'orange',label='$f>300$')
		pl.plot(e5[0:-1],h5,'r',label='$f>400$')	
		pl.plot(e6[0:-1],h6,'grey',linestyle='dashed',label='Intrusions')
		lg = pl.legend(loc=0,frameon=False,prop={'size':12},ncol=2,columnspacing=1,handletextpad=0.2,title=None)
		pl.ylabel('Frequency')
		pl.xlabel('Net surface longwave radiation')
		pl.xlim((-150,100))
                #pl.xlabel('Precipitable water [kg m$^{-2}$]')
                #pl.xlim((0,12))
		pl.savefig('figs/%s_%sN.pdf' % (Field,self.blat),format='pdf')
		pl.show()	

        def interpTrack(self,xs,n):
                xold = np.linspace(0,1,len(xs))
                xnew = np.linspace(0,1,len(xs)*(n+1) - n)
                f  = interpolate.interp1d(xold,xs)
                xs = f(xnew)
                return xs

        def trackDensity(self,Xs,Ys,x_ax,y_ax):
		N = np.zeros((len(x_ax),len(y_ax)))
		for i in range(len(Xs)):
			xs,ys = Xs[i,:],Ys[i,:]
               		# Interpolate track to higher temporal res
                	xs,ys = self.interpTrack(xs,10),self.interpTrack(ys,10)  
                	# Data holders 
                	ind = []
                	for i in range(len(xs)):
				xi,yi = np.argmin(abs(x_ax-xs[i])),np.argmin(abs(y_ax-ys[i]))
				#if ind.count([yi,xi]) == 0:
				N[yi,xi]   = N[yi,xi] + 1    
				ind.append([yi,xi])
                return N

	def regression2(self,Season='DJF'):
		# back1 tracks were computed using blat,res = 80,200
		# LambertProjectors and DataServers
		proj      = LambertProjector(boundinglat=80,resolution=200.)
		npro      = LambertProjector(boundinglat=20,resolution=400.)
		x_ax,y_ax = npro.x[0,:],npro.y[:,0]
		ds        = reDataServer(Field='pw',LevType='surface_forecast',LatRange=(50,90))
		vq,Dates  = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1980-2015.30-1000hPa.70N.%s.p' % (Season))
		vqlons    = np.arange(0,360,1)
		N         = np.zeros((len(x_ax),len(y_ax)))
		Years     = range(2000,2000+1,1)
		for year in Years:
			print year
			# Dates and data
			dates     = ds.getDateList(Year=year,Month=1,Day=1)
			data      = proj(ds.getDataSnaps(Year=year,Month=1,Day=1),ds.lon,ds.lat)
			tracks    = self.td.getDataSnaps(Year=year,Month=1,Day=1)
			lat,lon,p = tracks[:,:,:,0],tracks[:,:,:,1],tracks[:,:,:,2]
			# Extract vq time slice
			#x0,x1 = Dates.index(dates[0]),Dates.index(dates[-1])
			#vqslc = vq[x0:x1+1,:]
			# Extract tracks above -20 W m**-2
			lat0s,lon0s,p0s = lat[:,:,0],lon[:,:,0],p[:,:,0]
			xs,ys           = proj.m(lon0s,lat0s)
			LON,LAT,P       = [],[],[]
			flux,sflux,llon	= [],[],[]
			for i in range(len(xs)):
				x0 = Dates.index(dates[i])
				for j in range(len(xs[i])):
					x,y   = xs[i,j],ys[i,j]
					ix,iy = np.argmin((proj.x[0,:]-x)**2),np.argmin((proj.y[:,0]-y)**2)
					f     = data[i,iy,ix]
					if (f > -20) and (proj.lat[iy,ix]>=80):
						x,flon  = self.interpolateTrackToLat(lon[i,j,:],lat[i,j,:],70)
						if flon != -999:
							ilon = np.argmin((vqlons-flon)**2)
							inds = np.arange(-4,4+1,1)+ilon
							for ij in range(len(inds)):
								if inds[ij] < 0:            inds[ij] = inds[ij]+len(vqlons)
								if inds[ij] >= len(vqlons): inds[ij] = inds[ij]-len(vqlons)+1
							print Dates[x0-x],inds,vq[x0-x,inds]
							flux.append(vq[x0-x,inds].mean())
							sflux.append(f)
							llon.append(flon)
		h,e = np.histogram(llon,40)
		pl.plot(e[0:-1],h,'k')
		pl.show()
		#pl.plot(sflux,flux,'k.',alpha=0.3)
		#pl.title(np.mean(flux))
		#pl.show()
		"""
						LON.append(lon[i,j,:])
						LAT.append(lat[i,j,:])
						P.append(p[i,j,:])
			LON,LAT,P = np.array(LON),np.array(LAT),np.array(P)

			# Get density of trajectories
			Xs,Ys = npro.m(LON,LAT)
			N    += self.trackDensity(Xs,Ys,x_ax,y_ax)
		# Plot
		pl.contourf(x_ax,y_ax,N,20,cmap=pl.cm.OrRd)
		pl.colorbar()
		npro.m.drawcoastlines()
		npro.m.drawparallels([70,80],latmax=90)
		pl.show()
		"""

	def interpolateTrackToLat(self,lons,lats,blat):
		# Interpolate track to the given latitude.
		# If the track does not reach given latitide reutrn None
		# return lat and long
		try:
			x   = np.where(lats<70)[0][0]
			lon = lons[x]
			if lon < 0: lon = lon + 360
		except:
			x   = -999
			lon = -999
		return x,lon

	def associateEvents(self,D,ds,dates,mind=-8,maxd=-2):
		# Associate intrusions to warm arctic events
		days0b = np.array([ds.getDays(*D[i][0]) for i in range(len(D))])       # Beginning days of injections
		days0e = np.array([ds.getDays(*D[i][-1]) for i in range(len(D))])      # Ending days of injections
		days   = np.array([ds.getDays(*date) for date in dates])                         # Warm event days
		xs,ls  = [],[]
		for i in range(len(days)):
		        day0  = days[i]                                                 # Day of warm event
		        dtEnd = days0e - day0                                           # Day difference between intrusion start dates and warm event
		        dtBeg = days0b - day0                                           # Day difference between intrusion end dates and warm event
		        xxs   = list(np.where(((dtBeg<mind)&(dtEnd>maxd))|((dtBeg>=mind)&(dtBeg<=maxd))|((dtEnd>=mind)&(dtEnd<=maxd))==True)[0])
		        xs   = xs + xxs
		        ls.append(xxs)
		xs = list(set(xs))                                                      # Exclude intrusions associated multiple times
		return xs,ls

	def regressFields(self,Season='DJF'):
			# Setup attributes
			proj    = LambertProjector(boundinglat=85,resolution=80.)
			#Years1 = range(1980,2007+1,1) + range(2009,2012+1,1)
			#Years1,Years2 = range(1980,2015+1,1),range(1998,2015+1,1)
			#Years1  = range(2006,2006+1,1)
			Years1 = range(1980,2015+1,1)
			# DataServers etc.
			f1,f2,f3 = 'fls','pw','tcc'
			d1       = reDataServer(Field=f1,Source='ERAInt',LatRange=(80,90),LonRange=(0,360),LevType='surface_forecast')
			d2       = reDataServer(Field=f2,Source='ERAInt',LatRange=(80,90),LonRange=(0,360),LevType='surface_analysis')
			#d3       = reDataServer(Field=f3,Source='ERAInt',LatRange=(80,90),LonRange=(0,360),LevType='surface_analysis')
			#clim     = np.array([d3.getDataSnaps(year,Season='DJF',dailymean=False,step=6).mean(axis=0) for year in range(1980,2015+1,1)]).mean(axis=0)
			# Intrusions associated
			#G,Q,D   = unpick('/qnap/cian/sheba/intrusionfiles/ERAInt_intrusions.DJF.6x6hr.9deg.200.6dt.20.5.filtered.23steps.80N.0.4.-100.0.full.p')
			G,Q,D   = unpick('../intrusions/ERAInt_intrusions.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p')
			#for year in Years:
			#Years1 = range(year,year+1,1)
			dates   = [d1.getDateList(Year=year,Season='DJF') for year in Years1]
			dates   = [dates[i][j] for i in range(len(dates)) for j  in range(len(dates[i]))]
			xs,ls1  = self.associateEvents(D,d1,dates,-8,-2)
			ls1     = (6/24.)*np.array([ [ [np.sum([len(D[xxsi]) for xxsi in xxs]) for jj in range(proj.ny)] for kk in range(proj.nx)] for xxs in ls1 ]).reshape(-1)
			ls1mean = ls1.mean()
                        #dates   = [d1.getDateList(Year=year,Season='DJF') for year in Years2]
                        #dates   = [dates[i][j] for i in range(len(dates)) for j  in range(len(dates[i]))]
                        #xs,ls2  = self.associateEvents(D,d1,dates,-8,-2)
                        #ls2     = np.array([ [ [np.sum([len(D[xxsi]) for xxsi in xxs]) for jj in range(proj.ny)] for kk in range(proj.nx)] for xxs in ls2 ]).reshape(-1)
			#ls2mean = ls2.mean()
			print ls1mean#,ls2mean
			# Get data
			s1,s2,s3 = [],[],[]
			t1,t2,t3 = [],[],[]
			for year in Years1:
				print year
			        s1.append(proj(d1.getDataSnaps(Year=year,Season=Season),d1.lon,d1.lat))
				s2.append(proj(d2.getDataSnaps(Year=year,Season=Season),d2.lon,d2.lat))
                                #s3.append(proj(d3.getDataSnaps(Year=year,Season=Season)-clim[np.newaxis,:,:],d3.lon,d3.lat)) 
                        #for year in Years2:
                        #        print year
                        #        t1.append(proj(d1.getDataSnaps(Year=year,Season=Season),d1.lon,d1.lat))
                        #        t2.append(proj(d2.getDataSnaps(Year=year,Season=Season),d2.lon,d2.lat))
			s1  = np.array([s1[i][j] for i in range(len(s1)) for j in range(len(s1[i]))]).reshape(-1)
			s2  = np.array([s2[i][j] for i in range(len(s2)) for j in range(len(s2[i]))]).reshape(-1)
			#s3  = np.array([s3[i][j] for i in range(len(s3)) for j in range(len(s3[i]))]).reshape(-1)
                        #t1  = np.array([t1[i][j] for i in range(len(t1)) for j in range(len(t1[i]))]).reshape(-1)
                        #t2  = np.array([t2[i][j] for i in range(len(t2)) for j in range(len(t2[i]))]).reshape(-1)
			xedges  = np.linspace(0,10,41)
			#xedges  = np.linspace(200,280,41)
			yedges  = np.linspace(-100,75,41)
			#yedges = np.linspace(200,280,41)
			#yedges  = np.linspace(0,350,41)
			#x,y,H1  = self.bivarPDF(s1,s2,xedges,yedges,norm=True)
			x,y,H1,N1,Nf1 = self.bivar2d(s2,s1,ls1,xedges,yedges)
			#x,y,H2,N2 = self.bivar2d(t2,t1,ls2,xedges,yedges)
			fend = 'all'
                        pl.figure(1) 
                        pl.pcolor(x,y,H1-ls1mean,cmap=pl.cm.RdBu_r,vmin=-6,vmax=6)
                        cbar = pl.colorbar(extend='both')
                        cbar.set_label('Mean injection association anomaly [days]')
                        pl.xlabel('%s [%s]' % (d2.long_name,d2.units))
                        pl.ylabel('%s [%s]' % (d1.long_name,d1.units))
                        pl.title('All points north of 80N [%s]' % (Season))
                        pl.savefig('figs/regress/%s_%s.%s.anom.H.pdf' % (f1,f2,fend),format='pdf')
                        #pl.close()
			pl.figure(2)
			pl.pcolor(x,y,H1,cmap=pl.cm.OrRd,vmin=0,vmax=16)
			cbar = pl.colorbar(extend='both')
			cbar.set_label('Mean injection association anomaly [days]')
			pl.xlabel('%s [%s]' % (d2.long_name,d2.units))
			pl.ylabel('%s [%s]' % (d1.long_name,d1.units))
			pl.title('All points north of 80N [%s]' % (Season))
			pl.savefig('figs/regress/%s_%s.%s.full.H.pdf' % (f1,f2,fend),format='pdf')
			#pl.close()
                	pl.figure(3)
			#pl.pcolor(x,y,N1,cmap=pl.cm.RdBu_r,vmin=None, vmax=None)
			pl.contourf(x,y,Nf1,40,cmap=pl.cm.RdBu_r)
			xx,yy = np.meshgrid(x,y)
			sx,sy = np.ma.masked_where(Nf1>10,xx),np.ma.masked_where(Nf1>10,yy)
			#pl.plot(sx,sy,'k.',alpha=0.5)
                	cbar = pl.colorbar()
                	cbar.set_label('Number frequency')
                	pl.xlabel('%s [%s]' % (d2.long_name,d2.units))
                	pl.ylabel('%s [%s]' % (d1.long_name,d1.units))
                	pl.title('All points north of 80N [%s]' % (Season))
                	pl.savefig('figs/regress/%s_%s.%s.N.pdf' % (f1,f2,fend),format='pdf')
			pl.show()
			#pl.close()

	def bivarPDF(self,x,y,xedges,yedges,norm=True):
                H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=None,normed=norm,weights=None)
                xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
		return xedges,yedges,H

	def bivar2d(self,x,y,z,xedges,yedges):
		H,N = np.zeros((len(xedges),len(yedges))),np.zeros((len(xedges),len(yedges)))
		for i in range(len(x)):
			ix,iy = np.argmin((xedges-x[i])**2),np.argmin((yedges-y[i])**2)
			H[iy,ix] = H[iy,ix] + z[i]
			N[iy,ix] = N[iy,ix] + 1
		H,N = np.ma.masked_where(N==0,H),np.ma.masked_where(N==0,N)
		return xedges,yedges,H/N,N/N.sum(),N
			
	def plotTracks(self,lats1,lons1,lats2,lons2):
		# Centroids
		Xs1,Ys1 = [],[]
		Xs2,Ys2 = [],[]	
		for t in range(len(lats1)):
		        xs,ys = self.proj.m(lons1[t],lats1[t])
		        xs,ys = xs.mean(axis=0)[:21],ys.mean(axis=0)[:21]
		        Xs1.append(xs)
		        Ys1.append(ys)
		Xs1,Ys1 = np.array(Xs1),np.array(Ys1)
                for t in range(len(lats2)):
                        xs,ys = self.proj.m(lons2[t],lats2[t])
                        xs,ys = xs.mean(axis=0)[:21],ys.mean(axis=0)[:21]
                        Xs2.append(xs)
                        Ys2.append(ys)
                Xs2,Ys2 = np.array(Xs2),np.array(Ys2)
		# Plot
		for i in range(len(Xs1)):
		        pl.plot(Xs1[i],Ys1[i],'r')
                for i in range(len(Xs1)):
                        pl.plot(Xs2[i],Ys2[i],'b')
		self.proj.m.drawcoastlines()
		self.proj.m.drawparallels([60,70,80],latmax=90)
		pl.show()

if __name__ == "__main__":

	from toPick import *

	rf = RegressFlux(Type='back')

	#rf.regression2()
	rf.regression1('fls','surface_forecast')
	#rf.regressFields()	
