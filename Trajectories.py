from cmipDataServer import DataServer as cmipDataServer
from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from scipy import interpolate,stats,mgrid,signal
from scipy import ndimage as nd
from UnPickle import *
from toPick import *
from datesToN import *
from stipling import stipling
from trajNetCDF import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers

import matplotlib.pyplot as pl
import sys,glob,os
import numpy as np
import math

class Trajectories:

	def __init__(	self,
			Source    = 'CMCC-CESM',
			dT        =  1,
			ff        =  1,
			blat	  =  70,
			res       =  200	):

		# Class attributes
		self.Source    = Source
		self.dT        = dT
		self.ff        = ff
		# Dependants
		if self.ff ==    1: self.fstr      = 'fwrd'
		if self.ff ==   -1: self.fstr      = 'back'
		if self.dT ==    1: self.dailymean = True
		if self.dT == 0.25: self.dailymean = False
		# U,V and W DataServers
		if Source == 'ERAInt':
                        self.uds  = reDataServer(Field='U',LevType='plev',Source=Source)
                        self.wds  = reDataServer(Field='W',LevType='plev',Source=Source)
                        self.vds  = reDataServer(Field='V',LevType='plev',Source=Source)
                        self.PVds = reDataServer(Field='pv',LevType='plev',Source=Source)
		elif Source == 'NCEP1':
                        self.uds  = reDataServer(Field='uplev',LevType='plev',Source=Source)
                        self.wds  = reDataServer(Field='omega',LevType='plev',Source=Source)
                        self.vds  = reDataServer(Field='vplev',LevType='plev',Source=Source)
                        #self.qdsh = reDataServer(Field='shum',LevType='plev',Source=Source)
                elif Source[0:3] == 'CAM':
                        self.uds  = reDataServer(Field='U',Source=Source)
                        self.wds  = reDataServer(Field='W',Source=Source)
                        self.vds  = reDataServer(Field='V',Source=Source)
                        #self.qdsh = reDataServer(Field='q',Source=Source)
		else:
			self.wds  = cmipDataServer(Field='wap',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day')
			self.uds  = cmipDataServer(Field='ua' ,LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day')
			self.vds  = cmipDataServer(Field='va' ,LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day')
			#self.qdsh = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='historical',DataFreq='day')
			#self.qdsr = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day')
			#self.qdss = [self.qdsh,self.qdsr]
		# LambertProjector
		self.res        = res
		self.proj       = LambertProjector(boundinglat=blat,resolution=self.res)
		self.res        = self.proj.res
		print self.res
		self.x,self.y   = self.proj.x[0,:],self.proj.y[:,0]
		self.nx,self.ny = len(self.x),len(self.y)
		# Mass of girdbox on each level
		self.dP = []
		dP = np.diff(self.uds.lev)/2.
		self.dP.append(dP[0])
		for i in range(len(dP)-1):
	        	self.dP.append(dP[i]+dP[i+1])
		self.dP.append(dP[-1])
		self.dP = np.array(self.dP)
		self.dM = 100.*self.dP/9.80665

	def data(self,datelist):
		# Datelist spans injection lifetime plus N days for trajectories
		# dt: data is daily, interpolate to dt hour time frequency	

		# Get wind fields
		udata  = self.uds.getDataSnaps(datelist=datelist,dailymean=self.dailymean)
		vdata  = self.vds.getDataSnaps(datelist=datelist,dailymean=self.dailymean)
		wdata  = self.wds.getDataSnaps(datelist=datelist,dailymean=self.dailymean)
		PVdata = self.PVds.getDataSnaps(datelist=datelist,dailymean=self.dailymean)/(1e-06)
		# Map to projection
		udata  = self.proj(udata,self.uds.lon,self.uds.lat)
		vdata  = self.proj(vdata,self.vds.lon,self.vds.lat)
		wdata  = self.proj(wdata,self.wds.lon,self.wds.lat)
		PVdata = self.proj(PVdata,self.PVds.lon,self.PVds.lat)
		# Rotate horizonal components
		H = []
		for j in range(len(datelist)):
			h = [self.proj.m.rotate_vector(udata[j][i],vdata[j][i],self.proj.lon,self.proj.lat,returnxy=False) for i in range(self.uds.nlev)]
			H.append(h)
		H = np.array(H)
		udata,vdata = H[:,:,0,:,:],H[:,:,1,:,:]
		if self.ff == -1: udata,vdata,wdata = -1*udata[:,:,:,:],-1*vdata[:,:,:,:],-1*wdata[:,:,:,:],
		return udata,vdata,wdata,PVdata

	def interpolateField(self,x0,y0,p0,field,levels):

		# Vertical interpolation
		if p0 in self.uds.lev:
			l = np.argmin(abs(levels-p0))
			field = field[l,:,:]
		else:
			l1 = np.argmin(np.abs(levels - p0))
			if levels[l1] < p0: l2 = l1 + 1
			if levels[l1] > p0: l2 = l1 - 1
			pb    = np.array([levels[l1],levels[l2]]).max()
			pt    = np.array([levels[l1],levels[l2]]).min()
			lb    = np.argmin(np.abs(levels-pb))
			lt    = np.argmin(np.abs(levels-pt))
			dP    = pb - pt
			dF    = field[lb,:,:] - field[lt,:,:]
			field = field[lt,:,:] + (p0-pt)*dF/dP		
		# Horizontal interpolation
		rbs = interpolate.RectBivariateSpline(self.x,self.y,field)
		d   = rbs.ev(y0,x0)
		return float(d)

	def track0(self,x0,y0,p0,u,v,w,pv,dt=24):
		# Input particle position and current wind fields

		# Find winds at particle position
		u0  = self.interpolateField(x0,y0,p0,u,self.uds.lev)
		v0  = self.interpolateField(x0,y0,p0,v,self.vds.lev)
		w0  = self.interpolateField(x0,y0,p0,w,self.wds.lev)
		# Intergrate in time
		x1 = x0 + u0*dt*60*60		# m
		y1 = y0 + v0*dt*60*60		# m
		p1 = p0 + w0*dt*60*60/100.	# hPa
		# Set verticle coordinate bounds
		if p1 > self.uds.lev[-1]:
			p1 = self.uds.lev[-1]
		elif p1 < self.uds.lev[0]:
			p1 = self.uds.lev[0]
		pv1 = self.interpolateField(x1,y1,p1,pv,self.PVds.lev)
		return x1,y1,p1,pv1

	def trackN(self,x0,y0,p0,u,v,w,pv,dt=1):

		N = int(self.dT*24./dt)
		xs,ys,ps = [x0],[y0],[p0]
		pv0      = self.interpolateField(x0,y0,p0,pv[0],self.PVds.lev)
		pvs      = [pv0]
		for i in range(len(u)-1):
			du  =   (u[i+1] - u[i])/N
			dv  =   (v[i+1] - v[i])/N
			dw  =   (w[i+1] - w[i])/N
			dpv = (pv[i+1] - pv[i])/N
			for j in range(0,N):
				un           =   u[i] + du*j
				vn           =   v[i] + dv*j
				wn           =   w[i] + dw*j
				pvn          =  pv[i] + dpv*(j+1)
				x0,y0,p0,pv0 = self.track0(x0,y0,p0,un,vn,wn,pvn,dt=dt)
				xs.append(x0)
				ys.append(y0)
				ps.append(p0)
				pvs.append(pv0)
		return xs,ys,ps,pvs

	def tracksN(self,x0s,y0s,p0s,u,v,w,pv,dt=1):
		# Takes in set of initial positions and chunk of data for making trajectories
		X,Y,P,PV = [],[],[],[]
		for i in range(len(x0s)):
			x0,y0,p0     = x0s[i],y0s[i],p0s[i]
			xs,ys,ps,pvs = self.trackN(x0,y0,p0,u,v,w,pv,dt=dt)
			X.append(xs)
			Y.append(ys)
			P.append(ps)
			PV.append(pvs)
		return X,Y,P,PV

	def makeDateList(self,date0,Ndays=5):
		# Extend date0 forward by number of days
		day0     = self.uds.getDays(*date0)
		datelist = [self.uds.getDate(day0+self.ff*i,nd=24) for i in np.arange(0,Ndays+self.dT,self.dT)]	
		return datelist

	def intrusionTracks(self,g,d,dt=1,Ndays=5):

		# Ndays for data
		N = len(d)*self.dT + Ndays
		# Get datelist for computation
		datelist = self.makeDateList(d[0],Ndays=N)	
		# Get data and particle positions at t=0
		# 70N and 900hPa
		u,v,w = self.data(datelist)

		TX,TY,TP = [],[],[]
		for i in range(len(d)):
			x0s,y0s = self.proj.m(g[i],[70 for kk in g[i]])
			p0s     = [900 for kk in g[i]]
			# Compute trajectories for first timestep of intrusion
			X,Y,P     = self.tracksN(x0s,y0s,p0s,u[i:Ndays/self.dT+1+i,:,:,:],v[i:Ndays/self.dT+1+i,:,:,:],w[i:Ndays/self.dT+1+i,:,:,:],dt=dt)
			lons,lats = self.proj.m(np.array(X),np.array(Y),inverse=True)
			TX.append(lons)
			TY.append(lats)
			TP.append(P)
		return TX,TY,TP

	def fill(self,data,invalid=None):
		if invalid is None: invalid = np.isnan(data)
		ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
		return data[tuple(ind)]

	def allModelField(self,Sector=(40,60)):
		Models = [g[9:] for g in glob.glob('../rcp/*')][:2]
		date0h,date0r = [],[]
		for Model in Models:
			dh,dr = self.fieldAlongTracks(self,Model=Model,Sector=(40,60))
			date0h.append(dh.mean(axis=1))
			date0r.append(dr.mean(axis=1))
		date0h = np.array(date0h).mean(axis=0)
		date0r = np.array(date0r).mean(axis=0)
		# Plot
		days  = range(-5,1)
		pl.plot(days,date0h,'b',linewidth=2)
		pl.plot(days,date0h-datash,'b--',linewidth=1)
		pl.plot(days,date0h+datash,'b--',linewidth=1)
		pl.plot(days,date0r,'r',linewidth=2)
		pl.plot(days,date0r-datasr,'r--',linewidth=1)
		pl.plot(days,date0r+datasr,'r--',linewidth=1)
		pl.show()

	def allField(self):
		Models = [g[9:] for g in glob.glob('../rcp85/*')]
		Dh,Dr  = [],[]
		for Model in Models:
			try:
				dh1,dr1 = unpick('alongtrack/%s..p' % (Model))
				dh2,dr2 = unpick('alongtrack/%s.back_.p' % (Model))
				dh,dr   = list(dh2[::-1])+list(dh1)[1:],list(dr2[::-1])+list(dr1)[1:]
				Dh.append(dh)
				Dr.append(dr)
			except:
				pass
		Dh,Dr = np.array(Dh).mean(axis=0),np.array(Dr).mean(axis=0)
		pl.plot(Dh,'b')
		pl.plot(Dr,'r')
#		pl.plot(Dr-Dh,'b')
#		pl.plot([0,6],[0,0],'k--')
		pl.show()

	def fieldAlongTracks(self,Model='CMCC-CESM',Sector=(40,60)):
		# Trajectories and injections
		Gh,Qh,Dh     = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
		LONh,LATh,Ph = unpick('/qnap/cian/cmip/scripts/intrusiontracks/%s_%s_tracks.6hour.p' % (Model,self.fstr))
		Gh,Qh,Dh,LONh,LATh,Ph = self.filterInjections(Gh,Qh,Dh,LONh,LATh,Ph,Sector=Sector)
                Gr,Qr,Dr     = unpick('../intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
                LONr,LATr,Pr = unpick('intrusiontracks/%s_%s_tracks.6hour.p' % (Model,self.fstr))
		Gr,Qr,Dr,LONr,LATr,Pr = self.filterInjections(Gr,Qr,Dr,LONr,LATr,Pr,Sector=Sector)
		# Get field
		data0h = self.getField(Gh,Qh,Dh,LONh,LATh,Ph,0)
		data0r = self.getField(Gr,Qr,Dr,LONr,LATr,Pr,1)
		data0h = data0h.mean(axis=1).squeeze()
		data0r = data0r.mean(axis=1).squeeze()
		return data0h,data0r

	def getField(self,G,Q,D,LON,LAT,P,ind=0):
#		if ind==0: ds = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='historical',DataFreq='day')
#               if ind==1: ds = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day')
		ds        = self.qdss[ind]
		data0     = [[] for i in range(6)]
		for i in range(len(LON)):	
			for t in range(len(LON[i])):
				date0        = D[i][t]	
				day          = 0
				lons,lats,ps = LON[i][t],LAT[i][t],P[i][t]
				xs,ys        = self.proj.m(lons,lats)
				xs,ys,ps     = xs[:,:21],ys[:,:21],np.array(ps)[:,:21]
				xs           = np.array(xs).mean(axis=0)[::4] # centroid
				ys           = np.array(ys).mean(axis=0)[::4] # centroid
				ps           = np.array(ps).mean(axis=0)[::4] # centroid
				for k in range(len(xs)):
					x0,y0,p0 = xs[k],ys[k],ps[k]	
					date     = ds.getDate(ds.getDays(*date0) + day)		
					data     = ds.snapshot(*date)
                                	if (type(data)==np.ma.core.MaskedArray) and (data.mask.shape == data.data.shape):
                                        	data = self.fill(data.data,invalid=data.mask)
					data = (data*self.dM[:,np.newaxis,np.newaxis]).sum(axis=0)
					data = self.proj(data,ds.lon,ds.lat)
					data_,x,y = self.interp2dN(data,self.x,self.y,x0,y0,kind='linear')
#					data_     = self.interpolateField(x0,y0,p0,data)
					data0[k].append(data_)
					day       = day + self.ff
		data0 = np.array(data0)
		return data0

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

	def interpTrack(self,xs,n):

		xold = np.linspace(0,1,len(xs))
		xnew = np.linspace(0,1,len(xs)*(n+1) - n)

		f  = interpolate.interp1d(xold,xs)
		xs = f(xnew)
		return xs

	def trackDensity(self,xs,ys):

		# Centroid trajectories
		xs  = np.array(xs).mean(axis=0) # centroid
		ys  = np.array(ys).mean(axis=0) # centroid
		# Interpolate tracks to higher temporal res
		#xs,ys = self.interpTrack(xs,10),self.interpTrack(ys,10)
		# Tangents
		yts,xts = self.getTangents(xs,ys,mag=1)
		# Data holders
		N   = np.zeros((self.nx,self.ny))
		T   = np.zeros((self.nx,self.ny,2))
		ind = []
		for i in range(len(xs)):
				xi,yi = np.argmin(abs(self.x-xs[i])),np.argmin(abs(self.y-ys[i]))
			#if ind.count([yi,xi]) == 0:
				N[yi,xi]   = N[yi,xi] + 1
				T[yi,xi,0] = T[yi,xi,0] + xts[i]
				T[yi,xi,1] = T[yi,xi,1] + yts[i]
			#	ind.append([yi,xi])
		Ntile = np.tile(N[:,:,np.newaxis],(1,1,2))
		T     = np.ma.masked_where(Ntile==0,T)/np.ma.masked_where(Ntile==0,Ntile)
		T     = T.data
		return N,T

	def density(self,LON,LAT,nend=None):
		# Intrusions trajectories
		nn,Ttot = [],[]
		for i in range(len(LON)):
#			print '%s of %s ...' % (i+1,len(LON))
			N,T = [],[]
			for t in range(len(LON[i])):
				lons,lats = LON[i][t],LAT[i][t]
				xs,ys     = self.proj.m(lons,lats)
				xs,ys     = xs[:,:nend],ys[:,:nend]
				Nit,Tit   = self.trackDensity(xs,ys)
				N.append(Nit)
				Ttot.append(Tit)
			N = np.array(N).sum(axis=0)
			#N = N/len(LON[i])
			#N[np.where(N>=1)] = 1
			nn.append(N)
		nn   = np.array(nn)
		Ntot = nn.sum(axis=0)
		Ttot = np.array(Ttot)
		Ttot = np.ma.masked_where(Ttot==0,Ttot).mean(axis=0)
		return nn,Ntot,Ttot

	def interp2dN(self,field,x,y,new_x,new_y,kind='linear'):
		# 2d interp function
		f = interpolate.interp2d(x,y,field,kind=kind,bounds_error=True)
		# Interpolate to new grid
		field = f(new_x,new_y)
		return field,new_x,new_y

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

	def filterTracks(self,LON,LAT,P,blat=80.):

		lon,lat,p = [],[],[]
		for i in range(len(LAT)):
			n0,n1 = 0,0
			for t in range(len(LAT[i])):
				for j in range(len(LAT[i][t])):
					n0 = n0 + 1
					if (LAT[i][t][j] >= blat).any():
						n1 = n1 + 1
			if 1.*n1/n0 >= 0.4:
				lon.append(LON[i])
				lat.append(LAT[i])
				p.append(P[i])
		return lon,lat,p

	def sectorProj(self,data,LatRange=(70,90),LonRange=(-30,60),plon=None,plat=None,x=None,y=None):
		# LonRange should be -180 --> 180 style
		if (x != None) and (y != None):
			xx = np.array([x for xii in range(len(x))])
			yy = np.array([np.zeros(len(x)) + x[xii] for xii in range(len(x))])
			plon,plat = self.proj.m(xx,yy,inverse=True)
		if (plon == None) and (plat == None):
			plon,plat = self.proj.lon,self.proj.lat
		data = np.ma.masked_where(plat<LatRange[0],data)
		data = np.ma.masked_where(plat>LatRange[1],data)
		if LonRange[0] < LonRange[1]:
			data = np.ma.masked_where(plon<LonRange[0],data)
			data = np.ma.masked_where(plon>LonRange[1],data)
		if LonRange[1] < LonRange[0]:
			data = np.ma.masked_where(plon<LonRange[0],data)
			data = np.ma.masked_where(plon<LonRange[1],data)
		return data

	def plotTrack(self,lats,lons,dt=0.25):
		# Centroids
		Xs,Ys = [],[]
		for t in range(len(lats)):
			xs,ys = self.proj.m(lons[t],lats[t])
			xs,ys = xs.mean(axis=0)[:20],ys.mean(axis=0)[:20]
			Xs.append(xs)
			Ys.append(ys)
		Xs,Ys = np.array(Xs),np.array(Ys)
                # Start times
                starts = np.arange(0,len(lats)*dt,dt)
		# Plot
		days = [0,1,2,3,4,5]
		for ii in range(len(days)):
			day = days[ii]
			pl.figure(ii+1)
			for t in range(len(starts)):
				if starts[t] <= day:
					timex = int((day-starts[t])/dt)
					pl.plot(Xs[t,:timex+1],Ys[t,:timex+1],'LightSlateGray',linewidth=1.5)
					pl.plot(Xs[t,:timex+1][::4],Ys[t,:timex+1][::4],'o',color='LightSlateGray',mec='LightSlateGrey',markersize=7)
			pl.title(day)
			self.proj.m.drawcoastlines()
			self.proj.m.drawparallels([70],latmax=90)
			# Focused
#			pl.ylim(0e06,5e06)
#			pl.xlim(1e06,6e06)
		pl.show()

	def plotTracks(self,lats,lons):
		# Centroids
		Xs,Ys = [],[]
		color = ['r','b']
		for t in range(len(lats)):
			xs,ys = self.proj.m(lons[t],lats[t])
			xs,ys = xs.mean(axis=0)[:],ys.mean(axis=0)[:]
			Xs.append(xs)
			Ys.append(ys)
		Xs,Ys = np.array(Xs),np.array(Ys)
		# Plot
		for i in range(len(Xs)):
			pl.plot(Xs[i],Ys[i],'r')
		self.proj.m.drawcoastlines()
		self.proj.m.drawparallels([60,70,80],latmax=90)
		pl.show()

	def densityToSeason(self,nn,dates,YearRange=(1990,2012),Season='DJ'):

		years,N,d,ind = datesToSeason(dates,YearRange=YearRange,Season=Season)
		ns = []
		for i in range(len(ind)):
			n = np.zeros(nn.shape[-2:])
			for j in ind[i]:
				n = n + nn[j]
			ns.append(n)
		ns = np.array(ns)
		return years,ns

	def detrend2d(self,years,field):

		n0,n1,n2 = field.shape
		m,c,p   = np.zeros((n1,n2)),np.zeros((n1,n2)),np.zeros((n1,n2))
		for i in range(n1):
			for j in range(n2):
				slope, intercept, r_value, p_value, std_err = stats.linregress(years,field[:,i,j])
				m[i,j] = slope
				c[i,j] = intercept
				p[i,j] = p_value
		line  = np.array([m*year + c for year in years])
		field = field - line
		return field,10*m,p

	def datesAtMaxInten(self,G,Q,Dates):
		dates,lons = [],[]
		for i in range(len(Q)):
			ymax    = [max(Q[i][j]) for j in range(len(Q[i]))]
			yargmax = [np.argmax(Q[i][j]) for j in range(len(Q[i]))]
			k       = np.argmax(ymax)
			l       = yargmax[k]
			dates.append(Dates[i][k])
			lons.append(G[i][k][l])
		return dates,lons

	def filterInjections(self,G,Q,Dates,LON,LAT,P,Sector=(330,105),YearRange=(1949,2016)):
		dates,lons    = self.datesAtMaxInten(G,Q,Dates)
		inds ,c       = [],0
		g,q,d,ln,lt,p = [],[],[],[],[],[]
		for i in range(len(lons)):
			if Sector[0] > Sector[1]:
				if (Sector[0] <= lons[i] <= 360) or (0 <= lons[i] < Sector[1]):
					if YearRange[0] <= Dates[i][0][0] <= YearRange[1]:
						g.append(G[i])
						q.append(Q[i])
						d.append(Dates[i])
						ln.append(LON[i])
						lt.append(LAT[i])
						p.append(P[i])
                        elif Sector[0] < Sector[1]:
                                if (Sector[0] <= lons[i] < Sector[1]):
					if YearRange[0] <= Dates[i][0][0] <= YearRange[1]:
						g.append(G[i])
						q.append(Q[i])
						d.append(Dates[i])
						ln.append(LON[i])
						lt.append(LAT[i])
						p.append(P[i])
		return g,q,d,ln,lt,p

	def regressModel(self,Model='CMCC-CESM',Season='DJF',ExpType='rcp85',Sector=(0,360),Field=['flds','rlds']):
		# T2 clim
		if ExpType == 'rcp85':
			years = range(2075,2100+1,1)
			Dir   = '/mnt/climstorage/cian'
		if ExpType == 'historical':
			years = range(1980,2005+1,1)
			Dir   = '/qnap/cian/cmip'
		tds   = cmipDataServer(Field=Field[1],LevType='surface',Source=Model,ExpType=ExpType,DataFreq='mon')
		f     = np.array([self.proj(tds.getDataSnaps(Year=year,Season=Season).mean(axis=0),tds.lon,tds.lat) for year in years])
		fclim = f.mean(axis=0)
		# Trajectory density
		LON,LAT,P       = unpick('%s/scripts/intrusiontracks/%s_tracks.6hour.p' % (Dir,Model))
		G,Q,D           = unpick('%s/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Dir,Model))
		# Filter by Sector of origin
		G,Q,D,LON,LAT,P = self.filterInjections(G,Q,D,LON,LAT,P,Sector=Sector)
		nn,Ntot,Ttot    = self.density(LON,LAT)
		dates           = [d[0] for d in D]
		years,ns        = self.densityToSeason(nn,dates,YearRange=(years[0],years[-1]),Season=Season)
		# Detrend ns and f
		ns,nsm,nsp = self.detrend2d(years,ns)
		f,fm,fmp = self.detrend2d(years,f)
		# Regression with n
		n0,n1 = f.shape[-2:]
		R,P   = np.zeros((n0,n1)),np.zeros((n0,n1))
		for i in range(n0):
			for j in range(n1):
				slope, intercept, r_value, p_value, std_err = stats.linregress(ns[:,i,j],f[:,i,j])
				if math.isnan(slope):
					slope   = 0
					p_value = 1
				R[i,j],P[i,j] = slope,p_value
		return R,P,Ntot,fclim

	def regressSeason(self,Sectors=[(330,105),(105,280),(280,330)],Field=['flds','rlds']):
		mdNLin      = [[] for isi in Sectors]
		dNLin       = [[] for isi in Sectors]	
		Models      = [g[9:] for g in glob.glob('../rcp85/*') if (g[9:]!='CMCC-CESM') and (g[9:]!='inmcm4') and (g[9:]!='MIROC-ESM') and (g[9:]!='FGOALS-g2')]
		if Field[0] == 'flds': cseq,units = np.arange(-10,10+1,1),'W m$^{-2}$'
		if Field[0] == 'T2':   cseq,units = np.arange(-1,1+0.2,0.2),'K'
		for ii in range(len(Sectors)):
			Sector = Sectors[ii]
			print Sector
			index  = []
			for Model in Models:
				print Model
				dT,mdN,R1,R2,Pc,NN = [],[],[],[],[],[]
				# CMIP5 historical regression
				regfileh = '/mnt/climstorage/cian/scripts/regressions/%s/%s.%s.%s-%s.p' % (Field[0],'historical',Model,Sector[0],Sector[1])
				if os.path.isfile(regfileh) == False:
					Rh,Ph,Ntoth,fclimh = self.regressModel(Model=Model,Season='DJF',ExpType='historical',Sector=Sector,Field=Field)
					toPick([Rh,Ph,Ntoth,fclimh],regfileh)
				else:
					Rh,Ph,Ntoth,fclimh = unpick(regfileh)
                                # CMIP5 rcp85 regression
                                regfiler = '/mnt/climstorage/cian/scripts/regressions/%s/%s.%s.%s-%s.p' % (Field[0],'rcp85',Model,Sector[0],Sector[1])
                                if os.path.isfile(regfiler) == False:
                                        Rr,Pr,Ntotr,fclimr = self.regressModel(Model=Model,Season='DJF',ExpType='rcp85',Sector=Sector,Field=Field)
                                        toPick([Rr,Pr,Ntotr,fclimr],regfiler)
                                else:
                                        Rr,Pr,Ntotr,fclimr = unpick(regfiler)
				# Smooth
				Rh,x,y    = self.interp2d(Rh,self.x,self.y,6,kind='linear')
				Ph,x,y    = self.interp2d(Ph,self.x,self.y,6,kind='linear')
				Ntoth,x,y = self.interp2d(Ntoth/25.,self.x,self.y,6,kind='linear')
                                Rr,x,y    = self.interp2d(Rr,self.x,self.y,6,kind='linear')
                                Pr,x,y    = self.interp2d(Pr,self.x,self.y,6,kind='linear')
                                Ntotr,x,y = self.interp2d(Ntotr/26.,self.x,self.y,6,kind='linear')
				dN        = Ntotr - Ntoth
				dF,x,y    = self.interp2d(fclimr-fclimh,self.x,self.y,6,kind='linear')
				# Append to list
				dT.append(dF)
				mdN.append(dN*Rr)
				R1.append(Rh)
				R2.append(Rr)
				Pc.append(Ph)
				NN.append(dN)
				dNLin[ii].append(dN)
				mdNLin[ii].append(dN*Rr)
			# Model mean biases
			dT,R1,R2,mdN,NN,Pc = np.array(dT).mean(axis=0),np.array(R1).mean(axis=0),np.array(R2).mean(axis=0),\
					     np.array(mdN).mean(axis=0),np.array(NN).mean(axis=0),np.array(Pc).mean(axis=0)
			if Field[0] == 'flds': cseq,units = np.arange(-40,40+5,5),'W m$^{-2}$'
			if Field[0] == 'T2':   cseq,units = np.arange(-18,18+2,2),'K'
			self.plotMap(x,y,dT,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Actual bias [%s]' % (units),\
				     title='all',savename='figs/bias/%s/d%s/%s-%s/%s.dT.bias.pdf' % (Field[0],Field[0],Sector[0],Sector[1],'all'))
			if Field[0] == 'flds': cseq,units = np.arange(-8,8+1,1),'W m$^{-2}$'
			if Field[0] == 'T2':   cseq,units = np.arange(-4,4+0.5,0.5),'K'
			self.plotMap(x,y,mdN,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Predicted bias [%s]' % (units),\
				     title='all',savename='figs/bias/%s/mdN/%s-%s/%s.mdN.bias.pdf' % (Field[0],Sector[0],Sector[1],'all'),sxy=None)
			if Field[0] == 'flds': cseq,units = np.arange(-6,6+1,1),'W m$^{-2}$'
			if Field[0] == 'T2':   cseq,units = np.arange(-1,1+0.2,0.2),'K'
			self.plotMap(x,y,R1-R2,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Regression coefficient [%s intrusion$^{-1}$]' % (units),\
				     title='all',savename='figs/bias/%s/R/%s-%s/%s.dR.pdf' % (Field[0],Sector[0],Sector[1],'all'))
			if Field[0] == 'flds': cseq,units = np.arange(-4,4+0.5,0.5),'W m$^{-2}$'
			if Field[0] == 'T2':   cseq,units = np.arange(-1,1+0.2,0.2),'K'
			self.plotMap(x,y,R1,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Regression coefficient [%s intrusion$^{-1}$]' % (units),\
				     title='all',savename='figs/bias/%s/R/%s-%s/%s.R.pdf' % (Field[0],Sector[0],Sector[1],'all'),P=Pc)
			cseq = np.arange(-2.5,2.5+0.5,0.5)
			self.plotMap(x,y,NN,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
				     title='all',savename='figs/bias/N/%s-%s/%s.dN.pdf' % (Sector[0],Sector[1],'all'))
		dNLin   = np.array(dNLin).sum(axis=0).mean(axis=0)
		mdNLin  = np.array(mdNLin).sum(axis=0).mean(axis=0)
		if Field[0] == 'flds': cseq,units = np.arange(-8,8+1,1),'W m$^{-2}$'
		if Field[0] == 'T2':   cseq,units = np.arange(-2,2+0.25,0.25),'K'
		self.plotMap(x,y,mdNLin,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Predicted bias [%s]' % (units),\
			title='all',savename='figs/bias/%s/%s.lin.mdN.bias.pdf' % (Field[0],'all'))       
		cseq = np.arange(-2.5,2.5+0.5,0.5)
		self.plotMap(x,y,dNLin,cseq,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
			title='all',savename='figs/bias/%s/%s.lin.dN.bias.pdf' % (Field[0],'all'))

	def intStats(self):
		# CMIP models
		Mh,Mr  = [],[]
		prop   = []
		Models = [g[9:] for g in glob.glob('../rcp85/*')]
		for Model in Models:
			print Model
			Gr,Qr,Dr = unpick('../intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
			Gh,Qh,Dh = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
			prop.append( 100.*(len(Gr)-len(Gh))/len(Gh) )
			dr,dh    = np.array([len(Dr[i]) for i in range(len(Dr))]),np.array([len(Dh[i]) for i in range(len(Dh))])
			qr       = np.array([np.array([np.sum(Qr[i][j]) for j in range(len(Qr[i]))]).sum() for i in range(len(Qr))])
			qh       = np.array([np.array([np.sum(Qh[i][j]) for j in range(len(Qh[i]))]).sum() for i in range(len(Qh))])
			ir,ih    = qr/dr,qh/dh
			ndi,idn  = len(Gh)*(qr.mean()-qh.mean()),qh.mean()*(len(Gr)-len(Gh))
			de       = qr.sum()-qh.sum()
			Mh.append(de/qh.sum())
			Mr.append((ndi+idn)/qh.sum())
		print np.array(prop).mean(),np.array(prop).std()
                pl.plot(Mr,Mh,'k+',markersize=10,mew=1.5)
                for i in range(len(Mr)):
                        pl.annotate(Models[i],xy=(Mr[i],Mh[i]),xytext=(-10,-15),size=10,textcoords='offset points',ha='right',va='bottom')
		pl.show()

		"""
			print ir.mean(),ir.std()
			Hh,edges = np.histogram(qh,bins=range(0,100000+5000,5000),range=None,normed=True,weights=None,density=None)
			Hr,edges = np.histogram(qr,bins=range(0,100000+5000,5000),range=None,normed=True,weights=None,density=None)
			Mh.append(Hh)
			Mr.append(Hr)
                Mh    = np.array(Mh).mean(axis=0)
                Mr    = np.array(Mr).mean(axis=0)
                PP,MM = [],[]
                bins  = (edges[1:]+edges[0:-1])/2.
                pp,   = pl.plot(bins,Hh,'b',linewidth=2,alpha=0.7)
                PP.append(pp)
                MM.append('Historical: 1981-2005')
                pp,   = pl.plot(bins,Hr,'r',linewidth=2,alpha=0.7)
                PP.append(pp)
                MM.append('RCP85: 2075-2100')
                pl.xlim(bins[0],bins[-1])
                pl.xlabel('Duration')
                pl.ylabel('Normed frequency')
                lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':15},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
                pl.show()
		"""

	def lonPDF(self,norm=True):
		if norm == True:  char = 'norm'
		if norm == False: char = 'full'
		# CMIP models
		Mh,Mr  = [],[]
#		Models = [g[9:] for g in glob.glob('../rcp85/*')]
		Models = ['jimbob']
		months = [10,11,12,1]
		bins          = 70
		for Model in Models:
			print Model
#			Gr,Qr,Dr = unpick('../intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
#			Gh,Qh,Dh = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
			Gh0,Qh0,Dh0 = unpick('/qnap/cian/sheba/intrusionfiles/ERAInt_intrusions.%s.6x6hr.9deg.200.6dt.20.5.filtered.23steps.80N.0.4.-100.0.full.p' % ('ONDJ'))
			Gh,Qh,Dh = Gh0[133:291],Qh0[133:291],Dh0[133:291]
			Gr,Qr,Dr = Gh0[291:],Qh0[291:],Dh0[291:]
			Gh = [Gh[i] for i in range(len(Dh)) if (Dh[i][0][1] in months)]
			Qh = [Qh[i] for i in range(len(Dh)) if (Dh[i][0][1] in months)]
			Dh = [Dh[i] for i in range(len(Dh)) if (Dh[i][0][1] in months)]
                        Gr = [Gr[i] for i in range(len(Dr)) if (Dr[i][0][1] in months)]
                        Qr = [Qr[i] for i in range(len(Dr)) if (Dr[i][0][1] in months)]
                        Dr = [Dr[i] for i in range(len(Dr)) if (Dr[i][0][1] in months)]
			print len(Gh),len(Qh),len(Dh)
			print len(Gr),len(Qr),len(Dr)
			lonsh    = np.array([k for i in Gh for j in i for k in j])
			lonsr    = np.array([k for i in Gr for j in i for k in j])
                        Hh,edges = np.histogram(lonsh,bins=bins,range=None,normed=norm,weights=None,density=None)
                        Hr,edges = np.histogram(lonsr,bins=bins,range=None,normed=norm,weights=None,density=None)
#			Hh,edges = np.histogram(lonsh,bins=range(0,360+5,5),range=None,normed=norm,weights=None,density=None)
#			Hr,edges = np.histogram(lonsr,bins=range(0,360+5,5),range=None,normed=norm,weights=None,density=None)	
			Mh.append(Hh)
			Mr.append(Hr)
		Mh    = np.array(Mh).mean(axis=0)
		Mr    = np.array(Mr).mean(axis=0)
		PP,MM = [],[]
		bins  = (edges[1:]+edges[0:-1])/2.
		pp,   = pl.plot(bins,Hh,'b',linewidth=2,alpha=0.7)
		PP.append(pp)
		MM.append('Historical: 1981-2005')
		pp,   = pl.plot(bins,Hr,'r',linewidth=2,alpha=0.7)
		PP.append(pp)
		MM.append('RCP85: 2075-2100')
#		pl.xlim(0,359)
		pl.xlabel('Longitude')
		pl.ylabel('Normed frequency')
		lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':15},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('figs/intpdf.%s.pdf' % (char), format='pdf')
		pl.show()

	def plotDensity(self,Sector=(0,360)):
		# CMIP models
		N,Nrm,T,M,Nh,Nhrm,Th,Mh  = [],[],[],[],[],[],[],[]
		Models = [g[9:] for g in glob.glob('../rcp85/*')]
		for Model in Models:
			fname  = 'intrusiontracks/%s_%s_tracks.6hour.p' % (Model,self.fstr)
			fnameh = '/qnap/cian/cmip/scripts/intrusiontracks/%s_%s_tracks.6hour.p' % (Model,self.fstr)
			if os.path.isfile(fname) and os.path.isfile(fnameh):
				G,Q,D           = unpick('../intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
				LONf,LATf,Pf    = unpick('intrusiontracks/%s_%stracks.6hour.p' % (Model,''))
				LONb,LATb,Pb    = unpick('intrusiontracks/%s_%stracks.6hour.p' % (Model,'back_'))
				for i in range(len(LONf)):
					for j in range(len(LONf[i])):
						LONf[i][j] = np.append(LONb[i][j][:,::-1][:,20:],LONf[i][j][:,1:20],axis=1)
						LATf[i][j] = np.append(LATb[i][j][:,::-1][:,20:],LATf[i][j][:,1:20],axis=1)
						Pf[i][j] = np.append(np.array(Pb[i][j])[:,::-1][:,20:],np.array(Pf[i][j])[:,1:20],axis=1)
				LON,LAT,P       = LONf,LATf,Pf
				G,D,Q,LON,LAT,P = self.filterInjections(G,Q,D,LON,LAT,P,Sector=Sector)
				nn,Ntot,Ttot    = self.density(LON,LAT)
				N.append(Ntot)
				Nrm.append(Ntot/Ntot.sum())
				T.append(Ttot)
				M.append(Model)

				Gh,Qh,Dh              = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
				LONhf,LAThf,Phf       = unpick('/qnap/cian/cmip/scripts/intrusiontracks/%s_%stracks.6hour.p' % (Model,''))
				LONhb,LAThb,Phb       = unpick('/qnap/cian/cmip/scripts/intrusiontracks/%s_%stracks.6hour.p' % (Model,'back_'))
                                for i in range(len(LONhf)):
                                        for j in range(len(LONhf[i])):
                                                LONhf[i][j] = np.append(LONhb[i][j][:,::-1][:,20:],LONhf[i][j][:,1:20],axis=1)
						LAThf[i][j] = np.append(LAThb[i][j][:,::-1][:,20:],LAThf[i][j][:,1:20],axis=1)
						Phf[i][j] = np.append(np.array(Phb[i][j])[:,::-1][:,20:],np.array(Phf[i][j])[:,1:20],axis=1)
				LONh,LATh,Ph          = LONhf,LAThf,Phf
				Gh,Dh,Qh,LONh,LATh,Ph = self.filterInjections(Gh,Qh,Dh,LONh,LATh,Ph,Sector=Sector)
				nnh,Ntoth,Ttoth       = self.density(LONh,LATh)
				Nh.append(Ntoth)
				Nhrm.append(Ntoth/Ntoth.sum())
				Th.append(Ttoth)
				Mh.append(Model)
				print Model
			else:
				print '%s failed to complete...' % (Model)
		# Compute bias and scale
		sf      = 1./26
		sfrm    = 1
		N,Nrm   = sf*np.array(N),sfrm*np.array(Nrm)
		Nh,Nhrm = sf*np.array(Nh),sfrm*np.array(Nhrm)
		Nt      = N - Nh
		Ntrm    = Nrm - Nhrm
		T       = np.array(T)
		Th      = np.array(Th)
		# Plot
		cseqb = np.arange(-2.5,2.5+0.5,0.5)
		cseq  = np.arange(0,6+0.5,0.5)
		cseqbrm,cseqrm = np.arange(-0.0024,0.0024+0.0003,0.0003),np.arange(0,0.015,0.001)
		for i in range(len(M)):
			Model = M[i]
			# Plot Model density bias
			n,x,y  = self.interp2d(Nt[i],self.x,self.y,6,kind='cubic')
			self.plotMap(x,y,n,cseqb,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
				     title=Model,savename='figs/bias/N/%s-%s/%s.%s.dN.pdf' % (Sector[0],Sector[1],Model,self.fstr))
			# Plot Model density absolute
			n,x,y  = self.interp2d(N[i],self.x,self.y,6,kind='cubic')
			self.plotMap(x,y,n,cseq,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
				     title=Model,savename='figs/bias/N/%s-%s/%s.%s.N.pdf' % (Sector[0],Sector[1],Model,self.fstr),Tang=T[i])
		# Stipling for > 80% of models with same bias sign
		Nt_,Ntrm_ = [],[]
		for i in range(len(Nt)):
			n_,x_,y_ = self.interp2d(Nt[i],self.x,self.y,1.5,kind='cubic')
			nrm_,x_,y_ = self.interp2d(Ntrm[i],self.x,self.y,1.5,kind='cubic')
			Nt_.append(n_)
			Ntrm_.append(nrm_)
		Nt_,Ntrm_       = np.array(Nt_),np.array(Ntrm_)
		stipx,stipy     = stipling(Nt_,xx=None,yy=None,x=x_,y=y_,thresh=0.8)
		stipxrm,stipyrm = stipling(Ntrm_,xx=None,yy=None,x=x_,y=y_,thresh=0.8)
		# Mean bias
		Nt,Ntrm = Nt.mean(axis=0),Ntrm.mean(axis=0)
		Nt,x,y  = self.interp2d(Nt,self.x,self.y,6,kind='cubic')
		Ntrm,x,y  = self.interp2d(Ntrm,self.x,self.y,6,kind='cubic')
		self.plotMap(x,y,Nt,cseqb,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
			     title='all',savename='figs/bias/N/%s-%s/%s.%s.dN.pdf' % (Sector[0],Sector[1],'all',self.fstr),stip=(stipx,stipy))
		self.plotMap(x,y,Ntrm,cseqbrm,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
                             title='all',savename='figs/bias/N/%s-%s/%s.%s.norm.dN.pdf' % (Sector[0],Sector[1],'all',self.fstr),stip=(stipxrm,stipyrm))
		# Mean absolute N
		N      = N.mean(axis=0)
		T      = T.mean(axis=0)
		N,x,y  = self.interp2d(N,self.x,self.y,6,kind='cubic')
		self.plotMap(x,y,N,cseq,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
			     title='all',savename='figs/bias/N/%s-%s/%s.%s.N.pdf' % (Sector[0],Sector[1],'all',self.fstr),Tang=T)
                Nh      = Nh.mean(axis=0)
                Th      = Th.mean(axis=0)
                Nh,x,y  = self.interp2d(Nh,self.x,self.y,6,kind='cubic')
                self.plotMap(x,y,Nh,cseq,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$ DJF}$^{-1}$',\
                             title='all',savename='figs/bias/N/%s-%s/%s.%s.hist.N.pdf' % (Sector[0],Sector[1],'all',self.fstr),Tang=Th)

	def plotMap(self,x,y,field,cseq,cmap,extend,cbarlabel,title,savename=None,P=None,sxy=None,Tang=None,stip=None):
		cf = pl.contourf(x,y,field,cseq,cmap=cmap,extend=extend)
		if P != None:
			cl = pl.contour(x,y,P,levels=[0.20],colors='k',linewidths=2,alpha=0.8)
		if sxy != None:
			sx,sy = sxy
			pl.plot(sx,sy,'k.',alpha=0.3)
		if Tang != None:
			#sx,x,y   = self.interp2d(self.proj.x,self.x,self.y,1,kind='linear')
			#sy,x,y   = self.interp2d(self.proj.y,self.x,self.y,1,kind='linear')
			#u,x,y    = self.interp2d(Tang[:,:,0],self.x,self.y,1,kind='linear')
			#v,x,y    = self.interp2d(Tang[:,:,1],self.x,self.y,1,kind='linear')
			sx,sy,u,v = self.proj.x,self.proj.y,Tang[:,:,0],Tang[:,:,1]
			Q = pl.quiver(sx,sy,u,v,units='inches',scale=2,\
				scale_units='inches',headwidth=3,headlength=5,headaxislength=4.5,pivot='tail')
			qk = pl.quiverkey(Q, 0.2, 1.02, 1, '%s%s' % (100,'%'), labelpos='W',fontproperties={'weight': 'bold'})
		if stip != None:
			pl.plot(stip[0],stip[1],'g.',markersize=8,alpha=0.5)
		cbar   = pl.colorbar(cf)
		cbar.set_label(cbarlabel)
		pl.title(title)
		#self.proj.m.drawcoastlines(color='0.4')
		drawCoastlinesNoRivers(self.proj.m,color='0.4')
		self.proj.m.drawparallels([70,80],latmax=90)
		if savename != None:
			pl.savefig(savename)
			pl.close()
		else:
			pl.show()

	def makeMonth(self,Year,Month,steps,lon0s,x0s,y0s,p0s,Ndays,dt,fname):
		if self.ff ==  1: ffj0 = None
		if self.ff == -1: ffj0 = -1
		# Make month file
		datelist  = self.uds.getDateList(Year,Month,Day=1)[::ffj0]	# Initialisation times (reverse if back trajectories)
		datelist_ = [self.uds.getDate( tr.uds.getHours(*datelist[-1]) + self.ff*6*i ) for i in range(1,Ndays*4+1,1)]
		hourlist  = [self.uds.getHours(*date) for date in datelist]
		u,v,w,pv  = self.data(datelist+datelist_)
		TX,TY,TP,TPV = [],[],[],[]
		for i in range(len(datelist)):
			print datelist[i]
		        X,Y,P,PV = self.tracksN(x0s,y0s,p0s,u[i:Ndays*4+i+1,:,:,:],v[i:Ndays*4+i+1,:,:,:],w[i:Ndays*4+i+1,:,:,:],pv[i:Ndays*4+i+1,:,:,:],dt=dt)
		        TX.append(X)
		        TY.append(Y)
		        TP.append(P)
			TPV.append(PV)
		TX,TY,TP,TPV = np.array(TX)[::ffj0,:,:],np.array(TY)[::ffj0,:,:],np.array(TP)[::ffj0,:,:],np.array(TPV)[::ffj0,:,:]
		lons,lats    = self.proj.m(TX,TY,inverse=True)
		Ps           = TP
		PVs	     = TPV
		trajNetCDF(hourlist[::ffj0],lons,lats,Ps,PVs,steps,lon0s,fname)

	def makeFiles(self,Years,Months,Ndays,dt,blat0=70,area_integration=True):
		# Make directory for saving files if it does not exist
		os.system('mkdir -p /mnt/climstorage/cian/Tracks/%s/%s/%s' % (self.Source,self.fstr,blat0)) 
		# Initialisation of particles (100 points north of 80N or 360 points along 70N)
		if area_integration:
			proj_tmp    = LambertProjector(boundinglat=80,resolution=200.)
			lon0s,lat0s = proj_tmp.lon.reshape(-1),proj_tmp.lat.reshape(-1)
		else:
			lon0s,lat0s = np.arange(360),[blat0 for kk in range(360)]
		steps       = self.ff*np.arange(24.*Ndays/dt + 1)*dt
		x0s,y0s     = self.proj.m(lon0s,lat0s)
		p0s         = [800 for kk in range(len(x0s))]
		# Loop trhough years and months
		for Year in Years:
			for Month in Months:
				#try:
					fname = '/mnt/climstorage/cian/Tracks/%s/%s/%s/%s_%s_%02d.nc' % (self.Source,self.fstr,blat0,self.fstr,Year,Month)	
					print fname
					if not os.path.isfile(fname):
						print 'Computing %s trajectories for %s %s ...' % (self.fstr,Year,Month)
						self.makeMonth(Year,Month,steps,lon0s,x0s,y0s,p0s,Ndays,dt,fname)
					else:
						print 'File %s already exists' % (fname)
				#except:
				#	pass

if __name__ == "__main__":

	from TrackDataServer import DataServer as TrackDataServer

#	td = TrackDataServer(Type='fwrd',steps=25)

	"""
	# Parameters
	YearRange = (int(sys.argv[1]),int(sys.argv[2]))
	Season    = str(sys.argv[3])
	Source    = str(sys.argv[4])
	ff        = int(sys.argv[5])
	G,Q,D     = unpick('/mnt/climstorage/cian/intrusions/%s_intrusions.%s-%s.%s.6x6hr.9deg.200.6dt.20.5.p' % (Source,YearRange[0],YearRange[1],Season))
	"""
	# Trajectory class setup
	tr = Trajectories(Source='ERAInt',dT=0.25,ff=1)

	tr.makeFiles(range(2018,2018+1,1),[3],Ndays=6,dt=6,blat0=70,area_integration=False)
#	tr.makeFiles(range(1900,1919+1,1),range(1,12+1,1),10,6,60,False)
#	tr.allField()

#	fname = 'alongtrack/%s.%s.p' % (Source,tr.fstr)
#	print fname
#	if os.path.isfile(fname)==False:
#		dh,dr = tr.fieldAlongTracks(Model=Source,Sector=(40,60))
#		toPick([dh,dr],fname)
#	else:
#		dh,dr = unpick(fname)
#	pl.plot(dh,'b')
#	pl.plot(dr,'r')
#	pl.show()

#	tr.intStats()
#	tr.lonPDF(norm=True)

	"""
	# Make trajectories
	fname = '/mnt/climstorage/intrusiontracks/%s_%s_tracks.%s-%s.%s.6hour.p' % (Source,tr.fstr,YearRange[0],YearRange[1],Season)
	print fname
	if os.path.isfile(fname) == False:
		LON,LAT,P     = [],[],[]
		for jj in range(len(D)):
			print '%s of %s ...' % (jj+1,len(D))
			lons,lats,ps = tr.intrusionTracks(G[jj],D[jj],dt=6,Ndays=10)
			LON.append(lons)
			LAT.append(lats)
			P.append(ps)
		toPick([LON,LAT,P],fname)
	else:
		print '%s already exists' % (fname)
		LON,LAT,P = unpick(fname)
	"""

#	for i in range(len(LAT)):
#		tr.plotTracks(LAT[i],LON[i])

#	tr.regressSeason(Sectors=[(280,330),(330,105),(105,280)],Field=['T2','tas'])
#	tr.regressSeason(Sectors=[(280,330),(330,105),(105,280)],Field=['flds','rlds'])
#	tr.plotDensity(Sector=(0,360))
