from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from UnPickle import *
from toPick import *
from scipy import interpolate
from scipy import ndimage as nd
from stipling import *

import glob,os
import numpy as np
import matplotlib.pyplot as pl

class VertComp:

	def __init__(	self,
			N        = 4,
			LatRange = (70,80),
			LonRange = (-180,180),
			Sector   = (330,105)	):

		# Attributes
		self.N        = N
		self.LatRange = LatRange
		self.LonRange = LonRange
		self.Sector   = Sector
		self.dt       = 24.
		self.plevs    = np.arange(50,1000+25,25)
#		self.days     = np.arange(-N,N+0.25,0.25)
		self.days     = np.arange(-N,N+1,1)
		# LambertProjector
#		self.proj      = LambertProjector(boundinglat=70,resolution=200.)
		self.proj      = LambertProjector(boundinglat=70,resolution=50.)
		self.x,self.y  = self.proj.x[0,:],self.proj.y[:,0]
		# Mask
		box = [LatRange,LonRange]
		if box[1][0] >= box[1][1]: lonA,latA = (self.proj.lon>box[1][1]) & (self.proj.lon<box[1][0]),(self.proj.lat>box[0][1]) | (self.proj.lat<box[0][0])
		if box[1][0] <  box[1][1]: lonA,latA = (self.proj.lon>box[1][1]) | (self.proj.lon<box[1][0]),(self.proj.lat>box[0][1]) | (self.proj.lat<box[0][0])
                self.mask = (lonA | latA)
#		x,y = np.ma.masked_array(self.proj.x,mask=self.mask),np.ma.masked_array(self.proj.y,mask=self.mask)
#		pl.plot(x,y,'k.')
#		self.proj.m.drawcoastlines()
#		self.proj.m.drawparallels([70,80],latmax=90)
#		pl.show()

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

	def intStats(self):
		Models  = [g[14:] for g in glob.glob('../historical/*') if g[14:]!='inmcm4']
		dur,flu = [],[]
		for Model in Models:
			print Model
			# Open trajectory and injection files
			LON,LAT,P = unpick('intrusiontracks/%s_tracks.6hour.p' % (Model))
			G,Q,Dates = unpick('../intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
			G,Q,Dates,LON,LAT,P = self.filterInjections(G,Q,Dates,LON,LAT,P)
			dur.append(np.array([len(i) for i in Dates]).mean())
			flu.append(np.array([np.mean(Q[l][k]) for l in range(len(Q)) for k in range(len(Q[l]))]).mean())
		dur = np.array(dur)
		flu = np.array(flu)
		print dur.mean(),flu.mean()

	def filterInjections(self,G,Q,Dates,LON,LAT,P):
		dates,lons    = self.datesAtMaxInten(G,Q,Dates)	
		inds ,c       = [],0
		g,q,d,ln,lt,p = [],[],[],[],[],[]
		for i in range(len(lons)):
			if self.Sector[0] > self.Sector[1]:
				if (self.Sector[0] <= lons[i] <= 360) or (0 <= lons[i] < self.Sector[1]):
					g.append(G[i])
					q.append(Q[i])
					d.append(Dates[i])
					ln.append(LON[i])
					lt.append(LAT[i])
					p.append(P[i])
			elif self.Sector[0] < self.Sector[1]:
				if (self.Sector[0] <= lons[i] < self.Sector[1]):
					g.append(G[i])
					q.append(Q[i])
					d.append(Dates[i])
					ln.append(LON[i])
					lt.append(LAT[i])
					p.append(P[i])
		return g,q,d,ln,lt,p

        def tracksNotMasked(self,Model='CMCC-CESM'):
                # Returns (time)x(coordinate) array of all trajectory points in unmasked region
                # If closest gridpoint in unmasked region then use trajectory point

		# DataServer for model
		if Model != 'ERAInt': qds = cmipDataServer(Field='ta',LevType='plev',Source=Model,ExpType='historical',DataFreq='day')
		if Model == 'ERAInt': qds = reDataServer(Field='T2',LevType='surface_analysis')
		# Open trajectory and injection files
		LON,LAT,P = unpick('/qnap/cian/cmip/scripts/intrusiontracks/%s_tracks.240Tg.6hour.p' % (Model))
		G,Q,Dates = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
		print len(LON),len(G)
		# Filter injections to self.Sector
		G,Q,Dates,LON,LAT,P = self.filterInjections(G,Q,Dates,LON,LAT,P)
		# Find trajectory points inside unmasked region
                time,points = [],[]
                for i in range(len(LON)):
                        for j in range(len(LON[i])):
				xs,ys = self.proj.m(LON[i][j],LAT[i][j])
				xs,ys = xs.mean(axis=0)[:21][::4],ys.mean(axis=0)[:21][::4]
                                for k in range(len(xs)):
                                        x,y   = xs[k],ys[k]
                                        ix,iy = np.argmin(np.abs(self.x-x)),np.argmin(np.abs(self.y-y))	
                                        if self.mask[iy,ix] == False:
#                                               hour = self.ds.getHours(*Dates[i][j]) + k*self.dt 
#						time.append(hour)	
						day  = qds.getDays(*Dates[i][j]) + k	
#						day  = qds.getDays(Dates[i][j][0],Dates[i][j][1],Dates[i][j][2],12) + k
						time.append(day)
						points.append([iy,ix])

#						pl.plot(self.proj.x,self.proj.y,'k.')
#						pl.plot(xs,ys,'r-')
#						pl.plot(xs,ys,'ro')
#						pl.plot(x,y,'bo')
#						self.proj.m.drawcoastlines()
#						self.proj.m.drawparallels([70,80],latmax=90)
#						pl.plot(self.proj.x[iy,ix],self.proj.y[iy,ix],'go')
#						pl.show()
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

 	def interpolateND(self,field,xold,xnew,axis,kind='linear'):
		# field in (time)x(lev)x(lon)
		f     = interpolate.interp1d(xold,field,kind=kind,axis=axis,bounds_error=True)
		field = f(xnew)
		return field

	def getClim2d(self,Model='CMCC-CESM',Field=['ta','T']):
		years = range(1981,2005+1,1)
#		filename = 'compfiles/%s/%s/%s.%s.%s-%s.clim.p' % (Field[1],Model,Model,Field[1],years[0],years[-1])
		filename = 'compfiles/%s/%s/%s.%s.%s-%s.%s-%sN.%s-%sE.clim.p' % \
			   (Field[1],Model,Model,Field[1],years[0],years[-1],self.LatRange[0],self.LatRange[1],self.LonRange[0],self.LonRange[1])
		print filename
		if os.path.isfile(filename) == False:
			C = []
			if Model != 'ERAInt':
				qds   = cmipDataServer(Field=Field[0],LevType='plev',Source=Model,ExpType='historical',DataFreq='day')
			if Model == 'ERAInt':
				qds   = reDataServer(Field=Field[1],LevType='plev')
			plevs = qds.lev
			for year in years:
				print year
				clim = qds.getDataSnaps(Year=year,Season='DJF',dailymean=True)
				if (type(clim)==np.ma.core.MaskedArray) and (clim.mask.shape == clim.data.shape):
					clim = self.fill(clim.data,invalid=clim.mask)
				clim = clim.mean(axis=0)
				clim = self.proj(clim,qds.lon,qds.lat)
				C.append(clim)
			C    = np.array(C).mean(axis=0)
			mask = np.array([self.mask for i in plevs])
			C    = np.ma.masked_array(C,mask=mask).mean(axis=1).mean(axis=1)
			toPick([C,plevs],filename)
		else:
			C,plevs = unpick(filename)
		return C,plevs

	def getClim1d(self,Model='CMCC-CESM',Field=['tas','T2']):
		years    = range(1981,2005+1,1)
		if Model != 'ERAInt':
			qds  = cmipDataServer(Field=Field[0],LevType='surface',Source=Model,ExpType='historical',DataFreq='mon')
			clim = np.array([qds.getDataSnaps(Year=year,Season='DJF').mean(axis=0) for year in years]).mean(axis=0)
		if Model == 'ERAInt':
			qds  = MMDataServer(Field=Field[1])
			clim = np.array([qds.getSeason(Year=year,Season='DJF') for year in years]).mean(axis=0)
		clim = self.proj(clim,qds.lon,qds.lat)
		clim = np.ma.masked_array(clim,mask=self.mask).mean(axis=0).mean(axis=0)
		return clim

        def getComp1d(self,hours,points,Model='CMCC-CESM',Field=['tas','T2']):
                filename = 'compfiles/%s/%s/%s.%s.%s-%sE.%sdays.comp.p' % (Field[1],Model,Model,Field[1],self.Sector[0],self.Sector[1],self.N)
                print filename          
                if os.path.isfile(filename) == False:
			# DataServer for model
			if Model != 'ERAInt':
				qds = cmipDataServer(Field=Field[0],LevType='surface',Source=Model,ExpType='historical',DataFreq='day')
			if Model == 'ERAInt':
				try:
					qds = reDataServer(Field=Field[1],LevType='surface_analysis')
				except:
					qds = reDataServer(Field=Field[1],LevType='surface_forecast')
			# Make the date list for the composite element
			days  = np.arange(-self.N,self.N+1,1)
			# Loop through all times and points
			comp = []
			for i in range(len(hours)):
				print '%s of %s ...' % (i+1,len(hours))
				dates = [self.ds.getDate(hours[i]+day*24) for day in days]
				# Data array
				data  = qds.getDataSnaps(datelist=dates,dailymean=True)
				if (type(data)==np.ma.core.MaskedArray) and (data.mask.shape == data.data.shape):
				        data = self.fill(data.data,invalid=data.mask)
				data = self.proj(data,qds.lon,qds.lat)
				# Extract the points    
				for iy,ix in points[i]:
					comp.append(data[:,iy,ix])
			N    = len(comp)
			comp = np.array(comp).mean(axis=0)
			toPick([days,comp,N],filename)
                else:
                        days,comp,N = unpick(filename)
                return days,comp,N

        def compAllModels1d(self):
		# Loop trhough Models
		dT,dF,CRD,CRU,CR,CS,CL = [],[],[],[],[],[],[]
		climT,climF            = [],[]

		# Reanalysis
		time,points  = self.tracksNotMasked(Model='ERAInt')
		hours,points = self.sortTimePoint(time,points)
		days,compRD,N = self.getComp1d(hours,points,Model='ERAInt',Field=['rlds','flds'])
		days,compR,N = self.getComp1d(hours,points,Model='ERAInt',Field=['rls','fls'])
		compRU       = compRD-compR
		days,compS,N = self.getComp1d(hours,points,Model='ERAInt',Field=['hfss','sshf'])
		days,compL,N = self.getComp1d(hours,points,Model='ERAInt',Field=['hfls','slhf'])
		CRD.append(compRD)#-compRD[0])
		CRU.append(compRU)
		CR.append(compR)#-compR[0])
		CS.append(compS)#-compS[0])
		CL.append(compL)#-compL[0])
		# Loop through models
		nModels = ['inmcm4','CSIRO-Mk-6-0','EC-Earth']
		Models  = [g[14:] for g in glob.glob('../historical/*') if g[14:] not in nModels]
		for Model in Models:
			print Model
			time,points   = self.tracksNotMasked(Model=Model)
			hours,points  = self.sortTimePoint(time,points)
			days,compRD,N = self.getComp1d(hours,points,Model=Model,Field=['rlds','flds'])
			days,compRU,N = self.getComp1d(hours,points,Model=Model,Field=['rlus','flus'])
			compRD,compRU = np.abs(compRD),np.abs(compRU)
			compR         = compRD-compRU
			days,compS,N  = self.getComp1d(hours,points,Model=Model,Field=['hfss','sshf'])
			days,compL,N  = self.getComp1d(hours,points,Model=Model,Field=['hfls','slhf'])
			x             = np.argmin((days-0)**2)
			CRD.append(compRD)#-compRD[0])
			CRU.append(compRU)
			CR.append(compR)#-compR[0])#-self.getClim1d(Model=Model,Field=['tas','T2']))
			CS.append(compS)#-compS[0])#-self.getClim1d(Model=Model,Field=['rlds','flds']))
			CL.append(compL)#-compL[0])
		Models = ['ERAInt']+Models
#			dT.append(compT[x]-compT[0])
#			dF.append(compF[x]-compF[0])
		# Regression of dF and dT
#		m = dF[0]/dT[0]
#		pl.figure(1)
#		pl.plot(dT,dF,'k+',markersize=10,mew=1.5)
#		for model,x0,y0 in zip(Models,dT,dF):
#			pl.annotate(model,xy=(x0, y0),xytext=(-10,-15),size=8,textcoords='offset points',ha='right',va='bottom',\
#			arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
#		pl.plot([0,8],[0,m*8],'k--',linewidth=1)
#		pl.xlabel('dT$_{SURF}$ [K]')
#		pl.ylabel('dF$_{SURF}$ [W m$^{-2}$]')
#		pl.xlim(1,7.5)
#		pl.ylim(10,35)
		# Composite
		CRD,CRU,CR,CS,CL = np.array(CRD),np.array(CRU),np.array(CR),np.array(CS),np.array(CL)
		CRD,CRU,CR,CS,CL = self.interpolateND(CRD,days,self.days,axis=1,kind='cubic'),self.interpolateND(CRU,days,self.days,axis=1,kind='cubic'),\
					self.interpolateND(CR,days,self.days,axis=1,kind='cubic'),self.interpolateND(CS,days,self.days,axis=1,kind='cubic'),\
					self.interpolateND(CL,days,self.days,axis=1,kind='cubic')
		CRDmean,CRUmean,CRmean,CSmean,CLmean = CRD[1:].mean(axis=0),CRU[1:].mean(axis=0),CR[1:].mean(axis=0),-1*CS[1:].mean(axis=0),-1*CL[1:].mean(axis=0)

		pl.figure(1)
		PP,MM = [],[]
		pp,   = pl.plot(self.days,CRDmean,'brown',linewidth=2)
		PP.append(pp)
		MM.append('Longwave down')
		pl.plot(self.days,CRD[0],'brown',linewidth=1,linestyle='--')

                pp, = pl.plot(self.days,CRUmean,'g',linewidth=2)
                PP.append(pp)
                MM.append('Longwave up')
                pl.plot(self.days,CRU[0],'g',linewidth=1,linestyle='--')

		lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.setp(lg.get_title(),fontsize=10)
		pl.xlim(-4,4)
		pl.ylabel('Energy flux [W m$^{-2}$]')
		pl.xlabel('Days')

		pl.figure(2)
		PP,MM = [],[]
		pp, = pl.plot(self.days,CRmean,'r',linewidth=2)
                PP.append(pp)
                MM.append('Net longwave')
		pl.plot(self.days,CR[0],'r--',linewidth=1)

		pp, = pl.plot(self.days,CSmean,'y',linewidth=2)
                PP.append(pp)
                MM.append('Sensible')
		pl.plot(self.days,CS[0],'y--',linewidth=1)

		pp, = pl.plot(self.days,CLmean,'b',linewidth=2)
                PP.append(pp)
                MM.append('Latent')
		pl.plot(self.days,CL[0],'b--',linewidth=1)

                pp, = pl.plot(self.days,CRmean+CSmean+CLmean,'k',linewidth=2)
                PP.append(pp)
                MM.append('Total')
                pl.plot(self.days,CR[0]+CS[0]+CL[0],'k--',linewidth=1)

		lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.setp(lg.get_title(),fontsize=10)
		pl.xlim(-4,4)
		pl.ylabel('Energy flux [W m$^{-2}$]')
		pl.xlabel('Days')
		pl.show()	

	def compAllModels1dPlev(self):
                # Loop trhough Models
                dT,dF  = [],[]
                Models = ['ERAInt']+[g[14:] for g in glob.glob('../historical/*') if ((g[14:] != 'inmcm4') and (g[14:] != 'CSIRO-Mk3-6-0') and (g[14:] != 'FGOALS-g2'))]
                for Model in Models:
			print Model
			time,points       = self.tracksNotMasked(Model=Model)
			hours,points      = self.sortTimePoint(time,points)	
			days,compF,N      = self.getComp1d(hours,points,Model=Model,Field=['rlds','flds'])
			days,plevs,comp,N = unpick('compfiles/q/%s/%s.q.0-360E.4days.comp.p' % (Model,Model))
			comp              = 1000.*self.interpolateND(comp,plevs,self.plevs,axis=1)
			x,z               = np.argmin((days-0)**2),np.argmin((self.plevs-850)**2)
			dT.append(comp[x,z]- comp[0,z])
			dF.append(compF[x] - compF[0])
		# Regression of dF and dT
		m = dF[0]/dT[0]
		pl.figure(1)
		pl.plot(dT,dF,'k+',markersize=10,mew=1.5)
		for model,x0,y0 in zip(Models,dT,dF):
		        pl.annotate(model,xy=(x0, y0),xytext=(-10,-15),size=8,textcoords='offset points',ha='right',va='bottom',\
		        arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
		pl.plot([0,8],[0,m*8],'k--',linewidth=1)
		pl.xlabel('dq$_{850hPa}$ [g kg$^{-1}$]')
		pl.ylabel('dF$_{SURF}$ [W m$^{-2}$]')
#		pl.xlim(1,4.5)
		pl.xlim(0.15,0.35)
		pl.ylim(10,35)
		pl.show()

	def getComp2d(self,hours,points,Model='CMCC-CESM',Field=['hus','q']):
#		filename = 'compfiles/%s/%s/%s.%s.%s-%sE.%sdays.comp.p' % (Field[1],Model,Model,Field[1],self.Sector[0],self.Sector[1],self.N)
		filename = 'compfiles/%s/%s/%s.%s.%s-%sE.%s-%sN.%s-%sE.%sdays.comp.p' % \
			   (Field[1],Model,Model,Field[1],self.Sector[0],self.Sector[1],self.LatRange[0],self.LatRange[1],self.LonRange[0],self.LonRange[1],self.N)
		print filename
		if os.path.isfile(filename) == False:
			# DataServer for model
			if Model != 'ERAInt':
				qds   = cmipDataServer(Field=Field[0],LevType='plev',Source=Model,ExpType='historical',DataFreq='day')
			if Model == 'ERAInt':
				qds = reDataServer(Field=Field[1],LevType='plev')
#				ids = reDataServer(Field='fls',LevType='surface_forecast')
			plevs = qds.lev
			# Make the date list for the composite element
			days  = np.arange(-self.N,self.N+1,1)
			i0    = self.N
			# Loop through all times and points
			comp = []
			for i in range(len(hours)):
				print '%s of %s ...' % (i+1,len(hours))
#				dates = [self.ds.getDate(hours[i]+day*24) for day in days]
				dates = [qds.getDate(hours[i]+day,nd=24) for day in days]
				print dates
				# Data arrays
#				indata = self.proj(ids.snapshot(*dates[i0]),ids.lon,ids.lat)
				data   = qds.getDataSnaps(datelist=dates,dailymean=True)
				if (type(data)==np.ma.core.MaskedArray) and (data.mask.shape == data.data.shape):
					data = self.fill(data.data,invalid=data.mask)
				data = self.proj(data,qds.lon,qds.lat)
				# Extract the points	
				for iy,ix in points[i]:
#					if indata[iy,ix] > -30:
#						self.plot(days,plevs,data[:,:,iy,ix],20,cmap=pl.cm.RdBu_r,Model=Model,Field=Field[1],stip=None)
						comp.append(data[:,:,iy,ix])
			N    = len(comp)
			comp = np.array(comp).mean(axis=0)
			toPick([days,plevs,comp,N],filename)
		else:
			days,plevs,comp,N = unpick(filename)
		return days,plevs,comp,N

	def fill(self,data,invalid=None):
		if invalid is None: invalid = np.isnan(data)
		ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
		return data[tuple(ind)]

	def compAllModels2d(self,Field=['hus','q']):
		# Contour levels
		if Field[1] == 'T': cseq,cseqb,cmap = np.arange(-4.5,4.5+0.5,0.5),np.arange(-1,1+0.2,0.2),pl.cm.RdBu_r
		if Field[1] == 'q': cseq,cseqb,cmap = np.arange(0,0.4+0.025,0.025),np.arange(-0.1,0.12,0.02),pl.cm.Blues
		# ERAInt
		print 'ERAInt'
		time,points       = self.tracksNotMasked(Model='ERAInt')
		hours,points      = self.sortTimePoint(time,points)
		days,plevs,comp,N = self.getComp2d(hours,points,Model='ERAInt',Field=Field)
		clim,plevs        = self.getClim2d(Model='ERAInt',Field=Field)
#		anom              = comp - comp[0,:][np.newaxis,:]
		anom              = comp - clim[np.newaxis,:]
		# Interpolate to new time and levels
		E    = self.interpolateND(anom,days,self.days,axis=0)
		E    = self.interpolateND(E,plevs,self.plevs,axis=1)
		self.plot(self.days,self.plevs,E,cseq,cmap=cmap,Model='ERAInt',Field=Field[1])
		# Loop through Models
		Models = [g[14:] for g in glob.glob('../historical/*') if (g[14:] != 'inmcm4')]# and (g[14:] != 'MIROC4h'))]
		C      = []	
		for Model in Models:
			print Model
			time,points       = self.tracksNotMasked(Model=Model)
			hours,points      = self.sortTimePoint(time,points)
			days,plevs,comp,N = self.getComp2d(hours,points,Model=Model,Field=Field)
			clim,plevs        = self.getClim2d(Model=Model,Field=Field)
			anom              = comp - clim[np.newaxis,:]
			# Interpolate to new time and levels
			anom            = self.interpolateND(anom,days,self.days,axis=0)
			anom            = self.interpolateND(anom,plevs,self.plevs,axis=1)
			self.plot(self.days,self.plevs,anom,cseq,cmap=cmap,Model=Model,Field=Field[1])
#			self.plot(self.days,self.plevs,anom-E,cseqb,cmap=pl.cm.RdBu_r,Model='%s.d%s' % (Model,Field[1]),Field=Field[1])
			C.append(anom)
		C    = np.array(C)
		stip = stipling(C,x=self.plevs,y=self.days,thresh=0.8)
		Cm   = C[1:].mean(axis=0)
		Cs   = C.std(axis=0)
		self.plot(self.days,self.plevs,Cm,cseq,cmap=cmap,Model='all',Field=Field[1],stip=stip)
#		self.plot(self.days,self.plevs,Cm-E,cseqb,cmap=pl.cm.RdBu_r,Model='all.diff',Field=Field[1])
#		self.plot(self.days,self.plevs,Cs,20,cmap=cmap,Model='all.std',Field=Field[1])

	def plot(self,days,plevs,anom,cseq,cmap,Model,Field,stip=None):
		# Rotate array
		anom = np.rollaxis(anom,1,0)
		if Field == 'q': sf,clabel = 1000.,'Specific humidity anomaly [g kg$^{-1}$]'
		if Field == 'T': sf,clabel = 1.,'Temperature anomaly [K]'
		pl.show()
		cf   = pl.contourf(days,plevs,sf*anom,cseq,cmap=cmap,extend='both')
		cbar = pl.colorbar(cf)
		cbar.set_label(clabel)
                if stip != None: pl.plot(stip[1][::2,::2],stip[0][::2,::2],'k.',alpha=0.5)
		pl.title('%s: %s' % (Model,anom[-1,:].mean()))
		pl.ylabel('Pressure [hPa]')
		pl.xlabel('Days')
		pl.ylim(1000,50)
		pl.xlim(-self.N,self.N)
		pl.show()
#		pl.savefig('figs/comp/%s/%s-%s/%s.pdf' % (Field,self.Sector[0],self.Sector[1],Model),format='pdf')
#		pl.close()

	def plotStates(self):
		dTcold,dTwarm = [],[]
		lx            = np.argmin((self.plevs-850)**2)
		# ERAInt clear and opaque
		filename1            = 'compfiles/T/%s/%s.T.1981-2005.clim.p' % ('ERAInt','ERAInt')
		filename2            = 'compfiles/T/%s/%s.T.%s-%sE.%sdays.comp.p' % ('ERAInt','ERAInt',self.Sector[0],self.Sector[1],self.N)
		ERAclim,plevs        = unpick(filename1)
		days,plevs,ERAcomp,N = unpick(filename2)
		ERAclim              = self.interpolateND(ERAclim,plevs,self.plevs,axis=0,kind='linear')
		ERAcomp              = self.interpolateND(ERAcomp,days,self.days,axis=0,kind='linear')
		ERAcomp              = self.interpolateND(ERAcomp,plevs,self.plevs,axis=1,kind='linear')
		self.plot(self.days,self.plevs,ERAcomp,cseq=20,cmap=pl.cm.Reds,Model='ERAInt',Field='T')
		x                    = np.argmax(ERAcomp[:,lx])
		maxwarmERA           = ERAcomp[x,:]
		dTcold.append(ERAcomp[0,lx]-ERAcomp[0,-1])
		dTwarm.append(ERAcomp[x,lx]-ERAcomp[x,-1])
#		dTcold.append(ERAclim[lx]-ERAclim[-1])
#		dTwarm.append(maxwarmERA[lx]-maxwarmERA[-1])
		# CMIP5 models
		Models = [g[14:] for g in glob.glob('../historical/*') if ((g[14:] != 'inmcm4') and (g[14:] != 'EC-Earth'))]
		for Model in Models:
			filename1         = 'compfiles/T/%s/%s.T.1981-2005.clim.p' % (Model,Model)
			filename2         = 'compfiles/T/%s/%s.T.%s-%sE.%sdays.comp.p' % (Model,Model,self.Sector[0],self.Sector[1],self.N)
			clim,plevs        = unpick(filename1)
			days,plevs,comp,N = unpick(filename2)
			clim              = self.interpolateND(clim,plevs,self.plevs,axis=0,kind='linear')
			comp              = self.interpolateND(comp,days,self.days,axis=0,kind='linear')
			comp              = self.interpolateND(comp,plevs,self.plevs,axis=1,kind='linear')
#			self.plot(self.days,self.plevs,comp,cseq=20,cmap=pl.cm.Reds,Model=Model,Field='T')
			x                 = np.argmax(comp[:,lx])
			maxwarm           = comp[x,:]
			dTcold.append(comp[0,lx]-ERAcomp[0,-1])
			dTwarm.append(comp[x,lx]-ERAcomp[x,-1])
			# Plot profiles
#			dTcold.append(clim[lx]-clim[-1])
#			dTwarm.append(maxwarm[lx]-maxwarm[-1])
		Models = ['ERAInt']+Models
		pl.plot(dTcold,dTwarm,'k+',markersize=10,mew=1.5)
		for model,x0,y0 in zip(Models,dTcold,dTwarm):
			pl.annotate(model,xy=(x0, y0),xytext=(-10,-15),size=8,textcoords='offset points',ha='right',va='bottom',\
				    arrowprops=dict(arrowstyle='->',connectionstyle='arc3,rad=0'))
		pl.plot(np.arange(-2,8+1,1),np.arange(-2,8+1,1),'k--')
		pl.xlabel('dT cold [K]')
		pl.ylabel('dT warm [K]')
		pl.xlim(-2,8)
		pl.ylim(-2,8)
		pl.show()
		"""
			pl.plot(maxwarm,self.plevs,'r',linewidth=1)
			pl.plot(maxwarmERA,self.plevs,'r--',linewidth=1.5,alpha=0.5)
			pl.plot(clim,self.plevs,'b',linewidth=1)
			pl.plot(ERAclim,self.plevs,'b--',linewidth=1.5,alpha=0.5)
			pl.title(Model)
			pl.ylim(1000,650)
			pl.xlim(240,260)
			pl.savefig('/qnap/cian/cmip/scripts/figs/comp/T/profiles/%s.pdf' % (Model),format='pdf')
			pl.close()
		"""

if __name__ == "__main__":

	vc  = VertComp()
	t,p = vc.tracksNotMasked('ERAInt')
	t,p = vc.sortTimePoint(t,p)
	print t
	print p
