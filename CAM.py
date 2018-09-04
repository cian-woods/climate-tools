from TrackDataServer import DataServer as trDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from ReanalysisDataServer import DataServer as reDataServer
from LambertProjector import *
from netCDF4 import Dataset
from scipy.signal import butter, lfilter
from UnPickle import *
from toPick import *
from scipy import stats,interpolate
from datesToN import *

import matplotlib.pyplot as pl

class CAM:

	def __init__(	self,	
			Season     =       'DJF',
			blat       =          70,
			res        =         200,
			YearRange  = (1901,1919),
			YearRangeE = (1980,2000)	):
		# Attributes
		self.Season     = Season
		self.blat       = blat
		self.res        = res
		self.YearRange  = YearRange
		self.YearRangeE = YearRangeE
		self.years      = range(YearRange[0],YearRange[1]+1,1)
		self.yearsE     = range(YearRangeE[0],YearRangeE[1]+1,1)
		self.yn         = 1.*len(self.years)
		self.yEn        = 1.*len(self.yearsE)
		self.proj       = LambertProjector(boundinglat=self.blat,resolution=self.res)
		self.projO      = LambertProjector(boundinglat=30,resolution=self.res)

        def densityToSeason(self,nn,dates,YearRange):
		years,N,d,ind = datesToSeason(dates,YearRange=YearRange,Season=self.Season)
		ns = []
		for i in range(len(ind)):
		        n = np.zeros(nn.shape[-2:])
		        for j in ind[i]:
		                n = n + nn[j]
		        ns.append(n)
		ns = np.array(ns)
		return years,ns

	def intrusionDensity(self,lat0=70):
		# Thresholds
		#fE,f4,f8,f16 = 200,415,780,1205
		fE,f4,f8,f16 = 200,1130,1840,2560
		# DataServers for regressions
		#dsE        = MMDataServer(Field='Ts',Source='ERAInt',LevType='surface_analysis')
		ds4         = MMDataServer(Field='Ts',Source='CAM4xCO2')
		ds8         = MMDataServer(Field='Ts',Source='CAM8xCO2')
		ds16        = MMDataServer(Field='Ts',Source='CAM16xCO2')
		# TrackDataServer and injection file
		#trdsE      = trDataServer(Source='ERAInt',Type='fwrd',blat=lat0,steps=21)
		trds4       = trDataServer(Source='CAM4xCO2',Type='fwrd',blat=lat0,steps=21)
		trds8       = trDataServer(Source='CAM8xCO2',Type='fwrd',blat=lat0,steps=21)
		trds16      = trDataServer(Source='CAM16xCO2',Type='fwrd',blat=lat0,steps=21)
		#GE,QE,DE   = unpick('/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.%sN.filtered.p' % ( 'ERAInt',self.yearsE[0],self.yearsE[-1],fE,lat0))
		G4,Q4,D4    = unpick('/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.%sN.filtered.%sN.p' % ( 'CAM4xCO2',self.years[0],self.years[-1],f4,lat0,lat0+10))
		G8,Q8,D8    = unpick('/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.%sN.filtered.%sN.p' % ( 'CAM8xCO2',self.years[0],self.years[-1],f8,lat0,lat0+10))
		G16,Q16,D16 = unpick('/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.%sN.filtered.%sN.p' % ('CAM16xCO2',self.years[0],self.years[-1],f16,lat0,lat0+10))
		# Intrusion trajectories
		#LONE,LATE,PE   = trdsE.getIntrusionTrajectories(GE,DE)
		LON4,LAT4,P4    = trds4.getIntrusionTrajectories(G4,D4)
		LON8,LAT8,P8    = trds8.getIntrusionTrajectories(G8,D8)
		LON16,LAT16,P16 = trds16.getIntrusionTrajectories(G16,D16)
		# Surface fields
		#sE = self.proj(np.array([dsE.getSeason(Year=year,Season=self.Season) for year in self.yearsE]),dsE.lon,dsE.lat)
		s4 = self.proj(np.array([ds4.getSeason(Year=year,Season=self.Season) for year in self.years]),ds4.lon,ds4.lat)
		s8 = self.proj(np.array([ds8.getSeason(Year=year,Season=self.Season) for year in self.years]),ds8.lon,ds8.lat)
		s16 = self.proj(np.array([ds16.getSeason(Year=year,Season=self.Season) for year in self.years]),ds16.lon,ds16.lat)
		# Densities
		#nnE,NtotE,TtotE,x,y   = trds4.density(LONE,LATE,self.proj)
		nn4,Ntot4,Ttot4,x,y    = trds4.density(LON4,LAT4,self.proj)
		nn8,Ntot8,Ttot8,x,y    = trds4.density(LON8,LAT8,self.proj)
		nn16,Ntot16,Ttot16,x,y = trds4.density(LON16,LAT16,self.proj)
		# Seasonal densities
		#years,nsE  = self.densityToSeason( nnE,  [DE[i][0] for i in range(len(DE))],self.YearRangeE)
		years,ns4  = self.densityToSeason( nn4,  [D4[i][0] for i in range(len(D4))], self.YearRange)
		years,ns8  = self.densityToSeason( nn8,  [D8[i][0] for i in range(len(D8))], self.YearRange)
		years,ns16 = self.densityToSeason(nn16,[D16[i][0] for i in range(len(D16))], self.YearRange)
		# Regressions
		mE,m4,m8,m16 = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
				#slopeE , intercept, r_value, p_value, std_err = stats.linregress( nsE[:,i,j], sE[:,i,j])
				slope4 , intercept, r_value, p_value, std_err = stats.linregress( ns4[:,i,j], s4[:,i,j])
				slope8 , intercept, r_value, p_value, std_err = stats.linregress( ns8[:,i,j], s8[:,i,j])
				slope16, intercept, r_value, p_value, std_err = stats.linregress(ns16[:,i,j],s16[:,i,j])
		#		mE[i,j]  =  slopeE
				m4[i,j]  =  slope4
				m8[i,j]  =  slope8
				m16[i,j] = slope16
		# Smooth
		#NtotE,xn,yn  = trds4.interp2d(NtotE,x,y,6,kind='linear')
		Ntot4,xn,yn  = trds4.interp2d(Ntot4,x,y,6,kind='linear')
		Ntot8,xn,yn  = trds4.interp2d(Ntot8,x,y,6,kind='linear')
		Ntot16,xn,yn = trds4.interp2d(Ntot16,x,y,6,kind='linear')
		#mE,xn,yn     = trds4.interp2d(mE,x,y,6,kind='linear')
		m4,xn,yn     = trds4.interp2d(m4,x,y,6,kind='linear')
		m8,xn,yn     = trds4.interp2d(m8,x,y,6,kind='linear')
		m16,xn,yn    = trds4.interp2d(m16,x,y,6,kind='linear')
		# Plot
		cseqf  = np.arange(0,12+1,1)
		cseqa  = np.arange(-2,2+0.25,0.25)
		clabel = 'Number density {%sx%s km$^{2}$ %s}$^{-1}$' % (self.res,self.res,self.Season)
		#self.plotMap(1,xn,yn,NtotE/self.yEn,cseq=cseqf,cmap=pl.cm.OrRd,extend='max',cbarlabel=clabel,title='ERAInt',proj=self.proj,savename=None,P=None,sxy=None,Tang=TtotE,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/N.CAM4xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(2,xn,yn,Ntot4/self.yn,cseq=cseqf,cmap=pl.cm.OrRd,extend='max',cbarlabel=clabel,title='CAM4xCO2',proj=self.proj,savename=sname,P=None,sxy=None,Tang=Ttot4,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/N.CAM8xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(3,xn,yn,Ntot8/self.yn,cseq=cseqf,cmap=pl.cm.OrRd,extend='max',cbarlabel=clabel,title='CAM8xCO2',proj=self.proj,savename=sname,P=None,sxy=None,Tang=Ttot8,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/N.CAM16xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(4,xn,yn,Ntot16/self.yn,cseq=cseqf,cmap=pl.cm.OrRd,extend='max',cbarlabel=clabel,title='CAM16xCO2',proj=self.proj,savename=sname,P=None,sxy=None,Tang=Ttot16,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/N.CAM8-4xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(5,xn,yn,(Ntot8-Ntot4)/self.yn,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='8 - 4 xCO2',proj=self.proj,savename=sname,P=None,sxy=None,Tang=Ttot8-Ttot4,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/N.CAM16-8xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(6,xn,yn,(Ntot16-Ntot8)/self.yn,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='16 - 8 xCO2',proj=self.proj,savename=sname,P=None,sxy=None,Tang=Ttot16-Ttot8,stip=None)

                cseqa  = np.arange(-0.8,0.8+0.1,0.1)
                clabel = 'Regression coefficient [K intrusion$^{-1}$]'
		#self.plotMap(7,xn,yn,mE,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='ERAInt dT/dN',proj=self.proj,savename=None,P=None,sxy=None,Tang=None,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/m.CAM4xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(8,xn,yn,m4,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='CAM4xCO2 dT/dN',proj=self.proj,savename=sname,P=None,sxy=None,Tang=None,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/m.CAM8xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(9,xn,yn,m8,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='CAM8xCO2 dT/dN',proj=self.proj,savename=sname,P=None,sxy=None,Tang=None,stip=None)
		sname = '/mnt/climstorage/cian/scripts/figs/CAM/m.CAM16xCO2.%s-%s.%s.%sN.pdf' % (self.years[0],self.years[-1],self.Season,lat0)
		self.plotMap(10,xn,yn,m16,cseq=cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=clabel,title='CAM16xCO2 dT/dN',proj=self.proj,savename=sname,P=None,sxy=None,Tang=None,stip=None)

		pl.show()

	def plotMap(self,fign,x,y,field,cseq,cmap,extend,cbarlabel,title,proj,savename=None,P=None,sxy=None,Tang=None,stip=None):
		pl.figure(fign)
	        cf = pl.contourf(x,y,field,cseq,cmap=cmap,extend=extend)
	        if P != None:
	                cl = pl.contour(x,y,P,levels=[0.20],colors='k',linewidths=2,alpha=0.8)
	        if sxy != None:
	                sx,sy = sxy
	                pl.plot(sx,sy,'k.',alpha=0.3)
	        if Tang != None:
	                sx,sy,u,v = proj.x,proj.y,Tang[:,:,0],Tang[:,:,1]
	                Q = pl.quiver(sx,sy,u,v,units='inches',scale=2,\
	                        scale_units='inches',headwidth=3,headlength=5,headaxislength=4.5,pivot='tail')
	                qk = pl.quiverkey(Q, 0.2, 1.02, 1, '%s%s' % (100,'%'), labelpos='W',fontproperties={'weight': 'bold'})
	        if stip != None:
	                pl.plot(stip[0],stip[1],'g.',markersize=8,alpha=0.5)
	        cbar   = pl.colorbar(cf)
	        cbar.set_label(cbarlabel)
	        pl.title(title)
	        proj.m.drawparallels([70,80],latmax=90)
		#x0,y0 = proj.m(0,self.lat)
		#x1,y1 = proj.m(90,self.lat)
		#x2,y2 = proj.m(180,self.lat)
		#x3,y3 = proj.m(270,self.lat)
		#pl.ylim(y0,y2)
		#pl.xlim(x3,x1)
	        if savename != None:
	                pl.savefig(savename,format='pdf')

	def plotClim(self,Field):
		ds4  = MMDataServer(Field=Field,Source='%s' % ('CAM4xCO2'),LevRange=(500,500))
		ds8  = MMDataServer(Field=Field,Source='%s' % ('CAM8xCO2'),LevRange=(500,500))
		ds16 = MMDataServer(Field=Field,Source='%s' % ('CAM16xCO2'),LevRange=(500,500))
		s4   = np.array([ds4.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0)
                s8   = np.array([ds8.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0)
                s16  = np.array([ds16.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0)
		s4   = s4 - s4.mean(axis=-1)[:,np.newaxis]
		s8   = s8 - s8.mean(axis=-1)[:,np.newaxis]
		s16   = s16 - s16.mean(axis=-1)[:,np.newaxis]
		s4 = self.projO(s4,ds4.lon, ds4.lat)
		s8 = self.projO(s8,ds4.lon, ds4.lat)
		s16 = self.projO(s16,ds4.lon, ds4.lat)
		#s4   = self.projO(np.array([ds4.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0), ds4.lon, ds4.lat)
		#s8   = self.projO(np.array([ds8.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0), ds8.lon, ds8.lat)
		#s16  = self.projO(np.array([ds16.getSeason(Year=year,Season=self.Season).squeeze() for year in self.years]).mean(axis=0), ds16.lon, ds16.lat)

		cseq = 16#np.arange(-1.05,1.05+0.15,0.15)
		pl.figure(1)
		cf   = pl.contourf(self.projO.x,self.projO.y,s4,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cbar = pl.colorbar(cf)
		cbar.set_label('%s [%s]' % (ds4.long_name,ds4.units))
		self.projO.m.drawparallels([60,70,80],latmax=90)
		pl.title('CAM4xCO2')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.stationary.CAM4xCO2.%s-%s.%s.pdf' % (Field,self.years[0],self.years[-1],self.Season), format='pdf')
                pl.figure(2)
                cf   = pl.contourf(self.projO.x,self.projO.y,s8,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('%s [%s]' % (ds4.long_name,ds4.units))
                self.projO.m.drawparallels([60,70,80],latmax=90)
                pl.title('CAM8xCO2')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.stationary.CAM8xCO2.%s-%s.%s.pdf' % (Field,self.years[0],self.years[-1],self.Season), format='pdf')
                pl.figure(3)
                cf   = pl.contourf(self.projO.x,self.projO.y,s16,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('%s [%s]' % (ds4.long_name,ds4.units))
                self.projO.m.drawparallels([60,70,80],latmax=90)
                pl.title('CAM16xCO2')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.stationary.CAM16xCO2.%s-%s.%s.pdf' % (Field,self.years[0],self.years[-1],self.Season), format='pdf')
		pl.show()

	def GMST(self,Field='pw',case=''):
		# DataServers
                ds4  = MMDataServer(Field=Field,Source='CAM4xCO2%s' % (case))
                ds8  = MMDataServer(Field=Field,Source='CAM8xCO2%s' % (case))
                ds16 = MMDataServer(Field=Field,Source='CAM16xCO2%s' % (case))
		# Latitude weights
		w4  = np.cos( ds4.lat*np.pi/180.)
		w8  = np.cos( ds8.lat*np.pi/180.)
		w16 = np.cos(ds16.lat*np.pi/180.)
		# Data
		s4  = (np.array([ ds4.getSeason(Year=year,Season='Annual') for year in self.years]).mean(axis=0).mean(axis=-1)*w4).sum()/w4.sum()
		s8  = (np.array([ ds8.getSeason(Year=year,Season='Annual') for year in self.years]).mean(axis=0).mean(axis=-1)*w8).sum()/w8.sum()
		s16 = (np.array([ds16.getSeason(Year=year,Season='Annual') for year in self.years]).mean(axis=0).mean(axis=-1)*w16).sum()/w16.sum()

		print  '4xCO2 GMST = %s %s' % ( round(s4.mean(),0), ds4.units)
		print  '8xCO2 GMST = %s %s' % ( round(s8.mean(),0), ds8.units)
		print '16xCO2 GMST = %s %s' % (round(s16.mean(),0),ds16.units)

	def getEventIndexes(self,x,thresh,min_duration):
		# Given 1d array, finds all continuous segments > thresh and persisting for at least min_duration
		# returns indexes of the beginning and ends of each segment
		if x.max() < thresh: return []
		if (x >= thresh).all(): return [[0,len(x)-1]]
		# find indexes were x crosses threshold upward
		up = (np.nonzero( (x[:-1]<thresh) & (x[1:]>=thresh) )[0] + 1).tolist()
		if x[0]>=thresh: up = [0] + up
		# find indexes were x crosses threshold downward
		dn = np.nonzero( (x[:-1]>=thresh) & (x[1:]<thresh) )[0].tolist()
		if x[-1]>=thresh: dn = dn + [len(dn)-1]
		# if first crossing is downward, drop it
		if dn[0] < up[0]: dn.pop(0)
		if len(up) > len(dn): up.pop(-1)
		events = zip(up,dn)
		events = [e for e in events if (e[1]-e[0]+1) >= min_duration]
		return events

	def vert1d(self,Field='T',case=''):
		# Plots vertical climatologies for each experiment
		# DataServers
		ds4  = MMDataServer(Field=Field,Source='CAM4xCO2%s' % (case))
		ds8  = MMDataServer(Field=Field,Source='CAM8xCO2%s' % (case))
		ds16 = MMDataServer(Field=Field,Source='CAM16xCO2%s' % (case))
		# Vertical climatologies
		s4  = self.proj(np.array([ ds4.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0), ds4.lon, ds4.lat).mean(axis=-1).mean(axis=-1)
		s8  = self.proj(np.array([ ds8.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0), ds8.lon, ds8.lat).mean(axis=-1).mean(axis=-1)
		s16 = self.proj(np.array([ds16.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0),ds16.lon,ds16.lat).mean(axis=-1).mean(axis=-1)
		# Plot
		pl.plot( s4, ds4.lev,'b',linewidth=2,alpha=0.65,label='4xCO2%s' % (case))
		pl.plot( s8, ds8.lev,'g',linewidth=2,alpha=0.65,label='8xCO2%s' % (case))
		pl.plot(s16,ds16.lev,'r',linewidth=2,alpha=0.65,label='16xCO2%s' % (case))
		pl.xlabel('%s [%s]' % (ds4.long_name,ds4.units))
		pl.ylabel('Pressure [hPa]')
		pl.grid()
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.ylim(1000,0)
		pl.show()

	def interpolateFluxFile(self,vq):
	        # Interpolates file of shape (ntime,nlon) to (ntime,360)
	        n0,n1 = vq.shape
	        print 'Interpolating data to 1 degree resolution...'
	        # Add wrap value to vq
	        ends = vq[:,-1]
	        vq   = np.append(vq,ends[:,np.newaxis],1)
	        # Interpolate vq to 1 degree resolution
	        yold = range(n0)
	        ynew = range(n0)
	        xold = np.array(list(np.arange(0,360,360./n1))+[360])
	        xnew = np.arange(0,360,1.)
	        f    = interpolate.interp2d(xold,yold,vq,kind='linear',bounds_error=True)
	        vq   = f(xnew,ynew)
	        return vq,np.arange(360)

	def fluxPDF(self,blat=70,case='',fcase='moist',normed=True):
                if case == '': cstr = 'low'
                if case != '': cstr = 'high'
		sf,sfstr = 1,''
		if fcase == 'mass': sf,sfstr = 1000,'10$^{3}$ '
		# Flux data
		vq4 ,dates4  = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s%s/%s%s.%s.%s-%s.%s-%shPa.%sN.%s.p' % (blat, 'CAM4xCO2',case,'CAM4xCO2',case,fcase,1901,1919,0,1000,blat,self.Season))
		vq8 ,dates8  = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s%s/%s%s.%s.%s-%s.%s-%shPa.%sN.%s.p' % (blat, 'CAM8xCO2',case,'CAM8xCO2',case,fcase,1901,1919,0,1000,blat,self.Season))
		vq16,dates16 = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s%s/%s%s.%s.%s-%s.%s-%shPa.%sN.%s.p' % (blat,'CAM16xCO2',case,'CAM16xCO2',case,fcase,1901,1919,0,1000,blat,self.Season))
		# Reshape fluxes
		vq4,lons  =  self.interpolateFluxFile(vq4)
		vq8,lons  =  self.interpolateFluxFile(vq8)
		vq16,lons = self.interpolateFluxFile(vq16)
		# Histograms
                h4,e4   = np.histogram( vq4.reshape(-1)[np.where(vq4.reshape(-1)>0)[0]]/sf,45,normed=normed)
                h8,e8   = np.histogram( vq8.reshape(-1)[np.where(vq8.reshape(-1)>0)[0]]/sf,45,normed=normed)
                h16,e16 = np.histogram(vq16.reshape(-1)[np.where(vq16.reshape(-1)>0)[0]]/sf,45,normed=normed)
                edges4  = (e4[0:-1]+e4[1:])/2
                edges8  = (e8[0:-1]+e8[1:])/2
                edges16 = (e16[0:-1]+e16[1:])/2
		# Label values
		lv4  = round( stats.scoreatpercentile( vq4.reshape(-1)[np.where(vq4.reshape(-1)>0)[0]],85),2)
		lv8  = round( stats.scoreatpercentile( vq8.reshape(-1)[np.where(vq8.reshape(-1)>0)[0]],85),2)
		lv16 = round( stats.scoreatpercentile(vq16.reshape(-1)[np.where(vq16.reshape(-1)>0)[0]],85),2)
                #lv4  = round( vq4.reshape(-1).std()/sf,2)
                #lv8  = round( vq8.reshape(-1).std()/sf,2)
                #lv16 = round( vq16.reshape(-1).std()/sf,2)
                # Plot
                pl.plot(edges4 ,  h4, 'b'  , linewidth=2, alpha=0.65, label= '4xCO2%s; %s %s%s' % (case,lv4, sfstr,'Tg day$^{-1}$ deg$^{-1}$'))
                pl.plot(edges8 ,  h8, 'g'  , linewidth=2, alpha=0.65, label= '8xCO2%s; %s %s%s' % (case,lv8, sfstr,'Tg day$^{-1}$ deg$^{-1}$'))
                pl.plot(edges16, h16, 'r'  , linewidth=2, alpha=0.65, label='16xCO2%s; %s %s%s' % (case,lv16,sfstr, 'Tg day$^{-1}$ deg$^{-1}$'))
		pl.plot( edges4[-1], h4[-1],'bo')
		pl.plot( edges8[-1], h8[-1],'go')
		pl.plot(edges16[-1],h16[-1],'ro')
                pl.xlabel('%s %s [%s%s]' % (fcase,'flux',sfstr,'Tg day$^{-1}$ deg$^{-1}$'))
                pl.ylabel('Frequency')
                pl.grid()
                pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		#pl.xscale('log')
		#pl.yscale('log')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%sflux.%s.hist.%s-%s.%sN.%s.pdf' % (fcase,cstr,self.years[0],self.years[-1],blat,self.Season), format='pdf')
                pl.show()

	def pdf(self,Field='fls',case='',showRE=False,normed=True):
		if Field == 'Ts':   LevType,sf = 'surface_analysis',  1
		if Field == 'fls':  LevType,sf = 'surface_forecast', -1
		if Field == 'pw':   LevType,sf = 'surface_analysis',  1
		if Field == 'sshf': LevType,sf = 'surface_forecast', -1
		if Field == 'slhf': LevType,sf = 'surface_forecast', -1
		# Calculates and plots histograms for each experiment
		if case == '': cstr = 'low'
		if case != '': cstr = 'high'
                # DataServers
		if showRE: dsE  = reDataServer(Field=Field,Source='ERAInt',LevType=LevType,LatRange=(50,90))
                ds4  = reDataServer(Field=Field,Source= 'CAM4xCO2%s' % (case))
                ds8  = reDataServer(Field=Field,Source= 'CAM8xCO2%s' % (case))
                ds16 = reDataServer(Field=Field,Source='CAM16xCO2%s' % (case))
		# Get data
		if showRE:
			sE  = np.zeros((0,self.proj.nx,self.proj.ny))
			for year in self.yearsE:
				print year
				sE = np.append(sE,self.proj(dsE.getDataSnaps(Year=year,Season=self.Season),dsE.lon,dsE.lat),axis=0)	
			sE  = sE.reshape(-1) 
		s4  = sf*self.proj(np.array([ ds4.getDataSnaps(Year=year,Season=self.Season) for year in  self.years]),  ds4.lon,  ds4.lat).reshape(-1)
		s8  = sf*self.proj(np.array([ ds8.getDataSnaps(Year=year,Season=self.Season) for year in  self.years]),  ds8.lon,  ds8.lat).reshape(-1)
		s16 = sf*self.proj(np.array([ds16.getDataSnaps(Year=year,Season=self.Season) for year in  self.years]), ds16.lon, ds16.lat).reshape(-1)
		# Take only points above self.blat
		if showRE: latE = np.tile(self.proj.lat,(len(sE)/(self.proj.nx**2),1,1)).reshape(-1)
		lat4 = np.tile(self.proj.lat,( len(self.years),  len(s4)/(len(self.years)*(self.proj.nx**2)), 1, 1) ).reshape(-1)
		if showRE: sE   =  sE[np.where(latE>=self.blat)]
		s4   =  s4[np.where(lat4>=self.blat)]
		s8   =  s8[np.where(lat4>=self.blat)]
		s16  = s16[np.where(lat4>=self.blat)]
		# PDFs
		if showRE: hE,eE   = np.histogram( sE,45,normed=normed)
		h4,e4   = np.histogram( s4,45,normed=normed)
		h8,e8   = np.histogram( s8,45,normed=normed)
		h16,e16 = np.histogram(s16,45,normed=normed)
		if showRE: edgesE  = (eE[0:-1]+eE[1:])/2
		edges4  = (e4[0:-1]+e4[1:])/2
		edges8  = (e8[0:-1]+e8[1:])/2
		edges16 = (e16[0:-1]+e16[1:])/2
		# Plot
		if showRE: pl.plot(edges4 ,  hE, 'k--', linewidth=1, alpha=0.65, label=  'ERAInt: %s-%s' % (self.yearsE[0],self.yearsE[-1]))
		pl.plot(edges4 ,  h4, 'b'  , linewidth=2, alpha=0.65, label= '4xCO2%s; %s %s' % (case,round( s4.std(),2), ds4.units))
		pl.plot(edges8 ,  h8, 'g'  , linewidth=2, alpha=0.65, label= '8xCO2%s; %s %s' % (case,round( s8.std(),2), ds8.units))
		pl.plot(edges16, h16, 'r'  , linewidth=2, alpha=0.65, label='16xCO2%s; %s %s' % (case,round(s16.std(),2),ds16.units))
		pl.xlabel('%s [%s]' % (ds4.long_name,ds4.units))
		pl.ylabel('Frequency')
		pl.grid()
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.%s.hist.%s-%s.%sN.%s.pdf' % (Field,cstr,self.years[0],self.years[-1],self.blat,self.Season), format='pdf')
		pl.show()

	def comp2d(self,Source,Field,inField,thresh=0,min_duration=8):
		# Creates lagged compostie figures (full and anom) of a vertical field
		# centered on the time of maximum/minimum index field for each segment satisfying
		# thresh and min_duration; see self.getEventIndexes(). Also includes a composite of
		# SEB components and surface temperature.
		if inField == 'pw':  sf,inUnits =  1,'kgkg'
		if inField == 'fls': sf,inUnits = -1,'Wm2'
		days = np.arange(-4,4+0.25,0.25)
		# DataServers for index field and composite field
		dsCMM = MMDataServer(Field=Field  ,Source=Source)
		dsI   = reDataServer(Field=inField,Source=Source)
		dsC   = reDataServer(Field=Field  ,Source=Source)
		# Turbulent fluxes and surface temperature DataServers
		dsL   = reDataServer(Field='slhf',Source=Source)
		dsH   = reDataServer(Field='sshf',Source=Source)
		dsN   = reDataServer(Field='fls',Source=Source)
		dsT   = reDataServer(Field='Ts'  ,Source=Source)
		# Composite climatology
		clim        = self.proj( np.array([dsCMM.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0), dsC.lon, dsC.lat )
		Tf,Tc       = [],[]
		LH,SH,TS,NL = [],[],[],[]
		n           = 0
		for year in self.years:
			print year
			sI = sf*self.proj(dsI.getDataSnaps(Year=year,Season=self.Season),dsI.lon,dsI.lat)
			sL = -1*self.proj(dsL.getDataSnaps(Year=year,Season=self.Season),dsL.lon,dsL.lat)
			sH = -1*self.proj(dsH.getDataSnaps(Year=year,Season=self.Season),dsH.lon,dsH.lat)
			sT =    self.proj(dsT.getDataSnaps(Year=year,Season=self.Season),dsT.lon,dsT.lat)
			sN = -1*self.proj(dsN.getDataSnaps(Year=year,Season=self.Season),dsN.lon,dsN.lat)
			sC =    self.proj(dsC.getDataSnaps(Year=year,Season=self.Season),dsC.lon,dsC.lat)
			for i in range(self.proj.nx):
				for j in range(self.proj.ny):
					if self.proj.lat[i,j] >= self.blat:	
						events = self.getEventIndexes(sI[:,i,j],thresh,min_duration)
						for e0,e1 in events:	
							# Choose which portion of event to use as central lag
							# i.e. min,max,first last ... max: i0 = e0 + np.argmax(sI[e0:e1+1,i,j])
							i0 = e1
							if (i0>=16) and (i0<len(sI)-16):
								Tc.append(clim[:,i,j])
								Tf.append(sC[i0-16:i0+16+1,:,i,j])	
								LH.append(sL[i0-16:i0+16+1,i,j])
								SH.append(sH[i0-16:i0+16+1,i,j])
								TS.append(sT[i0-16:i0+16+1,i,j])
								NL.append(sN[i0-16:i0+16+1,i,j])
								n = n + 1
		# Composite events in vertical
		Tf = np.array(Tf).mean(axis=0)
		Tc = np.array(Tc).mean(axis=0)
		Ta = Tf - Tc[np.newaxis,:]
                # Roll axes for plotting
                Tf = np.rollaxis(Tf,1,0)
                Ta = np.rollaxis(Ta,1,0)
		# Composite events for surface
		LH = np.array(LH).mean(axis=0)
		SH = np.array(SH).mean(axis=0)
		TS = np.array(TS).mean(axis=0)
		NL = np.array(NL).mean(axis=0)
		# Plot
		cseqa = np.arange(-6,6+1,1)
		cseqf = 14
		pl.figure(1)
		cf   = pl.contourf(days,dsC.lev,Ta,cseqa,cmap=pl.cm.RdBu_r,extend='both')
		cbar = pl.colorbar(cf)
		cbar.set_label('%s [%s]' % (dsC.long_name,dsC.units))
		pl.xlabel('Days')
		pl.ylabel('Pressure [hPa]')
		pl.ylim(1000,0)
		pl.title(n)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.anom.%s.%s.%s%s.%sdt.%s-%s.%sN.%s.pdf' % (Field,Source,inField,thresh,inUnits,min_duration,self.years[0],self.years[-1],self.blat,self.Season), format='pdf')
                pl.figure(2)
                cf   = pl.contourf(days,dsC.lev,Tf,cseqf,cmap=pl.cm.OrRd,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('%s [%s]' % (dsC.long_name,dsC.units))
                pl.xlabel('Days')
                pl.ylabel('Pressure [hPa]')
                pl.ylim(1000,0)
                pl.title('%s: %s' % (Source,n))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.full.%s.%s.%s%s.%sdt.%s-%s.%sN.%s.pdf' % (Field,Source,inField,thresh,inUnits,min_duration,self.years[0],self.years[-1],self.blat,self.Season), format='pdf')

		fig,ax1 = pl.subplots(num=3)
		ax2     = ax1.twinx()
		ax1.plot(days,TS,'r--',linewidth=2.00,alpha=0.65,label='Ts')
		ax2.plot(days,LH,'b'  ,linewidth=1.25,alpha=0.65,label='lhf')
		ax2.plot(days,SH,'y'  ,linewidth=1.25,alpha=0.65,label='shf')
		ax2.plot(days,NL,'r'  ,linewidth=1.25,alpha=0.65,label='NetLW')
		ax1.set_ylabel('%s [%s]' % (dsT.long_name,dsT.units))
		ax2.set_ylabel('Surface energy flux [W/m2]')
		ax1.set_xlabel('Days')
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.title('%s: %s' % (Source,n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/SEB.full.%s.%s.%s%s.%sdt.%s-%s.%sN.%s.pdf' % (Source,inField,thresh,inUnits,min_duration,self.years[0],self.years[-1],self.blat,self.Season))

		pl.show()

	def comp_vert_expt(self,Field='T',Source='CAM4xCO2'):
		# Composite vertical profiles above and below arbitrary thresholds
		# of a gievn index field, for a given experiment
		# Index field
		dsI = reDataServer(Field='fls',Source=Source)
		sI  = -1*self.proj(np.array([ dsI.getDataSnaps(Year=year,Season=self.Season) for year in self.years]),  dsI.lon,  dsI.lat)
		sI  = sI.reshape(-1)
		# Comp field
		dsC = reDataServer(Field=Field,Source=Source)
		sC  = self.proj(np.array([dsC.getDataSnaps(Year=year,Season=self.Season) for year in self.years]).reshape(-1,dsC.nlev,dsC.nlat,dsC.nlon), dsC.lon, dsC.lat)
		sC  = np.rollaxis(sC,1,0).reshape((dsC.nlev,-1))
		sC0 = sC[:,np.where(sI<-30)[0]]
		sC1 = sC[:,np.where(sI>10)[0]]
		# Mean and standard deviation
		sC0m = sC0.mean(axis=1)
		sC0s = sC0.std(axis=1)
                sC1m = sC1.mean(axis=1)
                sC1s = sC1.std(axis=1)
		return sC0m,sC0s,sC1m,sC1s

	def comp_vert_all(self,Field='T'):
		# DataServer for attributes
		dsC = reDataServer(Field=Field,Source='CAM4xCO2')
		# Make composites
		s4C0m,s4C0s,s4C1m,s4C1s = self.comp_vert_expt(Field,Source='CAM4xCO2')	
		# Plot
		pl.plot(s4C0m , dsC.lev, 'b'  , linewidth=2.0,alpha=0.7)
		pl.plot(s4C1m , dsC.lev, 'b--', linewidth=1.2,alpha=0.7)
		pl.grid()
		pl.ylabel('Pressure [hPa]')
		pl.xlabel('%s [%s]' % (dsC.long_name,dsC.units))
		pl.ylim(1000,0)
		pl.show()

	def plotClims(self,Field='Ts',case=''):
		if case == '': cstr = 'low'
                if case != '': cstr = 'high'
		# DataServers
		ds4  = MMDataServer(Field=Field,Source= 'CAM4xCO2%s' % (case),LatRange=(0,90))
		ds8  = MMDataServer(Field=Field,Source= 'CAM8xCO2%s' % (case),LatRange=(0,90))
		ds16 = MMDataServer(Field=Field,Source='CAM16xCO2%s' % (case),LatRange=(0,90))
		# Clims
		s4  = np.array([ ds4.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0).mean(axis=-1)
		s8  = np.array([ ds8.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0).mean(axis=-1)
		s16 = np.array([ds16.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0).mean(axis=-1)
		#s4  =  s4 -  s4[-1]
		#s8  =  s8 -  s8[-1]
		#s16 = s16 - s16[-1]
		# Plot
		pl.plot( s4, ds4.lat,'b',linewidth=2,alpha=0.7,label='CAM4xCO2%s'  % (case))
		pl.plot( s8, ds8.lat,'g',linewidth=2,alpha=0.7,label='CAM8xCO2%s'  % (case))
		pl.plot(s16,ds16.lat,'r',linewidth=2,alpha=0.7,label='CAM16xCO2%s' % (case))
		pl.grid()
		pl.xlabel('Zonal-mean %s [%s]' % (ds4.long_name,ds4.units))
		pl.ylabel('Latitude')
		pl.ylim(ds4.lat[0],ds4.lat[-1])
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.zonalmean.%s.%s-%s.%s.pdf' % (Field,cstr,self.years[0],self.years[-1],self.Season),format='pdf')
		pl.show()

	def butter_bandpass(self, lowcut, highcut, fs, order=5):
        	nyq = 0.5 * fs
        	low = lowcut / nyq
        	high = highcut / nyq
        	b, a = butter(order, [low, high], btype='band')
        	return b, a

	def butter_bandpass_filter(self, data, lowcut, highcut, fs, order=5):
	        b, a = self.butter_bandpass(lowcut, highcut, fs, order=order)
	        y = lfilter(b, a, data, axis=0)
		return y

	def bandpassFilter(self,Source,lwl,hwl):
		# Band pass between lwl and hwl (in days)
	        # Make file
	        fname = '/mnt/climstorage/cian/scripts/synopfiles/%s.msl.%s-%sdays-pass.std.%s-%s.%s.p' % (Source,lwl,hwl,self.years[0],self.years[-1],self.Season)
	        print fname
	        if os.path.isfile(fname) == False:
	                # MSL DataServer
	                ds = reDataServer(Field='slp',Source=Source)
	                S1 = []
	                for year in self.years:
	                        x1 = ds.getDataSnaps(year,Season=self.Season)
	                        print year,x1.shape
	                        x1 = x1-x1.mean(axis=0)
	                        if (lwl!=0) and (hwl!=0):
	                                # Sample rate and desired cutoff frequencies (days per season).
	                                fs1 = len(x1)
	                                lowcut,highcut = 1.*fs1/(4*hwl),1.*fs1/(4*lwl)
	                                # Apply filter
	                                y1 = self.butter_bandpass_filter(x1, lowcut, highcut, fs1, order=6)
	                        else:
	                                y1 = x1
	                        S1.append(y1.std(axis=0))
	                S1 = np.array(S1)
	                lons1,lats1 = ds.lon,ds.lat
	                toPick([S1,lons1,lats1],fname)
	        else:
	                S1,lons1,lats1 = unpick(fname)
		return S1,lons1,lats1

	def fieldBandVar(self,lwl,hwl,case=''):
		if case == '': cstr =  'low'
		if case != '': cstr = 'high'
		# Get variation data
		s4 , lons4 , lats4  = self.bandpassFilter( 'CAM4xCO2%s' % (case),lwl,hwl)
		s8 , lons8 , lats8  = self.bandpassFilter( 'CAM8xCO2%s' % (case),lwl,hwl)
		s16, lons16, lats16 = self.bandpassFilter('CAM16xCO2%s' % (case),lwl,hwl)
		# Project to Lambert
		s4  = self.projO( s4.mean(axis=0), lons4, lats4)/100.
		s8  = self.projO( s8.mean(axis=0), lons8, lats8)/100.
		s16 = self.projO(s16.mean(axis=0),lons16,lats16)/100.
		# Plot
		cseqa = np.arange(-1.6,1.6+0.2,0.2)
		cseqf = np.arange(0,8+0.5,0.5)
		pl.figure(1)
		cf   = pl.contourf(self.projO.x,self.projO.y,s8-s4,cseqa,cmap=pl.cm.RdBu_r,extend='both')
		cbar = pl.colorbar(cf)
		cbar.set_label('Standard deviation of \nSea Level Pressure [hPa]')
		self.projO.m.drawparallels([60,70,80],latmax=90)
		pl.title('8xCO2 - 4xCO2: %s' % (cstr))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/synop.%s-%sdays.%s.8-4xCO2.%s-%s.%s.pdf' % (lwl,hwl,cstr,self.years[0],self.years[-1],self.Season))
                pl.figure(2)
                cf   = pl.contourf(self.projO.x,self.projO.y,s16-s8,cseqa,cmap=pl.cm.RdBu_r,extend='both')
                cbar = pl.colorbar()
		cbar.set_label('Standard deviation of \nSea Level Pressure [hPa]')
                self.projO.m.drawparallels([60,70,80],latmax=90)
		pl.title('16xCO2 - 8xCO2: %s' % (cstr))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/synop.%s-%sdays.%s.16-8xCO2.%s-%s.%s.pdf' % (lwl,hwl,cstr,self.years[0],self.years[-1],self.Season))
                pl.figure(3)
                cf   = pl.contourf(self.projO.x,self.projO.y,s4,cseqf,cmap=pl.cm.OrRd,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('Standard deviation of \nSea Level Pressure [hPa]')
                self.projO.m.drawparallels([60,70,80],latmax=90)
                pl.title('4xCO2: %s' % (cstr))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/synop.%s-%sdays.%s.4xCO2.%s-%s.%s.pdf' % (lwl,hwl,cstr,self.years[0],self.years[-1],self.Season))
                pl.figure(4)
                cf   = pl.contourf(self.projO.x,self.projO.y,s8,cseqf,cmap=pl.cm.OrRd,extend='both')
                cbar = pl.colorbar()
                cbar.set_label('Standard deviation of \nSea Level Pressure [hPa]')
                self.projO.m.drawparallels([60,70,80],latmax=90)
                pl.title('8xCO2: %s' % (cstr))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/synop.%s-%sdays.%s.8xCO2.%s-%s.%s.pdf' % (lwl,hwl,cstr,self.years[0],self.years[-1],self.Season))
                pl.figure(5)
                cf   = pl.contourf(self.projO.x,self.projO.y,s16,cseqf,cmap=pl.cm.OrRd,extend='both')
                cbar = pl.colorbar()
                cbar.set_label('Standard deviation of \nSea Level Pressure [hPa]')
                self.projO.m.drawparallels([60,70,80],latmax=90)
                pl.title('16xCO2: %s' % (cstr))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/synop.%s-%sdays.%s.16xCO2.%s-%s.%s.pdf' % (lwl,hwl,cstr,self.years[0],self.years[-1],self.Season))
		pl.show()

	def binVert(self,Field,inField,bins,Source,cabnd0,cabnd1,adbnd,cfbnd0,cfbnd1,fdbnd):
		# Average vertical structure for gripoints in each bin of the index field
		sf = 1
		if inField == 'fls': sf = -1
		# DataServers
		dsCMM = MMDataServer(Field=  Field,Source=Source)
		dsI   = reDataServer(Field=inField,Source=Source)
		dsC   = reDataServer(Field=  Field,Source=Source)
		# Lats for masking
		lat4  = np.tile(self.proj.lat,(len(dsC.getDateList(Year=1900,Season=self.Season)), 1, 1) ).reshape(-1)
		latxs = np.where(lat4>=self.blat)[0]
                # Composite climatology
                clim  = self.proj( np.array([dsCMM.getSeason(Year=year,Season=self.Season) for year in self.years]).mean(axis=0), dsC.lon, dsC.lat )
		clim  = np.tile(clim,(len(dsC.getDateList(Year=1900,Season=self.Season)),1,1,1))
		clim  = np.rollaxis(clim,1,0).reshape((dsCMM.nlev,-1))
		# Central bins
		binsc = (bins[0:-1] + bins[1:])/2.
		Fa,Ff = np.zeros((len(bins)-1,dsC.nlev)),np.zeros((len(bins)-1,dsC.nlev))
		n     = np.zeros(len(bins)-1)
		x0s   = []
		for year in self.years:
			print year
			sI    = sf*self.proj(dsI.getDataSnaps(Year=year,Season=self.Season),dsI.lon,dsI.lat).reshape(-1)[latxs]
			sC    = self.proj(dsC.getDataSnaps(Year=year,Season=self.Season),dsC.lon,dsC.lat)
			sC    = np.rollaxis(sC,1,0).reshape((dsC.nlev,-1))[:,latxs]
			x0s.append(stats.scoreatpercentile(sI,50))
			for i in range(len(bins)-1):	
				xs      = np.where((bins[i]<=sI)&(sI<bins[i+1]))[0]
				Fa[i,:] = Fa[i,:] + (sC[:,xs] - clim[:,xs]).sum(axis=1)
				Ff[i,:] = Ff[i,:] + (sC[:,xs]             ).sum(axis=1)
				n[i]    = n[i]   + len(xs)
		x0 = round(np.mean(x0s),1)
		x1 = round(np.mean(x0s),1)
		nt = np.tile(n,(dsC.nlev,1))
		nt = np.rollaxis(nt,1,0)
		Fa = np.ma.masked_where(nt==0,Fa)/np.ma.masked_where(nt==0,nt)
		Ff = np.ma.masked_where(nt==0,Ff)/np.ma.masked_where(nt==0,nt)
		Fa = np.rollaxis(Fa,1,0)
		Ff = np.rollaxis(Ff,1,0)
		# Plot
		cseqa = np.arange(cabnd0,cabnd1+adbnd,adbnd)
		cseqf = np.arange(cfbnd0,cfbnd1+fdbnd,fdbnd)
		pl.figure(1)
		cf    = pl.contourf(binsc,dsC.lev,Fa,cseqa,cmap=pl.cm.RdBu_r,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('%s [%s]' % (dsC.long_name,dsC.units))
		pl.plot([x0,x1],[1000,0],'k--',linewidth=1.15,alpha=0.6)
		pl.xlabel('%s [%s]' % (dsI.long_name,dsI.units))
		pl.ylabel('Pressure [hPa]')
		pl.ylim(1000,0)
		pl.title(Source)
		pl.xscale('log')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.anom.indexbins.%s.%s.%s-%s.%s.%sN.pdf' % (Field,inField,Source,self.years[0],self.years[-1],self.Season,self.blat),format='pdf')
                pl.figure(2)
                cf    = pl.contourf(binsc,dsC.lev,Ff,cseqf,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('%s [%s]' % (dsC.long_name,dsC.units))
                pl.xlabel('%s [%s]' % (dsI.long_name,dsI.units))
                pl.ylabel('Pressure [hPa]')
                pl.ylim(1000,0)
                pl.title(Source)
		pl.xscale('log')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.full.indexbins.%s.%s.%s-%s.%s.%sN.pdf' % (Field,inField,Source,self.years[0],self.years[-1],self.Season,self.blat),format='pdf')
		pl.figure(3)
		nx0,nx1 = np.argmin((binsc-x0)**2),np.argmin((binsc-x1)**2)
		Ffn0    = (Ff[:,0:nx0]*n[np.newaxis,0:nx0]).sum(axis=1)/n[0:nx0].sum()
		Ffn1    = (Ff[:,nx1:]*n[np.newaxis,nx1:]).sum(axis=1)/n[nx1:].sum()
		pl.plot(Ffn0,dsC.lev,'b',linewidth=2,alpha=0.7,label='%s < %s %s' % (inField,x0,dsI.units))
		pl.plot(Ffn1,dsC.lev,'r',linewidth=2,alpha=0.7,label='%s > %s %s' % (inField,x1,dsI.units))
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.xlabel('%s [%s]' % (dsC.long_name,dsC.units))
		pl.ylabel('Pressure [hPa]')
		pl.ylim(1000,0)
		pl.grid()
		pl.title(Source)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.states.indexbins.%s.%s.%s-%s.%s.%sN.pdf' % (Field,inField,Source,self.years[0],self.years[-1],self.Season,self.blat),format='pdf')
		pl.figure(4)
		pl.plot(binsc,1.*n,'k',linewidth=2,alpha=0.7)
		pl.grid()
		pl.xlabel('%s [%s]' % (dsI.long_name,dsI.units))
		pl.ylabel('Frequency')
		pl.yscale('log')
		pl.title(Source)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/N.freq.indexbins.%s.%s.%s-%s.%s.%sN.pdf' % (inField,Source,self.years[0],self.years[-1],self.Season,self.blat),format='pdf')
		pl.show()

if __name__ == "__main__":

	cam = CAM(	blat       =          60,
			res        =         300,
			Season     =     'NDJFM',
			YearRange  = (1901,1919),
			YearRangeE = (1980,2012)	)

	#cam.intrusionDensity(lat0=60)
	#cam.plotClim('Z')
	#cam.GMST()
	#cam.binVert('q','pw',np.arange(0,70+1,1),'CAM16xCO2','',-0.009,0.009,0.001)
	#cam.binVert('T','pw',np.arange(0,70+1,1),'CAM16xCO2_high_continent',-18,18,2,170,290,10)
	#cam.comp2d(Source='CAM16xCO2',Field='T',inField='pw',thresh=30,min_duration=8)
	#cam.comp_vert_all(Field='T')
	#cam.fluxPDF(blat=70,case='_high_continent',fcase='moist',normed=True)
	#cam.pdf(Field='Ts',case='_high_continent',showRE=False,normed=False)
	#cam.vert1d(case='_high_continent')
	#cam.plotClims(Field='pw',case='')

	#cam.fieldBandVar(2,6,'')
