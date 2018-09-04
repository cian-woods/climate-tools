from netCDF4 import Dataset
from LambertProjector import *
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from ReanalysisDataServer import DataServer as reDataServer
from TrackDataServer import DataServer as trDataServer
from scipy import stats
from drawCoastlinesNoRivers import *
from UnPickle import *
from matplotlib.colors import LogNorm

import matplotlib.ticker as ticker
import numpy as np
import matplotlib.pyplot as pl

class BarentsSeaBudget:

	def __init__(	self,
			blat     =      70,
			res      =     100.,
			LatRange =  (70,80),
			LonRange =  (15,60),
			LatRange0 = (65,75),
			LonRange0 = (60,100),
                        LatRange1 = (70,80),
                        LonRange1 = (340,15)	):

		# Attributes and LambertProjector
		self.blat         = blat
		self.res          = res
		self.LonRange     = LonRange
		self.LatRange     = LatRange
                self.LonRange0    = LonRange0
                self.LatRange0    = LatRange0
                self.LonRange1    = LonRange1
                self.LatRange1    = LatRange1
		self.Re           = 6370e03
		self.proj         = LambertProjector(boundinglat=self.blat,resolution=self.res)
		self.x,self.y     = self.proj.x[0,:],self.proj.y[:,0]
		self.gridA        = (self.proj.res*1000)**2	# m^2
		self.x70,self.y70 = self.proj.m(np.linspace(0,360,1000), np.zeros(1000)+70 )
		self.x80,self.y80 = self.proj.m(np.linspace(0,360,1000), np.zeros(1000)+80 )
		# Cropping info
		dim = 340e04
                self.x00,self.y00   = self.proj.m(340,62)
                self.x11,self.y11   = self.x00+dim,self.y00+dim
                self.ix00,self.ix11 = np.argmin((self.x - self.x00)**2),np.argmin((self.x - self.x11)**2)
                self.iy00,self.iy11 = np.argmin((self.y - self.y00)**2),np.argmin((self.y - self.y11)**2)
		if self.ix11 - self.ix00 > self.iy11 - self.iy00:
			self.ix11 = self.ix00 + (self.iy11 - self.iy00)
		elif self.ix11 - self.ix00 < self.iy11 - self.iy00:
			self.iy11 = self.iy00 + (self.ix11 - self.ix00)
		self.projx,self.projy = self.x[self.ix00:self.ix11+1],self.y[self.iy00:self.iy11+1]
		# Mask
		plon360    = np.zeros((self.proj.nx,self.proj.ny))
		plon360[np.where(self.proj.lon<0)] = 360
		self.plonf = self.proj.lon + plon360
		if LonRange[1]>LonRange[0]:
			self.mask = (self.proj.lat<LatRange[1])&(self.proj.lat>LatRange[0])&(self.plonf<LonRange[1])&(self.plonf>LonRange[0])==False
		else:	
			lonmask   = (self.plonf<LonRange[1])|(self.plonf>LonRange[0])
			self.mask = (self.proj.lat<LatRange[1])&(self.proj.lat>LatRange[0])&(lonmask)==False
                if LonRange0[1]>LonRange0[0]:
                        self.mask0 = (self.proj.lat<LatRange0[1])&(self.proj.lat>LatRange0[0])&(self.plonf<LonRange0[1])&(self.plonf>LonRange0[0])==False
                else:
                        lonmask0   = (self.plonf<LonRange0[1])|(self.plonf>LonRange0[0])
                        self.mask0 = (self.proj.lat<LatRange0[1])&(self.proj.lat>LatRange0[0])&(lonmask0)==False
                if LonRange1[1]>LonRange1[0]:
                        self.mask1 = (self.proj.lat<LatRange1[1])&(self.proj.lat>LatRange1[0])&(self.plonf<LonRange1[1])&(self.plonf>LonRange1[0])==False
                else:
                        lonmask1   = (self.plonf<LonRange1[1])|(self.plonf>LonRange1[0])
                        self.mask1 = (self.proj.lat<LatRange1[1])&(self.proj.lat>LatRange1[0])&(lonmask1)==False
		# DataServers
		self.trds = trDataServer(Source='ERAInt',Type='fwrd',blat=70,steps=25)
		self.reci = reDataServer(Field =   'ci',LevType='surface_analysis')
		self.spds = MMDataServer(Field =  'slp' )
		self.cids = MMDataServer(Field =   'ci'	)
		self.lwds = MMDataServer(Field =  'fls'	)
		self.swds = MMDataServer(Field = 'fsns'	)
		self.ldds = MMDataServer(Field = 'flds' )
		self.shds = MMDataServer(Field = 'sshf'	)
		self.lhds = MMDataServer(Field = 'slhf'	)
		self.stds = MMDataServer(Field =  'sst'	)
		self.tsds = MMDataServer(Field =  'T2'  )
		self.vsds = MMDataServer(Field='V',LevRange=(1000,1000))
                self.usds = MMDataServer(Field='U',LevRange=(1000,1000))
                self.ztds = MMDataServer(Field='Z',LevRange=  (500,500))
		self.vqds = MMDataServer(Field='vq')
		self.uqds = MMDataServer(Field='uq')
		# BoxBounds
		self.bounds  = self.boxBounds(self.LonRange,self.LatRange,34)
		self.bounds0 = self.boxBounds(self.LonRange0,self.LatRange0,34)
		self.bounds1 = self.boxBounds(self.LonRange1,self.LatRange1,34)
		# LandSea mask
		self.lsm = self.proj(self.cids.lsm,self.cids.lsm_lon,self.cids.lsm_lat)
		self.lsm[np.where(self.lsm>0)]  = 1
		self.lsm[np.where(self.lsm==0)] = 0
		self.areaocean = np.ma.masked_where(self.lsm==1,np.ones((self.proj.nx,self.proj.ny))*self.gridA)
		self.areaocean = np.ma.masked_array(self.areaocean,self.mask).sum()
                self.arealand  = np.ma.masked_where(self.lsm==0,np.ones((self.proj.nx,self.proj.ny))*self.gridA)
                self.arealand  = np.ma.masked_array(self.arealand,self.mask).sum()

		# Bathymetry data
		bd = Dataset('/mnt/climstorage/cian/Bathymetry_etopo2-cian.nc','r')
		bath_lon = bd.variables['lon'][:]
		bath_lat = bd.variables['lat'][:]
		self.bathy    = self.proj(-1*bd.variables['Bathymetry'][:].squeeze(),bath_lon,bath_lat)
		b_ = np.ma.masked_where(self.lsm==1,self.bathy)
		b_ = np.ma.masked_array(b_,mask=self.mask)

	def n_point_smooth(self,vf,n=1,iterations=1):
	        c         = 0
	        n0,n1     = vf.shape
	        vf_smooth = np.zeros((n0,n1)) 
		kernel    = np.array([[1,2,1],[2,3,2],[1,2,1]])
	        while c < iterations:
	                for i in range(n,n0-n):
	                        # Edges
	                        tmp00,tmp01     = vf[i-n:i+n+1,0],vf[i-n:i+n+1,-1]
	                        vf_smooth[i,0]  = (tmp00*kernel[0]).sum()/kernel[0].sum()
	                        vf_smooth[i,-1] = (tmp01*kernel[0]).sum()/kernel[0].sum()
	                        for j in range(n,n1-n):
	                                # Edges
	                                tmp10,tmp11     = vf[0,j-n:j+n+1],vf[-1,j-n:j+n+1]
	                                vf_smooth[0,j]  = (tmp10*kernel[0]).sum()/kernel[0].sum()
	                                vf_smooth[-1,j] = (tmp11*kernel[0]).sum()/kernel[0].sum()
	                                # Interior
	                                tmp = vf[i-n:i+n+1,j-n:j+n+1]
	                                vf_smooth[i,j] = (kernel*tmp).sum()/kernel.sum()
	                # Corners
	                vf_smooth[0,0]   = (vf[0,0] + vf[0,1] + vf[1,0]).mean()
	                vf_smooth[-1,-1] = (vf[-1,-1] + vf[-1,-2] + vf[-2,-1]).mean()
	                vf_smooth[0,-1]  = (vf[0,-1] + vf[0,-2] + vf[1,-1]).mean()
	                vf_smooth[-1,0]  = (vf[-1,0] + vf[-2,0] + vf[-1,1]).mean()
	                vf = vf_smooth
	                c  = c + 1
	        return vf

	def regress(self,YearRange,Season):

		years,hc,hc_detrend = self.OHC(YearRange,Season,plot=False)
		sst        = self.proj(np.array([self.stds.getSeason(Year=year,Season=Season) for year in years]),self.stds.lon,self.stds.lat)
		ci         = self.proj(np.array([self.cids.getSeason(Year=year,Season=Season) for year in years]),self.cids.lon,self.cids.lat)
		mask,mask_ = np.zeros((self.proj.nx,self.proj.ny)),ci>1.0
		for i in range(self.proj.nx):
                        for j in range(self.proj.ny):
				if (mask_[:,i,j].all()) or (self.lsm[i,j]==1) or (self.mask[i,j]==True):
					mask[i,j] = True
				else:
					mask[i,j] = False

		st      = np.ma.masked_array(sst,mask=np.tile(mask,(len(years),1,1))).mean(axis=1).mean(axis=1)
		sst,m,p = self.detrend2d(years,sst)
		r       = np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
				r[i,j] = np.corrcoef(hc_detrend,sst[:,i,j])[0][1]

		pl.figure(1)
		sx,sy = np.ma.masked_array(self.proj.x,mask=mask),np.ma.masked_array(self.proj.y,mask=mask)
		cseq  = np.arange(-0.6,0.6+0.1,0.1)
		cf    = pl.contourf(self.proj.x,self.proj.y,r,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Corr. coeff.')
		pl.plot(sx,sy,'k.',alpha=0.5)
		self.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)

                pl.figure(2)
		sx,sy = np.ma.masked_where(p>0.1,self.proj.x),np.ma.masked_where(p>0.1,self.proj.y)
                cseq  = 13#np.arange(-0.6,0.6+0.1,0.1)
                cf    = pl.contourf(self.proj.x,self.proj.y,m*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('%s [%s decade$^{-1}$]' % (self.stds.long_name,self.stds.units))
		pl.plot(sx,sy,'k.',alpha=0.5)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)

		pl.figure(3)
		pl.plot(years,st,'k.-',linewidth=1.5,alpha=0.7)
		pl.grid()
		pl.xlabel('Year')
		pl.ylabel('SST')
		pl.xlim(years[0],years[-1])

		pl.show()

	def OHC(self,YearRange,Season='Annual',plot=True):
		# Open file
		years = range(1955,2017+1,1)
		if Season == 'Annual':
			d = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_yearly.nc','r')
			HC   = d.variables['h18_hc'][:].squeeze()
		if Season == 'JFM':
			d  = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
			HC   = d.variables['h18_hc'][:].squeeze()[0::4]
                if Season == 'AMJ':
                        d  = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
                        HC   = d.variables['h18_hc'][:].squeeze()[1::4]
                if Season == 'JAS':
                        d  = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
                        HC   = d.variables['h18_hc'][:].squeeze()[2::4]
                if Season == 'OND':
                        d  = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
                        HC   = d.variables['h18_hc'][:].squeeze()[3::4]
		if Season == 'ONDJFM':
                        d   = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
                        HC0 = d.variables['h18_hc'][:].squeeze()[3::4][:-1,:,:]
			HC1 = d.variables['h18_hc'][:].squeeze()[4::4]
			HC  = (HC0+HC1)/2.
			years = range(1956,2017+1,1)
		if Season == '4xdecomp':
			d  = Dataset('/mnt/climstorage/cian/heat_content_anomaly_0-700_seasonal.nc','r')
			HC    = d.variables['h18_hc'][:].squeeze()[:]
			years = np.arange(1955,2017.75+0.25,0.25)
		try:
			xs0,xs1 = years.index(YearRange[0]),years.index(YearRange[1])
		except:
			xs0,xs1 = np.argmin((years-YearRange[0])**2),np.argmin((years-YearRange[1])**2)
		HC,years = HC[xs0:xs1+1],years[xs0:xs1+1]

		# Get data
		lats = d.variables['lat'][:]
		lons = d.variables['lon'][:]
		# Cut axes
		latx0,latx1 = np.argmin((lats-self.LatRange[0])**2),np.argmin((lats-self.LatRange[1])**2)
		lonx0,lonx1 = np.argmin((lons-self.LonRange[0])**2),np.argmin((lons-self.LonRange[1])**2)
		lats,lons   = lats[latx0:latx1+1],lons[lonx0:lonx1+1]
                llons,llats = np.meshgrid(lons,lats)
                weights     = np.cos(np.pi*lats/180)
                weights     = np.ma.masked_array(weights,mask=(lats<self.LatRange[0])|(lats>self.LatRange[1]))
		HC          = HC[:,latx0:latx1+1,lonx0:lonx1+1]

                # Mask
                plon360 = np.zeros((len(lats),len(lons)))
                plon360[np.where(llons<0)] = 360
                plonf = llons + plon360
                if self.LonRange[1]>self.LonRange[0]:
                        mask  = (llats<=self.LatRange[1])&(llats>=self.LatRange[0])&(plonf<=self.LonRange[1])&(plonf>=self.LonRange[0])==False
		HC = np.ma.masked_array(HC,mask=np.tile(mask,(len(HC),1,1)))
		HC_dt,trend,p = self.detrend2d(years,HC)
		trend = np.ma.masked_array(trend,mask=HC[0].mask)
		hc    = HC.sum(axis=1).sum(axis=1)
		hc    = hc - hc.mean()	# Make average heat content anomaly over period = 0
		slope,intercept,r_value,p_value,std_err = stats.linregress(years,hc)
		line = np.array([slope*x + intercept for x in years])

		if plot == True:
			# Plot
			pl.figure(1)
			pl.plot(years,hc,'k.-',linewidth=1.25,alpha=0.7)
			pl.plot(years,line,'k--',linewidth=1,alpha=0.6)
			pl.grid()
			pl.xlabel('Year')
			pl.ylabel('Ocean heat content anomaly [10$^{18}$ J]')
			pl.xlim(years[0],years[-1])
			pl.title('trend = %s 10$^{18}$ J year$^{-1}$' % (round(slope,3)))
			pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/OHC/ohc.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')
	
			#bmap = Basemap(projection='cyl',llcrnrlat=self.LatRange[0],urcrnrlat=self.LatRange[1],llcrnrlon=self.LonRange[0],urcrnrlon=self.LonRange[1],resolution='i')
			bmap = Basemap(projection='cyl',llcrnrlat=68.5,urcrnrlat=83.5,llcrnrlon=0.5,urcrnrlon=100.5,resolution='i')
			j,k  = bmap(llons,llats)
			fig  = pl.figure(2)
			ax   = fig.add_subplot(111)
			cseq = np.arange(-2,2+0.25,0.25)
			cf   = ax.contourf(j,k,trend*10,cseq,cmap=pl.cm.RdBu_r,extend='both')
			cbar = pl.colorbar(cf,orientation='horizontal')
			cbar.set_label('Ocean heat content trend [10$^{18}$ J decade$^{-1}$]')
			bmap.drawcoastlines(linewidth=0.5,color='0.5')
			bmap.fillcontinents(color='0.85')
			bmap.drawmeridians([0,20,40,60])
			bmap.drawparallels([70.5,75.5,80.5])
			ax.set_aspect('auto')
			#pl.xlim(self.LonRange[0],self.LonRange[1])
			#pl.ylim(self.LatRange[0],self.LatRange[1])
			pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/OHC/ohc.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

			pl.show()
		else:
			return years,hc,hc-line

        def boxBounds(self,LonRange,LatRange,n=35):
                if LonRange[1]>LonRange[0]:
			lons = np.linspace(LonRange[0],LonRange[1],n)
		else:
			lons = np.linspace(LonRange[0]-360,LonRange[1],n)
                lats = np.linspace(LatRange[0],LatRange[1],n)
                mer1 = self.proj.m([LonRange[0] for i in range(n)],lats)
                mer2 = self.proj.m([LonRange[1] for i in range(n)],lats)
                zon1 = self.proj.m(lons,[LatRange[0] for i in range(n)])
                zon2 = self.proj.m(lons,[LatRange[1] for i in range(n)])
                return [mer1,mer2,zon1,zon2]

	def areaLatLonBox(self,lat1,lat2,lon1,lon2):
		A = (np.pi/180)*(self.Re**2)*np.abs(np.sin(np.pi*lat1/180)-np.sin(np.pi*lat2/180))*np.abs(lon1-lon2)
		return A

	def detrend1d(self,years,data):
		m,c  = np.polyfit(years,data,1)	
		line = np.array([m*x + c for x in years])
		return data - line

        def detrend2d(self,years,field):
                n0,n1,n2 = field.shape
                m,c,p    = np.zeros((n1,n2)),np.zeros((n1,n2)),np.zeros((n1,n2))
                for i in range(n1):
                        for j in range(n2):
                                slope, intercept, r_value, p_value, std_err = stats.linregress(years,field[:,i,j])
                                m[i,j] = slope
                                c[i,j] = intercept
                                p[i,j] = p_value/2
                line  = np.array([m*year + c for year in years])
                field = field - line
                return field,m,p

	def surfacePDF(self,YearRange,Season,Freq='6hourly'):
		years = range(YearRange[0],YearRange[1]+1,1)
		# Sea ice
		if Freq == '6hourly':
			ci1 = self.proj(np.array([self.reci.getSeason(year,Season=Season1) for year in years]),self.reci.lon,self.reci.lat).reshape((-1,self.proj.nx,self.proj.ny))
		if Freq == 'Monthly':
			ci1 = self.proj(np.array([self.cids.getSeason(year,Season=Season1) for year in years]),self.cids.lon,self.cids.lat)
		tm0 = np.tile(self.mask,(len(ci1),1,1))
		tm1 = ((ci1<=1.00)&(tm0==False)&(np.tile(self.lsm,(len(ci1),1,1))==0))==False      # Surface fluxes
		# List for histogram
		H = []
		for i in range(len(ci1)):
			l   = ci1[i][np.where(tm1[i]==False)].reshape(-1)
			h,e = np.histogram(l,bins=np.linspace(0.01,1,50),normed=False)
			H.append(h)
		H   = np.array(H)
		Hm  = H.mean(axis=0)
		Hs  = H.std(axis=0)
		e   = (e[1:]+e[0:-1])/2

		pl.figure(1)
		ix0 = len(ci1)/2
		ix1 = len(years)/2
		pl.plot(e,Hm,'b',linewidth=1.6,alpha=0.5,label='%s-%s' % (years[0],years[-1]))
		pl.plot(e,Hm+Hs,'b--',linewidth=1.15,alpha=0.4)
		pl.plot(e,Hm-Hs,'b--',linewidth=1.15,alpha=0.4)
		pl.plot(e,H[0:ix0].mean(axis=0),'y',linewidth=1.15,alpha=0.4,label='%s-%s' % (years[0],years[ix1-1]))
		pl.plot(e,H[ix0:].mean(axis=0),'r',linewidth=1.15,alpha=0.4,label='%s-%s' % (years[ix1],years[-1]))
		pl.plot(e[0],Hm[0],'b.')
		pl.plot(e[-1],Hm[-1],'b.')
		pl.xlabel('%s [%s]' % (self.cids.long_name,self.cids.units))
		pl.ylabel('Frequency [%s$^{-1}$]' % (Season))
		pl.grid()
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.xlim(0,1)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/ci.hist.%s-%s.%s.%s.pdf' % (YearRange[0],YearRange[1],Season,Freq),format='pdf')
		pl.show()

	def Q_Aice(self,YearRange,Season1,Season2):
		years   = range(YearRange[0],YearRange[1]+1,1)
		# Intrusion data
		G,Q,D     = unpick('/mnt/climstorage/cian/newintrusions/ERAInt_intrusions.1980-2018.%s.6x6hr.9deg.200.6dt.20.5.70N.filtered.80N.p' % ('NDJF'))#(Season1))
		dates     = [D[i][0] for i in range(len(D))]
		LON,LAT,P = self.trds.getIntrusionTrajectories(G,D)
		nn,ns,years_ns,Ntot,Ttot,xx,yy = self.trds.density(LON,LAT,dates,self.proj,nend=None)
		ns = np.array([self.n_point_smooth(np.array(ns[i]).sum(axis=0),1,2) for i in range(len(ns)) if years_ns[i] in years])
		# BSO flux
		years_bso,BSO     = unpick('/mnt/climstorage/cian/BSO/hfTOT30d.bso.full.1998-2014.%s.p' % (Season1))
		years_bso,BSO_vt  = unpick('/mnt/climstorage/cian/BSO/fTOT30d.bso.full.1998-2014.%s.p'  % (Season1))
		BSO_T             = BSO/(1000*BSO_vt*1e06*4200)
		years_bso,BSO1    = unpick('/mnt/climstorage/cian/BSO/hfTOT30d.bso.full.1998-2014.Annual_wc.p')
		dBSO,dBSO_T,dBSO1 = self.detrend1d(years_bso,BSO),self.detrend1d(years_bso,BSO_T),self.detrend1d(years_bso,BSO1)
                # Mask and sea-ice series
                tm0  = np.tile(self.mask,(len(years),1,1))
		tm00 = np.tile(self.mask0,(len(years),1,1))
		tm01 = np.tile(self.mask1,(len(years),1,1))
		Season_con = 'NDJF'
		ci0  = self.proj(np.array([self.cids.getSeason(year,Season=Season_con) for year in years]),self.cids.lon,self.cids.lat)
                ci1  = self.proj(np.array([self.cids.getSeason(year,Season=Season1) for year in years]),self.cids.lon,self.cids.lat)
                ci2  = self.proj(np.array([self.cids.getSeason(year,Season=Season2) for year in years]),self.cids.lon,self.cids.lat)
		tm10 = ((ci0>=0.00)&( ci0<=0.15)&( tm0==False)&(np.tile(self.lsm,(len(years),1,1))==0))==False  # Surface fluxes (Season_con)
                tm1  = ((ci1>=0.00)&( ci1<=0.15)&( tm0==False)&(np.tile(self.lsm,(len(years),1,1))==0))==False  # Surface fluxes (Given Season1)
                tm2  = ((ci2<=1.00)&( ci2>=0.00)&( tm0==False)&(np.tile(self.lsm,(len(years),1,1))==0))==False  # Sea-ice area > 0%
		tm3  = ((ci1<=1.00)&( ci1>=0.15)&( tm0==False)&(np.tile(self.lsm,(len(years),1,1))==0))==False  # Sea-ice area > 15%
		tm1_ = ((np.tile(ci2.mean(axis=0),(len( years),1,1))>=0.00)&(np.tile(ci1.mean(axis=0),(len(years),1,1))<=0.00)&( tm0==False)&(np.tile(self.lsm,(len( years),1,1))==0))==False  # Surface fluxes
                # Fluxes and other fields
		tas = self.proj(np.array([self.tsds.getSeason(year,Season=Season1) for year in years]),self.tsds.lon,self.tsds.lat)
		sst = self.proj(np.array([self.stds.getSeason(year,Season=Season1) for year in years]),self.stds.lon,self.stds.lat)
		slp = self.proj(np.array([self.spds.getSeason(year,Season=Season1) for year in years]),self.spds.lon,self.spds.lat)/100.
		#ldf = self.proj(np.array([self.ldds.getSeason(year,Season=Season1) for year in years]),self.ldds.lon,self.ldds.lat)
                lwf = self.proj(np.array([self.lwds.getSeason(year,Season=Season1) for year in years]),self.lwds.lon,self.lwds.lat)
                swf = self.proj(np.array([self.swds.getSeason(year,Season=Season1) for year in years]),self.swds.lon,self.swds.lat)
                shf = self.proj(np.array([self.shds.getSeason(year,Season=Season1) for year in years]),self.shds.lon,self.shds.lat)
                lhf = self.proj(np.array([self.lhds.getSeason(year,Season=Season1) for year in years]),self.lhds.lon,self.lhds.lat)
		vs  = self.proj(np.array([self.vsds.getSeason(year,Season=Season1).squeeze() for year in years]),self.vsds.lon,self.vsds.lat)
		us  = self.proj(np.array([self.usds.getSeason(year,Season=Season1).squeeze() for year in years]),self.usds.lon,self.usds.lat)
		zt  = self.proj(np.array([self.ztds.getSeason(year,Season=Season1).squeeze() for year in years]),self.ztds.lon,self.ztds.lat)/9.81
		vq  = self.proj(np.array([self.vqds.getSeason(year,Season=Season1) for year in years]),self.vqds.lon,self.vqds.lat)
		uq  = self.proj(np.array([self.uqds.getSeason(year,Season=Season1) for year in years]),self.uqds.lon,self.uqds.lat)
                Thf = shf+lhf+lwf+swf
		Thf_clim = Thf.mean(axis=0)

		#for i in range(len(years)):
                #	pl.figure(1)
                #	cseq  = np.arange(-50,50+5,5)
                #	cf    = pl.contourf(self.proj.x,self.proj.y,Thf[i]-Thf_clim,cseq,cmap=pl.cm.coolwarm,extend='both')
                #	cbar  = pl.colorbar(cf)
                #	cbar.set_label('Total heat flux [W m$^{-2}$]')
                #	self.proj.m.drawcoastlines(color='0.85',linewidth=1,zorder=11)
                #	self.proj.m.fillcontinents(color='0.85',zorder=10)
                #	for x,y in self.bounds:
                #	        pl.plot(x,y,'g',linewidth=3,alpha=0.6,zorder=12)
                #	self.proj.m.drawparallels([70,80])
                #	drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
		#	pl.title(years[i])
		#	pl.show()

		dT  = tas - sst
		# Detrended fields
		#ldf_dt, m_ldf, p_ldf = self.detrend2d(years,ldf)
		Thf_dt, m_thf, p_thf = self.detrend2d(years,Thf)
		swf_dt, m_swf, p_swf = self.detrend2d(years,swf)
		lwf_dt, m_lwf, p_lwf = self.detrend2d(years,lwf)
		lhf_dt, m_lhf, p_lhf = self.detrend2d(years,lhf)
		shf_dt, m_shf, p_shf = self.detrend2d(years,shf)
		ci2_dt, m_ci2, p_ci2 = self.detrend2d(years,ci2)
		sst_dt, m_sst, p_sst = self.detrend2d(years,sst)
		slp_dt, m_slp, p_slp = self.detrend2d(years,slp)
		tas_dt, m_tas, p_tas = self.detrend2d(years,tas)
		ns_dt ,  m_ns,  p_ns = self.detrend2d(years,ns)
		dT_dt ,  m_dT,  p_dT = self.detrend2d(years,dT)
		vs_dt ,  m_vs,  p_vs = self.detrend2d(years,vs)
		us_dt ,  m_us,  p_us = self.detrend2d(years,us)
		zt_dt ,  m_zt,  p_zt = self.detrend2d(years,zt)
		vq_dt ,  m_vq,  p_vq = self.detrend2d(years,vq)
		uq_dt ,  m_uq,  p_uq = self.detrend2d(years,uq)

		# Downward longwave over ice
		#ldf_full = np.ma.masked_array(ldf,mask=tm3).mean(axis=1).mean(axis=1)
		#dldf     = np.ma.masked_array(ldf_dt,mask=tm3).mean(axis=1).mean(axis=1)
		# Intrusion density
                ns_full  = np.ma.masked_array(ns,mask=tm1).sum(axis=1).sum(axis=1)
                dns      = np.ma.masked_array(ns_dt,mask=tm1).mean(axis=1).mean(axis=1)
		# Sea-ice area
		ci2_tmp = ci2[:,:,:]
		ci2_tmp[np.where(ci2_tmp>=0.15)] = 1
		ci2_tmp[np.where(ci2_tmp<0.15)]  = 0
		ci_full  = (np.ma.masked_array(ci2_tmp,mask=tm2)*self.gridA).sum(axis=1).sum(axis=1)
		dci      = self.detrend1d(years,ci_full)
		#dci     = (np.ma.masked_array(ci2_dt,mask=tm2)*self.gridA).sum(axis=1).sum(axis=1)
                # Open ocean flux
                Thf_full0 = np.ma.masked_array(Thf,mask=tm10).mean(axis=1).mean(axis=1)
                Thf_full  = np.ma.masked_array(Thf,mask=tm1).mean(axis=1).mean(axis=1)		#*(self.areaocean-ci_full)
		dThf0     = np.ma.masked_array(Thf_dt,mask=tm10).mean(axis=1).mean(axis=1)	#*(self.areaocean-ci_full)
                dThf      = np.ma.masked_array(Thf_dt,mask=tm1).mean(axis=1).mean(axis=1)	#*(self.areaocean-ci_full)
                mThf      = np.polyfit(years,Thf_full,1)[0]
                mThf0     = np.polyfit(years,Thf_full0,1)[0]
		# SST
		sst_full = np.ma.masked_array(sst,mask=tm1).mean(axis=1).mean(axis=1)
		dsst     = np.ma.masked_array(sst_dt,mask=tm1).mean(axis=1).mean(axis=1)
                # SLP
                slp_full = np.ma.masked_array(slp,mask=tm00).mean(axis=1).mean(axis=1)
                dslp     = np.ma.masked_array(slp_dt,mask=tm00).mean(axis=1).mean(axis=1)
                # SLP (other region)
                slp_full1 = np.ma.masked_array(slp,mask=tm01).mean(axis=1).mean(axis=1)
                dslp1     = np.ma.masked_array(slp_dt,mask=tm01).mean(axis=1).mean(axis=1)
		# SLP difference
		slp_diff  = slp_full - slp_full1
		dslp_diff = self.detrend1d(years,slp_diff)
		# Z
                zt_full = np.ma.masked_array(zt,mask=tm00).mean(axis=1).mean(axis=1)
                dzt     = np.ma.masked_array(zt_dt,mask=tm00).mean(axis=1).mean(axis=1)
                # TAS
                tas_full = np.ma.masked_array(tas,mask=tm1).mean(axis=1).mean(axis=1)
		dtas     = np.ma.masked_array(tas_dt,mask=tm1).mean(axis=1).mean(axis=1)
		# Delta T
                dT_full  = np.ma.masked_array(dT,mask=tm1).mean(axis=1).mean(axis=1)
                ddT      = np.ma.masked_array(dT_dt,mask=tm1).mean(axis=1).mean(axis=1)
	
		# Regress series onto field(s)
		reg_thf_slp  , p_thf_slp  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_thf_vs   , p_thf_vs   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_thf_us   , p_thf_us   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_thf_ns   , p_thf_ns   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_dslp_thf , p_dslp_thf = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_thf_vq   , p_thf_vq   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		reg_thf_uq   , p_thf_uq   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
				m,c,r,p,er        = stats.linregress(dThf,slp_dt[:,i,j])
				reg_thf_slp[i,j]  = m
				p_thf_slp[i,j]    = p/2.
                                m,c,r,p,er        = stats.linregress(dThf,vs_dt[:,i,j])
                                reg_thf_vs[i,j]   = m
                                p_thf_vs[i,j]     = p/2.
                                m,c,r,p,er        = stats.linregress(dThf,us_dt[:,i,j])
                                reg_thf_us[i,j]   = m
                                p_thf_us[i,j]     = p/2.
                                m,c,r,p,er        = stats.linregress(dThf,ns_dt[:,i,j])
                                reg_thf_ns[i,j]   = m
                                p_thf_ns[i,j]     = p/2.
                                #m,c,r,p,er        = stats.linregress(dslp_diff,Thf_dt[:,i,j])
                                #reg_dslp_thf[i,j] = m
                                #p_dslp_thf[i,j]   = p/2.
                                m,c,r,p,er        = stats.linregress(dThf,vq_dt[:,i,j])
                                reg_thf_vq[i,j]   = m
                                p_thf_vq[i,j]     = p/2.
                                m,c,r,p,er        = stats.linregress(dThf,uq_dt[:,i,j])
                                reg_thf_uq[i,j]   = m
                                p_thf_uq[i,j]     = p/2.
		# Regress time series
		reg_dThf_dci,c,r,p,er = stats.linregress(dThf,dci)
		reg_dThf_dci0,c0,r0,p0,er0 = stats.linregress(dThf0,dci)

		trend_ci2 ,_c,_r,  ci_p,_er = stats.linregress(years,ci_full)
		trend_thf ,_c,_r, thf_p,_er = stats.linregress(years,Thf_full)
		trend_thf0,_c,_r,thf0_p,_er = stats.linregress(years,Thf_full0)

		#trend_ci2 = np.polyfit(years,ci_full,1)[0]
		#trend_thf = np.polyfit(years,Thf_full,1)[0]
		#trend_thf0 = np.polyfit(years,Thf_full0,1)[0]

		try:
                	print 'Corr. coeff HT_bso (Annual_wc) and ci = %s' % (np.corrcoef(dBSO1,dci)[0][1])
			#print 'Corr. coeff T_bso  and ci = %s' % (np.corrcoef(dBSO_T,dci)[0][1])
			print 'Corr. coeff HT_bso (Annual_wc) and HF = %s' % (np.corrcoef(dBSO1,dThf)[0][1])
			#print 'Corr. coeff T_bso and HF = %s' % (np.corrcoef(dBSO_T,dThf)[0][1])
		except:
			pass
		#print 'Corr. coeff. LWD and HF = %s' % (np.corrcoef(dldf,dThf)[0][1])
                #print 'Corr. coeff. SLP grad and HF = %s' % (np.corrcoef(dslp_diff,dThf)[0][1])

		print 'Corr. coeff. HF (%s mask) and ci = %s' % (Season_con,np.corrcoef(dci,dThf0)[0][1])
		print 'Reg.  coeff. HF (%s mask) and ci = %s 10**10 m**2 {W m**-2}**-1'   % (Season_con,reg_dThf_dci0/1e10)
		print 'p value of regression HF (%s mask) and ci = %s\n' % (Season_con,p0)

		print 'HF trend (%s mask) = %s W m**-2 year**-1' % (Season_con,trend_thf0)
		print 'HF climatology (%s mask) = %s W m**-2' % (Season_con,Thf_full0.mean())
		print 'p value of HF trend (%s mask) = %s\n' % (Season_con,thf0_p)

		print 'ci trend = %s 10**10 m**2 year**-1' % (trend_ci2/1e10)
		print 'ci climatology = %s 10**10 m**2' % (ci_full.mean()/1e10)
		print 'p value of ci trend = %s\n' % (ci_p)

		print 'Explained trend  = %s %s' % (100.*reg_dThf_dci*trend_thf/trend_ci2,'%')
		print 'Explained trend (%s mask) = %s %s' % (Season_con,100.*reg_dThf_dci0*trend_thf0/trend_ci2,'%')

		cc = 1
                fig,ax1  = pl.subplots(num=cc,figsize=(10,3))
		cc       = cc + 1
		Thf_anom = Thf_full - Thf_full.mean()
		m,c  = np.polyfit(years,Thf_anom,1)
		line = [m*x + c for x in years]
                ax1.plot(years,Thf_anom,'k-',linewidth=3.5,alpha=0.65)
		ax1.plot(years,line,'k--',linewidth=1,alpha=0.5)
		ax1.set_ylabel('Surface heat flux anomaly over open ocean [W m$^{-2}$]',fontsize=6)
		ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_ylim(-60,75)
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
		pl.title('trend = %s W m$^{-2}$ year$^{-1}$' % (mThf))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/hf.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')	

                fig,ax1   = pl.subplots(num=cc,figsize=(10,3))
                cc        = cc + 1
		Thf_anom0 = Thf_full0 - Thf_full0.mean()
                m,c       = np.polyfit(years,Thf_anom0,1)
                line0     = [m*x + c for x in years]
                ax1.plot(years,Thf_anom ,'k-',linewidth=3.5,alpha=0.65)
		ax1.plot(years,Thf_anom0,'b-',linewidth=3.0,alpha=0.65)
		ax1.plot(years,line ,'k--',linewidth=1,alpha=0.5)
		ax1.plot(years,line0,'b--',linewidth=1,alpha=0.5)
                ax1.set_ylabel('Ice free flux anomaly [W m$^{-2}$]')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
                ax1.set_ylim(-60,75)
                pl.title('trend = %s W m$^{-2}$ year$^{-1}$' % (mThf0))
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/hf_%smask.%s-%s.%s.pdf' % (Season_con,YearRange[0],YearRange[1],Season1),format='pdf')

                #fig,ax1 = pl.subplots(num=cc,figsize=(10,3))
                #cc      = cc + 1
                #ax1.plot(years,ldf_full,'k-',linewidth=1.5,alpha=0.65)
                #ax1.set_ylabel('Downward longwave over ice [W m$^{-2}$]')
                #ax1.set_xlabel('Year')
                #ax1.grid()
                #ax1.set_xlim(years[0],years[-1])
                #pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/lwd.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                fig,ax1 = pl.subplots(num=cc,figsize=(10,3))
                cc      = cc + 1
                ax1.plot(years,slp_full,'k-',linewidth=3.5,alpha=0.65)
                ax1.set_ylabel('Sea-level pressure [hPa]')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/slp.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                fig,ax1 = pl.subplots(num=cc,figsize=(10,3))
                cc      = cc + 1
                ax1.plot(years,slp_diff,'k-',linewidth=3.5,alpha=0.65)
                ax1.set_ylabel('Sea-level pressure difference [hPa]')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_ylim(20,-2.5)
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/dslp.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                fig,ax1  = pl.subplots(num=cc,figsize=(10,3))
		cc       = cc + 1
		ci_fulla = (ci_full - ci_full.mean())/1e10
		m,c      = np.polyfit(years,ci_fulla,1)
		line     = [m*x + c for x in years]
                ax1.plot(years,ci_fulla,'k-',linewidth=3.5,alpha=0.65)
		ax1.plot(years,line,'k--',linewidth=1,alpha=0.5)
                ax1.set_ylabel('March sea-ice area anomaly [10$^{10}$ m$^{2}$]')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_ylim(-40,40)
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/ci.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season2),format='pdf')

                fig,ax1    = pl.subplots(num=cc,figsize=(10,3))
		cc         = cc + 1
		BSOa,BSO1a = (BSO - BSO.mean())/1e12,(BSO1 - BSO1.mean())/1e12
                m,c        = np.polyfit(years_bso,BSOa,1)
                line       = [m*x + c for x in years_bso]
                m,c        = np.polyfit(years_bso,BSO1a,1)
                line1      = [m*x + c for x in years_bso]
                ax1.plot(years_bso,BSOa,'k-',linewidth=3.5,alpha=0.65)
		ax1.plot(years_bso,BSO1a,'b-',linewidth=3.0,alpha=0.65)
		ax1.plot(years_bso,line ,'k--',linewidth=1,alpha=0.5)
		ax1.plot(years_bso,line1,'b--',linewidth=1,alpha=0.5)
                ax1.set_ylabel('BSO heat transport anomaly [TW]')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_ylim(-25,40)
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/bso.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                fig,ax1 = pl.subplots(num=cc,figsize=(10,3))
		cc      = cc + 1
                ax1.plot(years,ns_full,'k-',linewidth=3.5,alpha=0.65)
                ax1.set_ylabel('Intrusion density')
                ax1.set_xlabel('Year')
                ax1.set_xlim(years[0],years[-1])
		ax1.set_ylim(400,0)
		ax1.set_xticks([i for i in range(1980,2018+1,1) if i not in [1980,1985,1990,1995,2000,2005,2010,2015]], minor=True)
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/N.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
		cc    = cc + 1
		cseq  = np.arange(-0.09,0.08+0.015,0.015)
		sx,sy = np.ma.masked_where(p_thf_slp>0.05,self.proj.x),np.ma.masked_where(p_thf_slp>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,reg_thf_slp,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('[hPa (W m$^{-2}$)$^{-1}$]')
                for x,y in self.bounds0:
                        pl.plot(x,y,'b',linewidth=3,zorder=10)
                for x,y in self.bounds1:
                        pl.plot(x,y,'r',linewidth=3,zorder=10)
		pl.plot(sx,sy,'k.',alpha=0.5)
		uv        = np.sqrt(reg_thf_us**2 + reg_thf_vs**2)
                urot,vrot = self.proj.m.rotate_vector(reg_thf_us,reg_thf_vs,self.proj.lon,self.proj.lat,returnxy=False)
                rotmask   = ((p_thf_us>0.05)&(p_thf_vs>0.05))|(uv<0.03)
                urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
                Q = pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],units='width',scale=1.25,width=0.003,headwidth=4.5,headlength=4,pivot='tail',alpha=0.85)
		qk = pl.quiverkey(Q, 0.4, 1.02, 0.05, '%s %s' % (0.05,'m s$^{-1}$ (W m$^{-2}$)$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/Thf-slp.reg.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                #pl.figure(cc)
                #cc    = cc + 1
                #cseq  = np.arange(-4,4+0.5,0.5)
                #sx,sy = np.ma.masked_where(p_dslp_thf>0.05,self.proj.x),np.ma.masked_where(p_dslp_thf>0.05,self.proj.y)
                #cf    = pl.contourf(self.proj.x,self.proj.y,reg_dslp_thf,cseq,cmap=pl.cm.coolwarm,extend='both')
                #cbar  = pl.colorbar(cf)
                #cbar.set_label('Reg coeff. [W m$^{-2}$ hPa$^{-1}$]')
                #for x,y in self.bounds0:
                #        pl.plot(x,y,'b',linewidth=3,zorder=10)
                #for x,y in self.bounds1:
                #        pl.plot(x,y,'r',linewidth=3,zorder=10)
                #pl.plot(sx,sy,'k.',alpha=0.5)
                #pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                #pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                #drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
                #pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/dslp-thf.reg.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
		cc    = cc + 1
                cseq  = np.arange(-1.2,1.2+0.2,0.2)
                sx,sy = np.ma.masked_where(p_slp>0.05,self.proj.x),np.ma.masked_where(p_slp>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,10*m_slp,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Trend [hPa decade$^{-1}$]')
                for x,y in self.bounds0:
                        pl.plot(x,y,'b',linewidth=3,zorder=10)
                for x,y in self.bounds1:
                        pl.plot(x,y,'r',linewidth=3,zorder=10)
                pl.plot(sx,sy,'k.',alpha=0.5)
		uv        = np.sqrt((10*m_us)**2 + (10*m_vs)**2)
                urot,vrot = self.proj.m.rotate_vector(10*m_us,10*m_vs,self.proj.lon,self.proj.lat,returnxy=False)
                rotmask   = ((p_us>0.05)&(p_vs>0.05))|(uv<0.4)
                urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
		Q = pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],units='width',scale=15,width=0.003,headwidth=4.5,headlength=4,pivot='tail',alpha=0.85)
                qk = pl.quiverkey(Q, 0.6, 1.02, 1, '%s %s' % (1,'m s$^{-1}$ decade$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/slp.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
		cc    = cc + 1
                cseq  = np.arange(0,0.05+0.005,0.005)
                sx,sy = np.ma.masked_where(p_thf_ns>0.05,self.proj.x),np.ma.masked_where(p_thf_ns>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,reg_thf_ns,cseq,cmap=pl.cm.OrRd,extend='max')
                cbar  = pl.colorbar(cf)
                cbar.set_label('[(W m$^{-2}$)$^{-1}$]')
                #pl.plot(sx,sy,'k.',alpha=0.5)
                uv        = np.sqrt(reg_thf_uq**2 + reg_thf_vq**2)
                urot,vrot = self.proj.m.rotate_vector(reg_thf_uq,reg_thf_vq,self.proj.lon,self.proj.lat,returnxy=False)
                rotmask   = ((p_thf_uq>0.05)&(p_thf_vq>0.05))|(uv<0.15)
                urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
                Q = pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],units='width',scale=8,width=0.003,headwidth=4.5,headlength=4,pivot='tail',alpha=0.85)
                qk = pl.quiverkey(Q, 0.4, 1.02, 0.3, '%s %s' % (0.3,'kg s$^{-1}$ m$^{-1}$ (W m$^{-2}$)$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
                self.proj.m.drawparallels([70,80])
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/Thf-N.reg.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
		cc    = cc + 1
                cseq  = np.arange(-0.6,0.6+0.1,0.1)
                sx,sy = np.ma.masked_where(p_ns>0.05,self.proj.x),np.ma.masked_where(p_ns>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,10*m_ns,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Trend [decade$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                uv        = np.sqrt((10*m_uq)**2 + (10*m_vq)**2)
                urot,vrot = self.proj.m.rotate_vector(10*m_uq,10*m_vq,self.proj.lon,self.proj.lat,returnxy=False)
                rotmask   = ((p_uq>0.05)&(p_vq>0.05))|(uv<1)
                urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
                Q = pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],units='width',scale=90,width=0.003,headwidth=4.5,headlength=4,pivot='tail',alpha=0.85)
                qk = pl.quiverkey(Q, 0.6, 1.02, 5, '%s %s' % (5,'kg s$^{-1}$ m$^{-1}$ decade$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/N.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                #pl.figure(cc)
                #cc    = cc + 1
                #cseq  = 13
                #sx,sy = np.ma.masked_where(p_ldf>0.05,self.proj.x),np.ma.masked_where(p_ldf>0.05,self.proj.y)
                #cf    = pl.contourf(self.proj.x,self.proj.y,10*m_ldf,cseq,cmap=pl.cm.coolwarm,extend='both')
                #cbar  = pl.colorbar(cf)
                #cbar.set_label('Downward longwave radiation [W m$^{-2}$ decade$^{-1}$]')
                #pl.plot(sx,sy,'k.',alpha=0.5)
                #self.proj.m.drawparallels([70,80])
                #drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.6')
                #pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/ldf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
		cc    = cc + 1
                cseq  = np.arange(-24,24+4,4)/2.
                cf    = pl.contourf(self.projx,self.projy,10*m_thf[self.iy00:self.iy11+1,self.ix00:self.ix11+1],cseq,cmap=pl.cm.coolwarm,extend='both')
		cl15  = pl.contour(self.projx,self.projy,100*ci2[0:15,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linestyles='dashed',linewidths=2.5,alpha=0.6)
		cl15  = pl.contour(self.projx,self.projy,100*ci2[-15:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Trend [W m$^{-2}$ decade$^{-1}$]')
		#self.proj.m.drawcoastlines(color='0.85',linewidth=1,zorder=11)
		self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3.0,zorder=12)
		pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
		pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
                pl.xlim(self.x[self.ix00],self.x[self.ix11])
                pl.ylim(self.y[self.iy00],self.y[self.iy11])
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/hf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-12,12+2,2)
                cf    = pl.contourf(self.projx,self.projy,10*m_swf[self.iy00:self.iy11+1,self.ix00:self.ix11+1],cseq,cmap=pl.cm.coolwarm,extend='both')
                cl15  = pl.contour(self.projx,self.projy,100*ci2[0:15,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linestyles='dashed',linewidths=2.5,alpha=0.6)
                cl15  = pl.contour(self.projx,self.projy,100*ci2[-15:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Latent flux trend [W m$^{-2}$ decade$^{-1}$]')
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3.0,zorder=12)
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
                pl.xlim(self.x[self.ix00],self.x[self.ix11])
                pl.ylim(self.y[self.iy00],self.y[self.iy11])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/swf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-12,12+2,2)
                cf    = pl.contourf(self.projx,self.projy,10*m_lhf[self.iy00:self.iy11+1,self.ix00:self.ix11+1],cseq,cmap=pl.cm.coolwarm,extend='both')
                cl15  = pl.contour(self.projx,self.projy,100*ci2[0:15,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linestyles='dashed',linewidths=2.5,alpha=0.6)
                cl15  = pl.contour(self.projx,self.projy,100*ci2[-15:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Latent flux trend [W m$^{-2}$ decade$^{-1}$]')
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3.0,zorder=12)
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
                pl.xlim(self.x[self.ix00],self.x[self.ix11])
                pl.ylim(self.y[self.iy00],self.y[self.iy11])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/lhf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-12,12+2,2)
                cf    = pl.contourf(self.projx,self.projy,10*m_shf[self.iy00:self.iy11+1,self.ix00:self.ix11+1],cseq,cmap=pl.cm.coolwarm,extend='both')
                cl15  = pl.contour(self.projx,self.projy,100*ci2[0:15,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linestyles='dashed',linewidths=2.5,alpha=0.6)
                cl15  = pl.contour(self.projx,self.projy,100*ci2[-15:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Sensible flux trend [W m$^{-2}$ decade$^{-1}$]')
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3.0,zorder=12)
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
                pl.xlim(self.x[self.ix00],self.x[self.ix11])
                pl.ylim(self.y[self.iy00],self.y[self.iy11])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/shf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-12,12+2,2)
                cf    = pl.contourf(self.projx,self.projy,10*m_lwf[self.iy00:self.iy11+1,self.ix00:self.ix11+1],cseq,cmap=pl.cm.coolwarm,extend='both')
                cl15  = pl.contour(self.projx,self.projy,100*ci2[0:15,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linestyles='dashed',linewidths=2.5,alpha=0.6)
                cl15  = pl.contour(self.projx,self.projy,100*ci2[-15:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Longwave trend [W m$^{-2}$ decade$^{-1}$]')
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3.0,zorder=12)
                pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
                pl.xlim(self.x[self.ix00],self.x[self.ix11])
                pl.ylim(self.y[self.iy00],self.y[self.iy11])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/lwf.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-24,24+4,4)
                cf    = pl.contourf(self.proj.x,self.proj.y,100*10*m_ci2,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf) 
                cbar.set_label('Trend [% decade$^{-1}$]')
                self.proj.m.drawcoastlines(color='0.85',linewidth=1,zorder=11)
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3,alpha=0.6,zorder=12)
		self.proj.m.drawparallels([70,80])
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/ci.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season2),format='pdf')

		# ORAS4 vectors
		#Hx,Hy,h_lats,h_lons = unpick('/mnt/climstorage/cian/ht.clim.reanal.1980-2016.NDJF.p')
		#Hx,Hy = Hx.data,Hy.data
		#Hx[np.where(Hx==3985.)] = 0
		#Hy[np.where(Hy==3985.)] = 0
		#Hx,Hy = self.proj(Hx,h_lons,h_lats),self.proj(Hy,h_lons,h_lats)
                #urot_h,vrot_h = self.proj.m.rotate_vector(Hx,Hy,self.proj.lon,self.proj.lat,returnxy=False)
                Ux,Uy,h_lats,h_lons = unpick('/mnt/climstorage/cian/vt.clim.reanal.1980-2016.%s.p' % (Season1))
                Ux,Uy = Ux.data,Uy.data 
                Ux,Uy = self.proj(Ux,h_lons,h_lats),self.proj(Uy,h_lons,h_lats)
		UU = np.sqrt(Ux**2 + Uy**2)
                urot_u,vrot_u = self.proj.m.rotate_vector(Ux,Uy,self.proj.lon,self.proj.lat,returnxy=False)
		urot_u,vrot_u = np.ma.masked_where(UU<0.015,urot_u),np.ma.masked_where(UU<0.015,vrot_u)

		# Crop
		Thf_clim    = Thf_clim[self.iy00:self.iy11+1,self.ix00:self.ix11+1]
		urot_u      =   urot_u[self.iy00:self.iy11+1,self.ix00:self.ix11+1]
		vrot_u      =   vrot_u[self.iy00:self.iy11+1,self.ix00:self.ix11+1]

                pl.figure(cc)
                cc    = cc + 1
                cseq  = np.arange(-200,200+25,25)/2.
                cf    = pl.contourf(self.projx,self.projy,Thf_clim,cseq,cmap=pl.cm.coolwarm,extend='both')
		cl15  = pl.contour(self.projx,self.projy,100*ci2[:,self.iy00:self.iy11+1,self.ix00:self.ix11+1].mean(axis=0),[15],colors='k',linewidths=2.5,alpha=0.6)
                cbar  = pl.colorbar(cf)
                cbar.set_label('Climatological surface heat flux [W m$^{-2}$]')
                self.proj.m.fillcontinents(color='0.85',zorder=10)
                for x,y in self.bounds:
                        pl.plot(x,y,'g',linewidth=3,zorder=12)
		#pl.quiver(self.proj.x[::2,::2],self.proj.y[::2,::2],urot_h[::2,::2],vrot_h[::2,::2],units='width',scale=1000,pivot='tail',alpha=0.85)
		Q = pl.quiver(self.projx[::1],self.projy[::1],urot_u[::1,::1],vrot_u[::1,::1],units='width',scale=1,width=0.003,headwidth=4.5,headlength=4,pivot='tail',alpha=0.85)
		#pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot_u[::1,::1],vrot_u[::1,::1],units='inches',scale=2,width=0.005,headwidth=4,headlength=6,pivot='tail',alpha=0.65)
		qk = pl.quiverkey(Q, 0.2, 1.02, 0.05, '%s %s' % (5,'cm s$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
		pl.plot(self.x70,self.y70,'k-',linewidth=1,alpha=0.65,zorder=9)
                pl.plot(self.x80,self.y80,'k-',linewidth=1,alpha=0.65,zorder=9)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.5,color='0.60')
		pl.xlim(self.x[self.ix00],self.x[self.ix11])
		pl.ylim(self.y[self.iy00],self.y[self.iy11])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/model/hf.clim.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

		pl.show()
		
        def sst_OHC(self,YearRange,Season):
		Tice,Tedge,cedge = 273.15-1.0,273.15-0.0,0.15
                years  = range(YearRange[0],YearRange[1]+1,1) 
                # Mask and sea-ice series
                tm0  = np.tile(self.mask,(len(years),1,1))
                ci1  = self.proj(np.array([self.cids.getSeason(year,Season=Season) for year in years]),self.cids.lon,self.cids.lat) 
                tmopen = ((ci1<=cedge)&( tm0==False)&(np.tile(self.lsm,(len( years),1,1))==0))==False  # SST
		tmice  = ((ci1 >cedge)&( tm0==False)&(np.tile(self.lsm,(len( years),1,1))==0))==False  # SST
		# Area of open ocean and ice
		Aopen = (np.ma.masked_array(np.ones((len(years),self.proj.nx,self.proj.ny)),mask=tmopen)*self.gridA).sum(axis=1).sum(axis=1)
		Aice  = self.areaocean - Aopen
                # SST
                sst = self.proj(np.array([self.stds.getSeason(year,Season=Season) for year in years]),self.stds.lon,self.stds.lat)
		#sst[np.where(sst<=Tice)] = Tice
		sst_open = np.ma.masked_array(sst,mask=tmopen)
		#sst_ice  = np.ma.masked_array(sst,mask=tmice)
		sst_ice  = np.ma.masked_array(Tice+np.zeros(sst.shape),mask=tmice)
		# OHC
		ohc_open = 1028*3985*(sst_open*self.bathy*self.gridA).sum(axis=1).sum(axis=1)/1e18
		ohc_ice  = 1028*3985*( sst_ice*self.bathy*self.gridA).sum(axis=1).sum(axis=1)/1e18
		ohc_open_anom = ohc_open - ohc_open.mean()
		ohc_ice_anom  = ohc_ice -  ohc_ice.mean()
		ohc_full_anom = ohc_open_anom + ohc_ice_anom

		# Plot
		mopen,copen,ropen,popen,eropen = stats.linregress(years,ohc_open_anom)
		mice,cice,rice,pice,erice = stats.linregress(years,ohc_ice_anom)
		mfull,cfull,rfull,pfull,erfull = stats.linregress(years,ohc_full_anom)
		lineopen = [mopen*x + copen for x in years]
		lineice  = [mice*x  + cice for x in years]
		linefull = [mfull*x + cfull for x in years]
                fig,ax1  = pl.subplots(num=1)
                #ax2      = ax1.twinx()
		#ax2.plot(years,ohc_open_anom,'r.-',linewidth=2,alpha=0.7)
		#ax2.plot(years,lineopen,'r--',linewidth=1.2,alpha=0.5)
		#ax2.plot(years, ohc_ice_anom,'b.-',linewidth=2,alpha=0.7)
		#ax2.plot(years, lineice,'b--',linewidth=1.2,alpha=0.5)
                ax1.plot(years,ohc_full_anom,'k.-',linewidth=2,alpha=0.7) 
                ax1.plot(years,linefull,'k--',linewidth=1.2,alpha=0.5)
		ax1.set_ylabel('Ocean heat content anomaly [10$^{18}$ J]')
		#ax2.set_ylabel('Ocean heat content anomaly [10$^{18}$ J]')
		ax1.set_xlabel('Year')
		ax1.grid()
		ax1.set_xlim(years[0],years[-1])
		pl.title('slope = %s 10$^{18}$ J year$^{-1}$; p = %s' % (round(mfull,3),round(pfull,5)))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/OHC.%s-%s.%s.pdf' % (years[0],years[-1],Season),format='pdf')

		pl.figure(2)
		cseq  = np.arange(272,278+0.5,0.5)
		cf    = pl.contourf(self.proj.x,self.proj.y,sst.mean(axis=0),cseq,cmap=pl.cm.OrRd,extend='max')

		pl.figure(999)
		cl1   = pl.contour(self.proj.x,self.proj.y,sst.mean(axis=0),[Tice] ,colors='g',linewidth=2)
		cl2   = pl.contour(self.proj.x,self.proj.y,sst.mean(axis=0),[Tedge],colors='k',linewidth=2)
		cl3   = pl.contour(self.proj.x,self.proj.y,ci1.mean(axis=0),[cedge],colors='b',linewidth=2)
		pl.close()
		pl.figure(2)

		PP,MM = [],[]
		# Tedge contour
		paths2 = cl2.collections[0].get_paths()
		for i in range(len(paths2)):
			p   = paths2[i]
			v   = p.vertices
			x,y = v[:,0],v[:,1]
			pp, = pl.plot(x,y,'k',linewidth=2,alpha=0.75)
		PP.append(pp)
		MM.append(r'$T_{sea} = 0^{\circ}C$')
		# Sea-ice contour
                paths3 = cl3.collections[0].get_paths()
		lens   = [len(jj) for jj in paths3]
		x      = lens.index(sorted(lens)[-1])
		p      = paths3[x]
		v      = p.vertices
		x,y    = v[:,0],v[:,1]
		pp,    = pl.plot(x,y,'b',linewidth=1.5,alpha=0.65)
                PP.append(pp)
                MM.append(r'$c_{ice} = 15\%$')
		pl.legend(PP,MM,loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)

		cbar  = pl.colorbar(cf)	
		cbar.set_label('%s [K]' % ('Sea surface temperature'))
		for x,y in self.bounds:
		        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		self.proj.m.drawparallels([70,80])
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		self.proj.m.fillcontinents(color='0.85')
		pl.title('%s' % (Season))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/sst.clim.%s-%s.%s.pdf' % (years[0],years[-1],Season),format='pdf')

		pl.figure(3)
		cseq = [10,25,50,75,100,250,500,750,1000,2500,5000]
		cf   = pl.contourf(self.proj.x,self.proj.y,self.bathy,cseq,cmap=pl.cm.YlGnBu_r,norm=LogNorm())
		cbar = pl.colorbar(cf)
		cbar.set_label('Depth [m]')
		for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		self.proj.m.drawparallels([70,80])
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.fillcontinents(color='0.85')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/bathy.pdf',format='pdf')

		pl.show() 

	def fluxesInBox(self,YearRange,Season1,Season2):
		years  = range(YearRange[0],YearRange[1]+1,1)
		years_ = range(1980,2012+1,1)
		sis    = len(self.reci.getDateList(1997,Season=Season1))*(60*60*6)
		# Mask and sea-ice series
		tm0  = np.tile(self.mask,(len(years),1,1))
		tm0_ = np.tile(self.mask,(len(years_),1,1))
		ci1_ = self.proj(np.array([self.cids.getSeason(year,Season=Season1) for year in years_]),self.cids.lon,self.cids.lat)
		ci1  = self.proj(np.array([self.cids.getSeason(year,Season=Season1) for year in years]),self.cids.lon,self.cids.lat)
		ci2  = self.proj(np.array([self.cids.getSeason(year,Season=Season2) for year in years]),self.cids.lon,self.cids.lat)
		tm1  = (           ( ci1<=1.00)&( tm0==False)&(np.tile(self.lsm,(len( years),1,1))==0))==False	# Surface fluxes
		tm1_ = (           (ci1_<=0.00)&(tm0_==False)&(np.tile(self.lsm,(len(years_),1,1))==0))==False  # SST
		tm2  = ((ci2<=1.0)&( ci2>=0.00)&( tm0==False)&(np.tile(self.lsm,(len( years),1,1))==0))==False	# Sea-ice area

		mask,mask_ = np.zeros((self.proj.nx,self.proj.ny)),ci1_>0.15
                for i in range(self.proj.nx):
                        for j in range(self.proj.ny):
                                if (mask_[:,i,j].any()) or (self.lsm[i,j]==1) or (self.mask[i,j]==True):
                                        mask[i,j] = True
                                else:
                                        mask[i,j] = False

		# Fluxes and SST
                lwf   = self.proj(np.array([self.lwds.getSeason(year,Season=Season1) for year in years]),self.lwds.lon,self.lwds.lat)
		swf   = self.proj(np.array([self.swds.getSeason(year,Season=Season1) for year in years]),self.swds.lon,self.swds.lat)
		shf   = self.proj(np.array([self.shds.getSeason(year,Season=Season1) for year in years]),self.shds.lon,self.shds.lat)
                lhf   = self.proj(np.array([self.lhds.getSeason(year,Season=Season1) for year in years]),self.lhds.lon,self.lhds.lat)
		sst   = self.proj(np.array([self.stds.getSeason(year,Season=Season1) for year in years_]),self.stds.lon,self.stds.lat)
		Thf   = lwf+swf+shf+lhf

		# Compute fluxes in masked region (W m^-2)
		ci = (np.ma.masked_array(ci2,mask=tm2)*self.gridA).sum(axis=1).sum(axis=1)
		lw = (np.ma.masked_array(lwf,mask=tm1)*self.gridA).sum(axis=1).sum(axis=1)#*sis
		sw = (np.ma.masked_array(swf,mask=tm1)*self.gridA).sum(axis=1).sum(axis=1)#*sis
		sh = (np.ma.masked_array(shf,mask=tm1)*self.gridA).sum(axis=1).sum(axis=1)#*sis
		lh = (np.ma.masked_array(lhf,mask=tm1)*self.gridA).sum(axis=1).sum(axis=1)#*sis
		st =  np.ma.masked_array(sst,mask=np.tile(mask,(len(sst),1,1))).mean(axis=1).mean(axis=1)
		Th = lw + sw + sh + lh
		# Compute trends
		lwf,lwm,lwp = self.detrend2d(years,lwf)
		swf,swm,swp = self.detrend2d(years,swf)
		shf,shm,shp = self.detrend2d(years,shf)
		lhf,lhm,lhp = self.detrend2d(years,lhf)
		Thf,Thm,Thp = self.detrend2d(years,Thf)
		sst,stm,stp = self.detrend2d(years_,sst)
		mean_trend  = np.ma.masked_array(Thm,mask=tm1[0]).mean()
		area        = np.ma.masked_array(np.ones((self.proj.nx,self.proj.ny))*self.gridA,mask=tm1[0]).sum()
                m0,c0,r0,p0,er0 = stats.linregress(years,Th)
                m1,c1,r1,p1,er1 = stats.linregress(years,ci)
		m2,c2,r2,p2,er2 = stats.linregress(years_,st)

		print 1000*(m1*1)*334e03

		# Plot
                pl.figure(1)
		cseq  = np.arange(-16,16+2,2)
                #sx,sy = np.ma.masked_where(Thp/2>0.01,self.proj.x),np.ma.masked_where(Thp/2>0.01,self.proj.y)
		sx,sy = np.ma.masked_array(self.proj.x,tm1[0]),np.ma.masked_array(self.proj.y,tm1[0])
                cf    = pl.contourf(self.proj.x,self.proj.y,10*Thm,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                pl.plot(sx,sy,'k.',alpha=0.5)
                cbar.set_label('%s [W m$^{-2}$ decade$^{-1}$]' % ('Total heat flux'))
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                self.proj.m.drawparallels([70,80])
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		print mean_trend
		pl.title('%s x10$^{18}$ J year$^{-1}$' % (round(mean_trend*area*sis/1e18,3)))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/HF.trend.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(2)
                cseq  = np.arange(-0.9,0.9+0.10,0.15) 
                sx,sy = np.ma.masked_array(self.proj.x,mask=mask),np.ma.masked_array(self.proj.y,mask=mask)
                cf    = pl.contourf(self.proj.x,self.proj.y,10*stm,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                pl.plot(sx,sy,'k.',alpha=0.5)
                cbar.set_label('%s [K decade$^{-1}$]' % ('Sea surface temperature'))
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                self.proj.m.drawparallels([70,80])
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                #pl.title()
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/sst.trend.%s-%s.%s.pdf' % (years_[0],years_[-1],Season1),format='pdf')

		pl.figure(3)
		pl.plot(years,lw/1e18,'r--',linewidth=1.5,alpha=0.7,label='longwave')
		pl.plot(years,lw/1e18,'r.',alpha=0.7)
                pl.plot(years,sw/1e18,'g--',linewidth=1.5,alpha=0.7,label='shortwave')
                pl.plot(years,sw/1e18,'g.',alpha=0.7)
		pl.plot(years,sh/1e18,'y.-',linewidth=1.2,alpha=0.4,label='sensible')
		pl.plot(years,lh/1e18,'b.-',linewidth=1.2,alpha=0.4,label='latent')
		pl.grid()
		pl.ylabel('Heat flux [10$^{18}$ J]')
		pl.xlabel('Year')
		pl.xlim(years[0],years[-1])
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.title('%s = %s$_{lw}$ + %s$_{sw}$ + %s$_{sh}$ + %s$_{lh}$ (x10$^{18}$ J year$^{-1}$)' % (round(m0/1e18,2),\
			  round(np.polyfit(years,lw/1e18,1)[0],2),round(np.polyfit(years,sw/1e18,1)[0],2),round(np.polyfit(years,sh/1e18,1)[0],2),round(np.polyfit(years,lh/1e18,1)[0],2)))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/HF.components.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season1),format='pdf')

                pl.figure(4)
		line    = [m2*x + c2 for x in years_]
		HCtrend = (1000)*(3985)*(area)*(200)*(m2)/1e18
                pl.plot(years_,st,'k--',linewidth=1.5,alpha=0.7)
                pl.plot(years_,st,'k.',alpha=0.7)
		pl.plot(years_,line,'k--',linewidth=0.5,alpha=0.6)
                pl.grid()
                pl.ylabel('%s [%s]' % (self.stds.long_name,self.stds.units))
                pl.xlabel('Year')
                pl.xlim(years_[0],years_[-1])
		pl.title('trend = %s 10$^{18}$ J year$^{-1}$; p = %s' % (round(HCtrend,2),round(p2/2,2)))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/sst.%s-%s.%s.pdf' % (years_[0],years_[-1],Season1),format='pdf')

		fig,ax1   = pl.subplots(num=5)
                ax2       = ax1.twinx()
		line0     = np.array([m0*x + c0 for x in years])
		line1     = np.array([m1*x + c1 for x in years])
		Thdt,cidt = Th-line0,ci-line1
		cc        = np.corrcoef(Thdt[:],cidt[:])[0][1]
                ax1.plot(years,Th/1e12,'k.-',linewidth=2,alpha=0.7)
		ax2.plot(years,ci/1e10,'b',linewidth=1.2,alpha=0.35)
                ax1.set_ylabel('Heat flux [TW]')
		ax2.set_ylabel('Sea-ice area [10$^{10}$ m$^{2}$]')
                ax1.set_xlabel('Year')
		ax1.grid()
		pl.title('trend = %s 10$^{18}$ J year$^{-1}$; p = %s; r = %s' % (round(m0/1e18,4),round(p0/2.,2),round(cc,3)))
                ax1.set_xlim(years[0],years[-1])	
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/heatflux/HF_%s.ci_%s.%s-%s.pdf' % (Season1,Season2,YearRange[0],YearRange[1]),format='pdf')

		pl.show()

	def fieldToFlux_regress(self,Field,YearRange,Season):
		years = range(YearRange[0],YearRange[1]+1,1)
		ds = MMDataServer(Field=Field)
		# Field series
		field    = self.proj(np.array([ds.getSeason(year,Season=Season) for year in years]),ds.lon,ds.lat)
		field_t  = np.ma.masked_array(field,mask=np.tile(self.mask,(len(years),1,1))).mean(axis=1).mean(axis=1)/100.
		# Detrend
		m_field,c_field = np.polyfit(years,field_t,1)
		line_f   = np.array([m_field*x + c_field for x in years])
		field_dt = field_t - line_f
		# Surface fluxes
                lwf   = self.proj(np.array([self.lwds.getSeason(year,Season=Season) for year in years]),self.lwds.lon,self.lwds.lat)
                swf   = self.proj(np.array([self.swds.getSeason(year,Season=Season) for year in years]),self.swds.lon,self.swds.lat)
                shf   = self.proj(np.array([self.shds.getSeason(year,Season=Season) for year in years]),self.shds.lon,self.shds.lat)
                lhf   = self.proj(np.array([self.lhds.getSeason(year,Season=Season) for year in years]),self.lhds.lon,self.lhds.lat)
                Thf   = lwf+shf+lhf+swf
		Thf_t = np.ma.masked_array(Thf,mask=np.tile(self.mask,(len(years),1,1))).mean(axis=1).mean(axis=1)/100.
		# Detrend surface flux data
                m_thf,c_thf = np.polyfit(years,Thf_t,1)
                line_thf    = np.array([m_thf*x + c_thf for x in years])
                Thf_dt      = Thf_t - line_thf
		# Detrend
		field,m,p = self.detrend2d(years,field)
		Thf,m,p   = self.detrend2d(years,Thf)
		# Regress Field onto surface fluxes
		m_thf,p_thf = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		m_fld,p_fld = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
                                slope,intercept,r_value,p_value,std_err = stats.linregress(field_dt,Thf[:,i,j])
                                m_thf[i,j] = slope
                                p_thf[i,j] = p_value/2.
                                slope,intercept,r_value,p_value,std_err = stats.linregress(Thf_dt,field[:,i,j])
                                m_fld[i,j] = slope
                                p_fld[i,j] = p_value/2.

                pl.figure(1)
		cseq  = np.arange(-5,5+1,1)
                sx,sy = np.ma.masked_where(p_thf>0.025,self.proj.x),np.ma.masked_where(p_thf>0.025,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,m_thf,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
		cbar.set_label('Reg. coeff. [W m$^{-2}$ hPa$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.title('Total')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/slp.HF.regress.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season))

                pl.figure(2)
                cseq  = 13#np.arange(-5,5+1,1)
                sx,sy = np.ma.masked_where(p_fld>0.025,self.proj.x),np.ma.masked_where(p_fld>0.025,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,m_fld,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [hPa {W m$^{-1}$}$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                for x,y in self.bounds:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                pl.title('Total')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/HF.slp.regress.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season))

		pl.figure(3)
		pl.plot(years,field_t,'k')
		pl.plot(years,line_f,'k--')
		pl.xlabel('Year')
		pl.ylabel('Sea level pressure [hPa]')
		pl.title('trend = %s hPa decade$^{-1}$' % (round(m_field*10,2)))
		pl.xlim(YearRange[0],YearRange[1])
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/slp.series.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season))

		pl.show()

if __name__ == "__main__":

	YearRange = (int(sys.argv[1]),int(sys.argv[2]))
	Season1   = str(sys.argv[3])
	Season2   = str(sys.argv[4])

	#LonRange = (0,360)
	#LatRange = (60,90)      
	LonRange  = (20,60)
	#LonRange = (15,90)
	LatRange  = (70,80)
	BSB = BarentsSeaBudget(	blat     =        65,
				res      =       150,
				LonRange =  LonRange,
				LatRange =  LatRange	)

	BSB.Q_Aice(YearRange,Season1,Season2)
	#BSB.regress(YearRange,Season1)
	#BSB.OHC(YearRange,Season1,plot=True)
	#BSB.surfacePDF(YearRange,Season1,Freq='Monthly')
	#BSB.fluxesInBox(YearRange,Season1,Season2)
	#BSB.sst_OHC(YearRange,Season1)
	#BSB.fieldToFlux_regress(Field='slp',YearRange=YearRange,Season=Season1)
