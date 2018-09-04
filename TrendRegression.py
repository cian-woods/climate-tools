from cmipDataServer import DataServer as cmipDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from scipy import stats
from drawCoastlinesNoRivers import *
from toPick import *
from UnPickle import *
from stipling import *

import matplotlib.pyplot as pl
import sys
import glob

class TrendRegression:

	def __init__(	self,	
			Source    = 'CMCC-CESM',
			ExpType   = 'piControl',
			Field1    = 'tas',
			Field2    = 'psl',
			Season1   = 'DJF',
			Season2   = 'DJF',
			LatRange1 = (75,83),
			LonRange1 = (20,80),
                        LatRange2 = (65,80),
                        LonRange2 = (60,140),
                        LatRange3 = (70,80),
                        LonRange3 = (340,360)	):

                # Lambert Projector
                self.proj = LambertProjector(boundinglat=65,resolution=150.)
		# Attributes
		self.Source    = Source
		self.ExpType   = ExpType
		self.Field1    = Field1
		self.LatRange1 = LatRange1
		self.LonRange1 = LonRange1
                self.Field2    = Field2
                self.LatRange2 = LatRange2
                self.LonRange2 = LonRange2
                self.LatRange3 = LatRange3
                self.LonRange3 = LonRange3
		self.Season1   = Season1
		self.Season2   = Season2
		self.dd        = 0
		if ((Season1 == 'DJF') or (Season2 == 'DJF')) and (Season1 != Season2): self.dd = -1
		print self.dd
		self.bounds1   = self.boxBounds(LonRange1,LatRange1,n=35)
		self.bounds2   = self.boxBounds(LonRange2,LatRange2,n=35)
		self.bounds3   = self.boxBounds(LonRange3,LatRange3,n=35)
		# DataServers
		if Source == 'ERAInt':
			self.ds1 = MMDataServer(Field=self.Field1,LevType='surface_analysis')
			self.ds2 = MMDataServer(Field=self.Field2,LevType='surface_analysis')
		elif Source[0:7] == 'NCAR20C':
                        self.ds1 = MMDataServer(Field=self.Field1,LevType='surface',Source=self.Source)
                        self.ds2 = MMDataServer(Field=self.Field2,LevType='surface',Source=self.Source)
		else:
			self.ds1 = cmipDataServer(Field=self.Field1,LevType='surface',Source=Source,ExpType=ExpType,DataFreq='mon')
			self.ds2 = cmipDataServer(Field=self.Field2,LevType='surface',Source=Source,ExpType=ExpType,DataFreq='mon')
		# Mask for LatRange and LonRange
		plon360    = np.zeros((self.proj.nx,self.proj.ny))
		plon360[np.where(self.proj.lon<0)] = 360
		self.plonf = self.proj.lon + plon360
		self.mask1 = (self.proj.lat<LatRange1[1])&(self.proj.lat>LatRange1[0])&(self.plonf<LonRange1[1])&(self.plonf>LonRange1[0])==False
		self.mask2 = (self.proj.lat<LatRange2[1])&(self.proj.lat>LatRange2[0])&(self.plonf<LonRange2[1])&(self.plonf>LonRange2[0])==False
		#self.mask3 = (self.proj.lat<LatRange3[1])&(self.proj.lat>LatRange3[0])&(self.plonf<LonRange3[1])&(self.plonf>LonRange3[0])==False
                if LonRange3[1]>LonRange3[0]:
                        self.mask3 = (self.proj.lat<LatRange3[1])&(self.proj.lat>LatRange3[0])&(self.plonf<LonRange3[1])&(self.plonf>LonRange3[0])==False
                else:
                        lonmask3   = (self.plonf<LonRange3[1])|(self.plonf>LonRange3[0])
                        self.mask3 = (self.proj.lat<LatRange3[1])&(self.proj.lat>LatRange3[0])&(lonmask3)==False

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

	def reGrid(F,levold,latold,lonold,levnew,latnew,lonnew):
	        # Extend lon axis
	        lonres = np.abs(lonold[1]-lonold[0])
	        dlon   = np.append(lonold,lonold[-1]+lonres)
	        dlon   = np.append(lonold[0]-lonres,dlon)
	        if (type(F)==np.ma.core.MaskedArray) and (F.mask.shape == F.data.shape):
	                F = fill(F.data,invalid=F.mask)
	        F   = extendField(F)
	        F   = interpolateND(F,  dlon,lonnew,axis=2)
	        F   = interpolateND(F,latold,latnew,axis=1)
	        F   = interpolateND(F,levold,levnew,axis=0)
	        return F

	def OLR_Field_regression(self):
		bmap      = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360)
		years     = np.arange(1990,2012+1,1)
		#LatRange = (-15,5)
		#LonRange = (280,310)
		LatRange = (1,15)
		LonRange = (110,130)
                #LatRange = (-8,3)
                #LonRange = (155,190)
		#LatRange  = (50,75)
		#LonRange  = (40,90)
		#LatRange  = (1,15)
		#LonRange  = (80,104)
		# DataServers
		ds2 = MMDataServer(Field='slp',Source='ERAInt')
		ds1 = MMDataServer(Field='ttr',Source='ERAInt')
		# Meshed grid and mask
		lons,lats = np.meshgrid(ds1.lon,ds1.lat)
		xx,yy     = bmap(lons,lats)
		mask      = (lats<LatRange[1])&(lats>LatRange[0])&(lons<LonRange[1])&(lons>LonRange[0])==False
		# Get data and trend
		data1 = np.array([ds1.getSeason(Year=year,Season=self.Season1) for year in years])
		data2 = np.array([ds2.getSeason(Year=year,Season=self.Season2) for year in years+self.dd])
		m1,m2,i1,i2 = [],[],[],[]
		for i in range(data1.shape[1]):	
			m,c = np.polyfit(years,data1[:,i,:],1)
			m1.append(m)
			i1.append(c)
                        m,c = np.polyfit(years+self.dd,data2[:,i,:],1)
                        m2.append(m)
                        i2.append(c)
		m1,m2,i1,i2 = np.array(m1),np.array(m2),np.array(i2),np.array(i2)

                # Plot
                pl.figure(1)
                sx,sy = np.ma.masked_where(mask==True,xx),np.ma.masked_where(mask==True,yy)
                cseq  = 13#np.arange(-12,12+2,2)
                cf    = pl.contourf(xx,yy,m1*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf,orientation='horizontal')
                cbar.set_label('%s trend [%s decade $^{-1}$]' % (ds1.long_name,ds1.units))
                pl.plot(sx[::3,::3],sy[::3,::3],'k.',markersize=3,alpha=0.5)
                bmap.drawcoastlines()
		pl.show()

		# Detrend data
                lines1 = np.array([m1*x + i1 for x in years])
                lines2 = np.array([m2*x + i2 for x in years+self.dd])
		data1  = data1 - lines1
		data2  = data2 - lines2
		# Extract index in lon and lat range
		data1i = np.ma.masked_array(data1,mask=np.tile(mask,(len(years),1,1))).mean(axis=1).mean(axis=1)
		# Regress detrended data
		M1,P1 = np.zeros((ds1.nlat,ds1.nlon)),np.zeros((ds1.nlat,ds1.nlon))
		for i in range(len(ds1.lat)):
			for j in range(len(ds1.lon)):
				slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(data1i ,data2[:,i,j])
				M1[i,j],P1[i,j] = slope1,p_value1/2.
		#M1,P1 = self.proj(M1,ds1.lon,ds1.lat),self.proj(P1,ds1.lon,ds1.lat)

                # Plot
		pl.figure(1)
		sx,sy = np.ma.masked_where(mask==True,xx),np.ma.masked_where(mask==True,yy)
                cseq  = 13#np.arange(-12,12+2,2)
                cf    = pl.contourf(xx,yy,m1*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf,orientation='horizontal')
                cbar.set_label('%s trend [%s decade $^{-1}$]' % (ds1.long_name,ds1.units))
		pl.plot(sx[::3,::3],sy[::3,::3],'k.',markersize=3,alpha=0.5)
                bmap.drawcoastlines()

                pl.figure(2)
                cseq  = 13#np.arange(-12,12+2,2)
                cf    = pl.contourf(xx,yy,m2*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf,orientation='horizontal')
                cbar.set_label('%s trend [%s decade $^{-1}$]' % (ds2.long_name,ds2.units))
                bmap.drawcoastlines()

                # Plot
                pl.figure(3)
                sx,sy = np.ma.masked_where(P1>0.05,xx),np.ma.masked_where(P1>0.05,yy)
                cseq  = 13#np.arange(-12,12+2,2)
                cf    = pl.contourf(xx,yy,M1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf,orientation='horizontal')
                cbar.set_label('Reg. coeff. [%s/%s]' % (ds2.units,ds1.units))
                pl.plot(sx[::3,::3],sy[::3,::3],'k.',markersize=3,alpha=0.5)
                bmap.drawcoastlines()

		#pl.figure(2)
		#sx,sy = np.ma.masked_where(P1>0.05,self.proj.x),np.ma.masked_where(P1>0.05,self.proj.y)
		#cseq  = 13
		#cf    = pl.contourf(self.proj.x,self.proj.y,M1,cseq,cmap=pl.cm.coolwarm,extend='both')
		#cbar  = pl.colorbar(cf)
		#cbar.set_label('Reg. coeff. [%s/%s]' % (ds2.units,ds1.units))
		#pl.plot(sx,sy,'k.',alpha=0.5)
                #self.proj.m.drawparallels([70,80],latmax=90)
                #drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')

                pl.show()

	def timeseries(self,YearRange):
		years   = np.arange(YearRange[0],YearRange[1]+1,1)
		yearse  = np.arange(1980,2016+1,1)
		yearsn  = np.arange(1852,2014+1,1)
		Models  = [g[33:] for g in glob.glob('/mnt/climstorage/cian/historical/*')]
		D,n     = [],0
		# CMIP5 data
		for Source in Models:
			try:
                                ds1   = cmipDataServer(Field='tas',LevType='surface',Source=Source,ExpType=ExpType,DataFreq='mon') 
				data1 = self.proj(np.array([ds1.getSeason(Year=year,Season=self.Season1) for year in years]),ds1.lon,ds1.lat)
				data1 = np.ma.masked_array(data1,mask=np.tile(self.mask1,(len(years),1,1))).mean(axis=1).mean(axis=1)
				D.append(data1)
				n = n + 1
			except:
				pass
		D  = np.array(D)
		Dm = D.mean(axis=0)
		Ds = D.std(axis=0)
		# ERAInt data
		dse   = MMDataServer(Field='T2',LevType='surface_analysis',Source='ERAInt')
		datae = self.proj(np.array([dse.getSeason(Year=year,Season=self.Season1) for year in yearse]),dse.lon,dse.lat)
		datae = np.ma.masked_array(datae,mask=np.tile(self.mask1,(len(yearse),1,1))).mean(axis=1).mean(axis=1)
                # NCAR20C data
                dsn   = MMDataServer(Field='T2',LevType='surface',Source='NCAR20C_V2c')
                datan = self.proj(np.array([dsn.getSeason(Year=year,Season=self.Season1) for year in yearsn]),dsn.lon,dsn.lat)
                datan = np.ma.masked_array(datan,mask=np.tile(self.mask1,(len(yearsn),1,1))).mean(axis=1).mean(axis=1)

		pl.plot(years,Dm,'k-',linewidth=1.75,alpha=0.7)
		pl.plot(years,Dm+Ds,'k--',linewidth=1,alpha=0.6)
		pl.plot(years,Dm-Ds,'k--',linewidth=1,alpha=0.6)
		pl.plot(yearse,datae,'r-',linewidth=1.5,alpha=0.55)
		pl.plot(yearsn,datan,'b-',linewidth=1.5,alpha=0.65)
		pl.grid()
		pl.ylabel('%s [%s]' % (ds1.long_name,ds1.units))
		pl.xlabel('Year')
		pl.title(n)
		pl.show()

	def regression(self,YearRange):
		years   = np.arange(YearRange[0],YearRange[1]+1,1)
		Models  = ['ERAInt'] + [g[33:] for g in glob.glob('/mnt/climstorage/cian/historical/*')]

		M1,M2,T1,T2,n = [],[],[],[],0
		P1,P2,PT1,PT2 = [],[],[],[]
		for Source in Models:
			try:
				# DataServers
				if Source != 'ERAInt':
					ds1 = cmipDataServer(Field='tas',LevType='surface',Source=Source,ExpType=self.ExpType,DataFreq='mon')
					ds2 = cmipDataServer(Field='psl',LevType='surface',Source=Source,ExpType=self.ExpType,DataFreq='mon')
				else:
					ds1   = MMDataServer(Field='T2',Source='ERAInt')
					ds2   = MMDataServer(Field='slp',Source='ERAInt')
				# Get data
				data1   = self.proj(np.array([ds1.getSeason(Year=year,Season=self.Season1) for year in years]),ds1.lon,ds1.lat)
				data2   = self.proj(np.array([ds2.getSeason(Year=year,Season=self.Season2) for year in years+self.dd]),ds2.lon,ds2.lat)
				# Detrend
				data1,trend1,p_value1 = self.detrend2d(years        ,data1)
				data2,trend2,p_value2 = self.detrend2d(years+self.dd,data2)
				# Extract indexes in Lat and Lon ranges
				data1i  = np.ma.masked_array(data1,mask=np.tile(self.mask1,(len(data1),1,1))).mean(axis=1).mean(axis=1)
				data2i  = np.ma.masked_array(data2,mask=np.tile(self.mask2,(len(data2),1,1))).mean(axis=1).mean(axis=1)
				# Regress the various arrays
				m1,p1   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
				m2,p2   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
				for i in range(self.proj.nx):
					for j in range(self.proj.ny):
						slope1i, intercept1i, r_value1i, p_value1i, std_err1i = stats.linregress(data2i        ,data1[:,i,j])
						slope2i, intercept2i, r_value2i, p_value2i, std_err2i = stats.linregress(data1i        ,data2[:,i,j])
						m1[i,j],m2[i,j] = slope1i,slope2i
						p1[i,j],p2[i,j] = p_value1i/2.,p_value2i/2.
				M1.append(m1)
				M2.append(m2)
				P1.append(p1)
				P2.append(p2)
				T1.append(trend1)
				T2.append(trend2)
				PT1.append(p_value1/2.)
				PT2.append(p_value2/2.)
				n = n + 1
			except:
				pass
		# Trends and regressions
		M1,M2,P1,P2,T1,T2,PT1,PT2         = np.array(M1),np.array(M2),np.array(P1),np.array(P2),np.array(T1),np.array(T2),np.array(PT1),np.array(PT2)
		stipx1,stipy1                     = stipling(M1[1:],xx=self.proj.x,yy=self.proj.y,x=None,y=None,thresh=0.8)
		stipx2,stipy2                     = stipling(M2[1:],xx=self.proj.x,yy=self.proj.y,x=None,y=None,thresh=0.8)
		m1,m2,p1,p2,t1,t2,pt1,pt2         = M1[1:].mean(axis=0),M2[1:].mean(axis=0),P1[1:].mean(axis=0),P2[1:].mean(axis=0),T1[1:].mean(axis=0),T2[1:].mean(axis=0),PT1[1:].mean(axis=0),PT2[1:].mean(axis=0)
		em1,em2,ep1,ep2,et1,et2,ept1,ept2 = M1[0],M2[0],P1[0],P2[0],T1[0],T2[0],PT1[0],PT2[0]
		# Mean trend in boxes
		t1i,t2i     = np.ma.masked_array(t1,mask=self.mask1).mean(),np.ma.masked_array(t2,mask=self.mask2).mean()
		et1i,et2i   = np.ma.masked_array(et1,mask=self.mask1).mean(),np.ma.masked_array(et2,mask=self.mask2).mean()

		print 10*t1i,10*t2i
		print 10*et1i,10*et2i

		pl.figure(1)
		cseq  = 13#np.arange(-0.006,0.006+0.001,0.001)	
		cf    = pl.contourf(self.proj.x,self.proj.y,m1,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Reg. coeff. [%s/%s]' % (ds1.units,ds2.units))
                for x,y in self.bounds2:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.plot(stipx1,stipy1,'k.',alpha=0.5)
		self.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.title('all (%s models)' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.reg1.pdf',format='pdf')

		pl.figure(2)
		sx,sy = np.ma.masked_where(ep1>0.05,self.proj.x),np.ma.masked_where(ep1>0.05,self.proj.y)
		cseq  = 13#np.arange(-0.006,0.006+0.001,0.001)
		cf    = pl.contourf(self.proj.x,self.proj.y,em1,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Reg. coeff. [%s/%s]' % (ds1.units,ds2.units))
		#cbar.set_label('Reg. coeff. [m s$^{-1}$/%sPa]' % (round(len(years)*et2i,0)))
		for x,y in self.bounds2:
		        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.plot(sx,sy,'k.',alpha=0.5)
		self.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.title('ERAInt')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/ERAInt.reg1.pdf',format='pdf')

                pl.figure(3)
                cseq  = 13#np.arange(-90,90+15,15) 
                cf    = pl.contourf(self.proj.x,self.proj.y,m2,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s/%s]' % (ds2.units,ds1.units))
                for x,y in self.bounds1:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                pl.plot(stipx2,stipy2,'k.',alpha=0.5)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('all (%s models)' % (n))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.reg2.pdf',format='pdf')

                pl.figure(4)
                cseq  = 13#np.arange(-90,90+15,15)
		sx,sy = np.ma.masked_where(ep2>0.02,self.proj.x),np.ma.masked_where(ep2>0.02,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,em2,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s/%s]' % (ds2.units,ds1.units))
                for x,y in self.bounds1:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.plot(sx,sy,'k.',alpha=0.5)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('ERAInt')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/ERAInt.reg2.pdf',format='pdf')

		pl.figure(5)
                cseq  = 13#np.arange(-3,3+0.5,0.5)
		sx,sy = np.ma.masked_where(ept1>0.05,self.proj.x),np.ma.masked_where(ept1>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,et1*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Observed trend [%s decade$^{-1}$]' % (self.ds1.units))
                for x,y in self.bounds1:
                        pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
		pl.plot(sx,sy,'k.',alpha=0.5)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('ERAInt')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/ERAInt.trend.T2.pdf',format='pdf')

                pl.figure(6)
                cseq = 13#np.arange(-3,3+0.5,0.5)
                cf   = pl.contourf(self.proj.x,self.proj.y,t1*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('Model trend [%s decade$^{-1}$]' % (self.ds1.units))
                for x,y in self.bounds1:
                        pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('all (%s models)' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.trend.tas.pdf',format='pdf')

                pl.figure(7)
                cseq = np.arange(-1.2,1.2+0.2,0.2)
                cf   = pl.contourf(self.proj.x,self.proj.y,t2*10/100.,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('CMIP5 trend [hPa decade$^{-1}$]')
                for x,y in self.bounds2:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('all (%s models)' % (n))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.trend.psl.pdf',format='pdf')

                pl.figure(8)
                cseq = 13#np.arange(-3,3+0.5,0.5)
                cf   = pl.contourf(self.proj.x,self.proj.y,t1*10 + em1*et2i*10,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label('Predicted trend [%s decade$^{-1}$]' % (self.ds1.units))
                for x,y in self.bounds2:
                        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('all (%s models)' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.predictedtrend.pdf',format='pdf')

		pl.figure(9)
		cseq  = np.arange(-1.2,1.2+0.2,0.2)
		cf    = pl.contourf(self.proj.x,self.proj.y,et2*10/100.,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('ERAInt trend [hPa decade$^{-1}$]')
		for x,y in self.bounds2:
		        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		self.proj.m.drawparallels([70,80],latmax=90)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.title('ERAInt')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/ERAInt.trend.slp.pdf',format='pdf')

		pl.show()

	def boxBounds(self,LonRange,LatRange,n=35):
		lons = np.linspace(LonRange[0],LonRange[1],n)
		lats = np.linspace(LatRange[0],LatRange[1],n)
		mer1 = self.proj.m([LonRange[0] for i in range(n)],lats)
		mer2 = self.proj.m([LonRange[1] for i in range(n)],lats)
		zon1 = self.proj.m(lons,[LatRange[0] for i in range(n)])
		zon2 = self.proj.m(lons,[LatRange[1] for i in range(n)])
		return [mer1,mer2,zon1,zon2]

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

	def maskInRanges(self,data1,data2):
		# Takes data in (time)x(lat)x(lon) format and averages
		# the segments within LatRange and LonRange. Returns 1d array	
		data1_ = np.ma.masked_array(data1,mask=np.tile(self.mask1,(len(data1),1,1))).mean(axis=1).mean(axis=1)
		data2_ = np.ma.masked_array(data2,mask=np.tile(self.mask2,(len(data2),1,1))).mean(axis=1).mean(axis=1)
		data3_ = np.ma.masked_array(data2,mask=np.tile(self.mask3,(len(data2),1,1))).mean(axis=1).mean(axis=1)
		return data1_,data2_,data3_

	def getAllMeans(self,years):
		fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (self.Source,self.ExpType,years[0],years[-1],N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                               self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
		data1,data2 = [],[]
		if not os.path.isfile(fname):
			if (self.ExpType == 'piControl') and (self.Source != 'ERAInt'): years = years + self.ds1.Yearmin
			for year in years:
				try:
					d1 = np.ma.masked_array(self.proj(self.ds1.getSeason(Year=year,Season=self.Season1),self.ds1.lon,self.ds1.lat),mask=self.mask1).mean()
				except:
					d1 = -999
				try:
					d2 = np.ma.masked_array(self.proj(self.ds2.getSeason(Year=year+self.dd,Season=self.Season2),self.ds2.lon,self.ds2.lat),mask=self.mask2).mean()
				except:
					d2 = -999
				if (d1 == -999) or (d2 == -999):
					data1.append(-999)
					data2.append(-999)
				else:
					data1.append(d1)
                                        data2.append(d2)
			toPick([years,data1,data2],fname)
		else:
			years,data1,data2 = unpick(fname)
		return data1,data2

	def getTimeSliceTrends(self,year0,N,grad):
		# Opens data array containing all N year windows beginning
		# on all year0s

		try:
			years = np.arange(year0,year0+N,1)
			# Get data
			data1 = self.proj(np.array([self.ds1.getSeason(Year=year,Season=self.Season1) for year in years]),self.ds1.lon,self.ds1.lat)
			data2 = self.proj(np.array([self.ds2.getSeason(Year=year,Season=self.Season2) for year in years+self.dd]),self.ds2.lon,self.ds2.lat)
			#self.plotMap(data1,data2,year0,N)
			# Mask and area average
			data1,data2,data3 = self.maskInRanges(data1,data2)
			# Trends
			m1,m2,m3 = np.polyfit(years,data1,1)[0],np.polyfit(years+self.dd,data2,1)[0],np.polyfit(years+self.dd,data2-data3,1)[0]
			if grad:
				return m1*10,m3*10
			else:
				return m1*10,m2*10
		except:
			return -999,-999

	def getAllTimeSliceTrends(self,year0s,N,grad):
		fname = '/mnt/climstorage/cian/trendfiles/%s/%s.%s.%s-%s.%s.%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (N,self.Source,self.ExpType,year0s[0],year0s[-1],self.Season1,self.Season2,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
														self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
		if grad:
			fname = '/mnt/climstorage/cian/trendfiles/%s/%s.%s.grad.%s-%s.%s.%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.%s-%sN_%s-%sE.p' % (N,self.Source,self.ExpType,year0s[0],year0s[-1],self.Season1,self.Season2,N,self.LatRange1[0],\
																		   self.LatRange1[1],self.LonRange1[0],self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],\
																		   self.LonRange2[0],self.LonRange2[1],self.LatRange3[0],self.LatRange3[1],self.LonRange3[0],self.LonRange3[1])
		print fname
		if not os.path.isfile(fname):
			M1,M2 = [],[]
			for year0 in year0s:
				if (self.ExpType == 'piControl') and (self.Source != 'ERAInt'): year0 = year0 + self.ds1.Yearmin
				m1,m2 = self.getTimeSliceTrends(year0,N,grad)
				print year0,m1,m2
				M1.append(m1)
				M2.append(m2)
			toPick([year0s,M1,M2],fname)
		else:
			print 'File %s aready exits ...' % (fname)
			year0s,M1,M2 = unpick(fname)
		return M1,M2

	def plotMap(self,data1,data2,year0,N):
		os.system('mkdir -p /mnt/climstorage/cian/scripts/figs/trends/%s' % (self.Source))
		trend,pval = self.getTrend(range(len(data1)),data1)
		mx,my      = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
		pl.figure(1)
		cseq = np.arange(-3,3+0.5,0.5)
		cf   = pl.contourf(self.proj.x,self.proj.y,trend,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar = pl.colorbar(cf)
		cl   = pl.contour(self.proj.x,self.proj.y,data1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
		pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
		cbar.set_label('%s [%s decade$^{-1}$]' % (self.ds1.long_name,self.ds1.units))
		for x,y in self.bounds1:
		        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		self.proj.m.drawparallels([70,80])
		pl.title(year0)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/%s.%s.%s.%syears.pdf' % (self.Source,self.Field1,self.Source,year0,N),format='pdf')
		pl.close()

		trend,pval = self.getTrend(range(len(data2)),data2)
		mx,my      = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
		pl.figure(2)
		cseq = np.arange(-480,480+60,60)
		cf   = pl.contourf(self.proj.x,self.proj.y,trend,cseq,cmap=pl.cm.coolwarm,extend='both')
		cbar = pl.colorbar(cf)
		cl   = pl.contour(self.proj.x,self.proj.y,data2.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
		pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
		cbar.set_label('%s [%s decade$^{-1}$]' % (self.ds2.long_name,self.ds2.units))
		for x,y in self.bounds2:
		        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
		pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		self.proj.m.drawparallels([70,80])
		pl.title(year0)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/%s.%s.%s.%syears.pdf' % (self.Source,self.Field2,self.Source,year0,N),format='pdf')
		pl.close()

	def bivar(self,N):
		xedges,yedges = np.linspace(-450,450,12),np.linspace(-2.5,5.5,13)

		Models    = [g[33:] for g in glob.glob('/mnt/climstorage/cian/historical/*')]
		PIM1,PIM2 = [],[]
		HIM1,HIM2 = [],[]
		RCM1,RCM2 = [],[]
		CNM1,CNM2 = [],[]
		for Model in Models:
			try:
				fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'piControl',1,75-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                        	                                                                                        	self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
				year0s,pim1,pim2 = unpick(fname)
				PIM1,PIM2        = PIM1+pim1,PIM2+pim2
			except:
				pass
                        try:
                                fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',1850,2005-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                                self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                                year0s,him1,him2 = unpick(fname)
                                HIM1,HIM2        = HIM1+him1,HIM2+him2
                        except:
                                pass
			try:
       				fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',2005,2095-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
																self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
        			year0s,rcm1,rcm2 = unpick(fname)
        			RCM1,RCM2        = RCM1+rcm1,RCM2+rcm2
			except:
        			pass
			try:
			        fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',1980,2016-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
			                                                                                                        self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
			        year0s,cnm1,cnm2 = unpick(fname)
			        CNM1,CNM2        = CNM1+cnm1,CNM2+cnm2
			except:
			        pass
		# Coordinates
		piM1,piM2 = [PIM1[i] for i in range(len(PIM1)) if (PIM1[i] != -999) and (PIM2[i] != -999)],[PIM2[i] for i in range(len(PIM1)) if (PIM1[i] != -999) and (PIM2[i] != -999)]
		hiM1,hiM2 = [HIM1[i] for i in range(len(HIM1)) if (HIM1[i] != -999) and (HIM2[i] != -999)],[HIM2[i] for i in range(len(HIM1)) if (HIM1[i] != -999) and (HIM2[i] != -999)]
		rcM1,rcM2 = [RCM1[i] for i in range(len(RCM1)) if (RCM1[i] != -999) and (RCM2[i] != -999)],[RCM2[i] for i in range(len(RCM1)) if (RCM1[i] != -999) and (RCM2[i] != -999)]
		cnM1,cnM2 = [CNM1[i] for i in range(len(CNM1)) if (CNM1[i] != -999) and (CNM2[i] != -999)],[CNM2[i] for i in range(len(CNM1)) if (CNM1[i] != -999) and (CNM2[i] != -999)]
		# Bivariate PDFs
		xe,ye,Hpi = self.bivarPDF(piM2,piM1,xedges,yedges,norm=True)
		xe,ye,Hhi = self.bivarPDF(hiM2,hiM1,xedges,yedges,norm=True)
		xe,ye,Hrc = self.bivarPDF(rcM2,rcM1,xedges,yedges,norm=True)
		Hpi       = np.rollaxis(Hpi,1,0)
		Hhi       = np.rollaxis(Hhi,1,0)
		Hrc       = np.rollaxis(Hrc,1,0)
		# ERAInt points
		fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % ('ERAInt','historical',1980,2016-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                 self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
		y0s,erm1,erm2 = unpick(fname)
		# NCAR20C points
                fname = '/mnt/climstorage/cian/trendfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % ('NCAR20C_V2c','historical',1850,2016-N+1,N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                 self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                y0s,nrm1,nrm2 = unpick(fname)
                # Regression lines
		xs = np.linspace(-600,600,20)
		slopeNR, interceptNR, r_valueNR, p_valueNR, std_errNR = stats.linregress(nrm2,nrm1)
                slopeER, interceptER, r_valueER, p_valueER, std_errER = stats.linregress(erm2,erm1)
		slopePI, interceptPI, r_valuePI, p_valuePI, std_errPI = stats.linregress(piM2,piM1)
		slopeCN, interceptCN, r_valueCN, p_valueCN, std_errCN = stats.linregress(cnM2,cnM1)
		lineNR = [x*slopeNR + interceptNR for x in xs]
		lineER = [x*slopeER + interceptER for x in xs]
		linePI = [x*slopePI + interceptPI for x in xs]
		lineCN = [x*slopeCN + interceptCN for x in xs]

		# Plot bivariate PDFs
		pl.figure(1)
		cf   = pl.contourf(xe,ye,Hhi,np.arange(0.0002,0.0009+0.00015,0.00015),cmap=pl.cm.OrRd,extend='both')
		cf.cmap.set_under('w')
		cbar = pl.colorbar(cf)
		cbar.set_label('Frequency')
		cl1  = pl.contour(xe,ye,Hpi,5,colors='g',linewidths=0.6,alpha=0.75,linestyles='--')
		cl2  = pl.contour(xe,ye,Hrc,5,colors='k',linewidths=0.6,alpha=0.3)
	        manual_locations2 = [(-227,3.3),(-165,3.3),(-125,2.3),(-45,2.5),(40,1.3)]
		manual_locations1 = [(0,-2),(-50,-1.4),(-100,-1),(-50,-0.8)]
	        pl.clabel(cl1,fmt='%0.4f',colors='g',fontsize=9,manual=manual_locations1)
		pl.clabel(cl2,fmt='%0.4f',colors='k',fontsize=9,manual=manual_locations2)
		# ERAInt points
		pl.plot(erm2,erm1,'k.-',markersize=7,alpha=0.65)
	        years = range(1980,2016-N+2,1)
	        for label, x, y in zip([years[i] for i in [0,-1]], [erm2[i] for i in [0,-1]], [erm1[i] for i in [0,-1]]):
	                plt.annotate(   label,xy=(x, y),xytext=(0, 0),color='k',
	                                textcoords='offset points', ha='center', va='bottom',alpha=0.65)
		# NCAR20C points
		pl.plot(nrm2,nrm1,'k.',markersize=6,alpha=0.6)
                #years = range(1980,2014-23+2,1)
                #for label, x, y in zip([years[i] for i in [0,-1]], [nrm2[i] for i in [0,-1]], [nrm1[i] for i in [0,-1]]):
                #        plt.annotate(   label,xy=(x, y),xytext=(0, 0),color='k',
                #                        textcoords='offset points', ha='center', va='bottom',alpha=0.65)
		# 1980 - 2016 model points
		pl.plot(cnM2,cnM1,'b.',markersize=4,alpha=0.3)
		# Crosses
		pl.plot(np.mean(erm2),np.mean(erm1),'k+',markersize=10,alpha=0.65,mew=2)
		pl.plot(np.mean(hiM2),np.mean(hiM1),'r+',markersize=10,alpha=0.95,mew=2)
		pl.plot(np.mean(rcM2),np.mean(rcM1),'k+',markersize=10,alpha=0.50,mew=2)
		pl.plot(np.mean(piM2),np.mean(piM1),'g+',markersize=10,alpha=1.00,mew=2)
		pl.plot(np.mean(cnM2),np.mean(cnM1),'b+',markersize=10,alpha=0.85,mew=2)
		# Regression lines
		pl.plot(xs,lineNR,'k--',linewidth=1,alpha=0.35)
		pl.plot(xs,lineER,'k--',linewidth=1,alpha=0.35)
		pl.plot(xs,linePI,'g--',linewidth=1,alpha=0.35)
		#pl.plot(xs,lineCN,'b--',linewidth=1,alpha=0.35)
		# Attributes
		pl.grid(color='b')
		pl.xlim(-440,440)
		pl.ylim(-2,4)
		#pl.xlim(-530,530)
		#pl.ylim(-2.75,5.25)
		pl.xlabel('Sea level pressure [Pa decade$^{-1}$]')
		pl.ylabel('Surface temperature [K decade$^{-1}$]')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/bivar.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.pdf' %\
			  (N,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1]), format='pdf')


		# Marginal PDF
		#Hps,edges = np.histogram(piM2,bins=np.arange(-400,400+50,50),range=None,normed=True,weights=None,density=None)
		Hhs,edges = np.histogram(hiM2,bins=np.arange(-900,900+45,45),range=None,normed=True,weights=None,density=None)
		#Hrs,edges = np.histogram(rcM2,bins=np.arange(-400,400+50,50),range=None,normed=True,weights=None,density=None)

		edges = (edges[0:-1] + edges[1:])/2.
		pl.figure(2)
		pl.plot(edges,Hhs,'k',linewidth=2,alpha=0.85)
		#pl.plot(edges[0],Hhs[0],'ko')
		#pl.plot(edges[-1],Hhs[-1],'ko')
		pl.ylabel('Normalised frequency')
		pl.xlabel('Sea-level pressure trend [Pa decade$^{-1}$]')
		pl.grid()

		pl.show()

        def bivar1(self):
                xedges,yedges = np.linspace(99500,103500,12),np.linspace(230,275,13)

                # ERAInt points
                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % ('ERAInt','historical',1990,2016,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                 self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                y0s,erm1,erm2 = unpick(fname)
                # NCAR20C points
                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % ('NCAR20C_V2c','historical',1990,2014,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                 self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                y0s,nrm1,nrm2 = unpick(fname)
		# Model points
                Models    = [g[33:] for g in glob.glob('/mnt/climstorage/cian/historical/*')]
                PIM1,PIM2 = [],[]
                HIM1,HIM2 = [],[]
                RCM1,RCM2 = [],[]
                CNM1,CNM2 = [],[]
                for Model in Models:
                        try:
                                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'piControl',1,75,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                                self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                                year0s,pim1,pim2    = unpick(fname)
				pim1,pim2           = np.array(pim1),np.array(pim2)
				xs1,xs2             = np.where(pim1!=-999)[0],np.where(pim2!=-999)[0]
				pim1[xs1],pim2[xs2] = pim1[xs1]-pim1[xs1].mean()+260,pim2[xs2]-pim2[xs2].mean()+101500
                                PIM1,PIM2           = PIM1+list(pim1),PIM2+list(pim2)
                        except:
                                pass
                        try:
                                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',1850,2005,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                                self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                                year0s,him1,him2 = unpick(fname)
                                HIM1,HIM2        = HIM1+him1,HIM2+him2
                        except:
                                pass
                        try:
                                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',2005,2100,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                                self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                                year0s,rcm1,rcm2 = unpick(fname)
                                RCM1,RCM2        = RCM1+rcm1,RCM2+rcm2
                        except:
                                pass
                        try:
                                fname = '/mnt/climstorage/cian/meanfiles/%s.%s.%s-%s.%syears.%s-%sN_%s-%sE_%s-%sN_%s-%sE.p' % (Model,'rcp85',2060,2085,23,self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],\
                                                                                                                                self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1])
                                year0s,cnm1,cnm2    = unpick(fname)
				cnm1,cnm2           = np.array(cnm1),np.array(cnm2)
				xs1,xs2             = np.where(cnm1!=-999)[0],np.where(cnm2!=-999)[0]
				cnm1[xs1],cnm2[xs2] = cnm1[xs1]-cnm1[xs1].mean()+240,cnm2[xs2]-cnm2[xs2].mean()+101500
                                CNM1,CNM2           = CNM1+list(cnm1),CNM2+list(cnm2)
                        except:
                                pass
                # Coordinates
                piM1,piM2 = [PIM1[i] for i in range(len(PIM1)) if (PIM1[i] != -999) and (PIM2[i] != -999)],[PIM2[i] for i in range(len(PIM1)) if (PIM1[i] != -999) and (PIM2[i] != -999)]
                hiM1,hiM2 = [HIM1[i] for i in range(len(HIM1)) if (HIM1[i] != -999) and (HIM2[i] != -999)],[HIM2[i] for i in range(len(HIM1)) if (HIM1[i] != -999) and (HIM2[i] != -999)]
                rcM1,rcM2 = [RCM1[i] for i in range(len(RCM1)) if (RCM1[i] != -999) and (RCM2[i] != -999)],[RCM2[i] for i in range(len(RCM1)) if (RCM1[i] != -999) and (RCM2[i] != -999)]
                cnM1,cnM2 = [CNM1[i] for i in range(len(CNM1)) if (CNM1[i] != -999) and (CNM2[i] != -999)],[CNM2[i] for i in range(len(CNM1)) if (CNM1[i] != -999) and (CNM2[i] != -999)]
                # Bivariate PDFs
                xe,ye,Hpi = self.bivarPDF(piM2,piM1,xedges,yedges,norm=True)
                xe,ye,Hhi = self.bivarPDF(hiM2,hiM1,xedges,yedges,norm=True)
                xe,ye,Hrc = self.bivarPDF(rcM2,rcM1,xedges,yedges,norm=True)
                Hpi       = np.rollaxis(Hpi,1,0)
                Hhi       = np.rollaxis(Hhi,1,0)
                Hrc       = np.rollaxis(Hrc,1,0)
                # Regression lines
                xs = np.linspace(99000,104000,20)
                slopeER, interceptER, r_valueER, p_valueER, std_errER = stats.linregress(erm2,erm1)
		slopeNR, interceptNR, r_valueNR, p_valueNR, std_errNR = stats.linregress(nrm2,nrm1)
                slopePI, interceptPI, r_valuePI, p_valuePI, std_errPI = stats.linregress(piM2,piM1)
                slopeCN, interceptCN, r_valueCN, p_valueCN, std_errCN = stats.linregress(cnM2,cnM1)
		print 'ERAInt; p = %s' % (p_valueER/2.)
		print 'NCAR20C; p = %s' % (p_valueNR/2.)
		print 'piControl; p = %s' % (p_valuePI/2.)
		print 'contemp; p = %s' % (p_valueCN/2.)
                lineER = [x*slopeER + interceptER for x in xs]
		lineNR = [x*slopeNR + interceptNR for x in xs]
                linePI = [x*slopePI + interceptPI for x in xs]
                lineCN = [x*slopeCN + interceptCN for x in xs]
                # Plot bivariate PDFs
                #cf   = pl.contourf(xe,ye,Hhi,np.arange(0.000008,0.000056+0.000008,0.000008),cmap=pl.cm.OrRd,extend='both')
                #cf.cmap.set_under('w')
                #cbar = pl.colorbar(cf)
                #cbar.set_label('Frequency')
                cl1  = pl.contour(xe,ye,Hpi,5,colors='g',linewidths=0.6,alpha=0.75,linestyles='--')
                #cl2  = pl.contour(xe,ye,Hrc,5,colors='k',linewidths=0.6,alpha=0.3)
                #manual_locations2 = [(-227,3.3),(-165,3.3),(-125,2.3),(-45,2.5),(40,1.3)]
                manual_locations1 = [(0,-2),(-50,-1.4),(-100,-1),(-50,-0.8)]
                pl.clabel(cl1,fmt='%0.4f',colors='g',fontsize=9,manual=manual_locations1)
                #pl.clabel(cl2,fmt='%0.4f',colors='k',fontsize=9,manual=manual_locations2)
                # ERAInt and NCAR20C points
                pl.plot(erm2,erm1,'k.',markersize=7,alpha=0.65)
		pl.plot(nrm2,nrm1,'r.',markersize=5,alpha=0.4)
                # 1980 - 2016 model points
                pl.plot(cnM2,cnM1,'b.',markersize=4,alpha=0.3)
                # Crosses
                pl.plot(np.mean(erm2),np.mean(erm1),'k+',markersize=10,alpha=0.65,mew=2)
                #pl.plot(np.mean(hiM2),np.mean(hiM1),'r+',markersize=10,alpha=0.95,mew=2)
                #pl.plot(np.mean(rcM2),np.mean(rcM1),'k+',markersize=10,alpha=0.50,mew=2)
                pl.plot(np.mean(piM2),np.mean(piM1),'g+',markersize=10,alpha=1.00,mew=2)
                pl.plot(np.mean(cnM2),np.mean(cnM1),'b+',markersize=10,alpha=0.85,mew=2)
                # Regression lines
                pl.plot(xs,lineER,'k--',linewidth=1,alpha=0.35)
		pl.plot(xs,lineNR,'r--',linewidth=1,alpha=0.35)
                pl.plot(xs,linePI,'g--',linewidth=1,alpha=0.35)
                pl.plot(xs,lineCN,'b--',linewidth=1,alpha=0.35)
                # Attributes
                pl.grid(color='b')
                pl.xlim(99500,103500)
                pl.ylim(230,272)
		pl.xticks([99500,100000,100500,101000,101500,102000,102500,103000,103500],[995,1000,1005,1010,1015,1020,1025,1030,1035])
                pl.xlabel('Sea level pressure [hPa]')
                pl.ylabel('Surface temperature [K]')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/bivar1.23years.%s-%sN_%s-%sE_%s-%sN_%s-%sE.pdf' %\
                          (self.LatRange1[0],self.LatRange1[1],self.LonRange1[0],self.LonRange1[1],self.LatRange2[0],self.LatRange2[1],self.LonRange2[0],self.LonRange2[1]), format='pdf')
                pl.show()

	def bivarPDF(self,x,y,xedges,yedges,norm=True):
		H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=None,normed=norm,weights=None)
		xedges,yedges   = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
		return xedges,yedges,H

if __name__ == "__main__":

	Source    = str(sys.argv[1])
	ExpType   = str(sys.argv[2])
	Field1    = str(sys.argv[3])
	Field2    = str(sys.argv[4])
	Season1   = str(sys.argv[5])
	Season2   = str(sys.argv[6])
	YearRange = (int(sys.argv[7]),int(sys.argv[8]))
	N         = int(sys.argv[9])

	tr = TrendRegression(	Source    =       Source,
				ExpType   =      ExpType,
				Field1    =       Field1,
				Field2    =       Field2,
				Season1   =      Season1,
				Season2   = 	 Season2,
				LatRange1 =      (75,83),
				LonRange1 =      (20,80),
	                        LatRange2 =      (65,75),
	                        LonRange2 =     (60,100),
				LatRange3 =      (70,80),
                                LonRange3 =     (340,15))

	#tr.timeseries(YearRange)
	#tr.OLR_Field_regression()
	tr.regression(YearRange)
	#tr.bivar(N)
	#tr.bivar1()
	#M1,M2 = tr.getAllTimeSliceTrends(range(YearRange[0],YearRange[1]-N+2,1),N,grad=True)
	#R1,R2 = tr.getAllMeans(np.arange(YearRange[0],YearRange[1]+1,1))

