import matplotlib.pyplot as pl
import glob

from toPick import *
from UnPickle import *
from cmipDataServer import DataServer as cmipDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from drawCoastlinesNoRivers import *
from scipy import stats
from UnPickle import *
from stipling import *

class FieldChange:

	def __init__(	self,
			Source  =  'ERAInt'	,
			ExpType =  'historical' ,
			Field1  = ['slp','psl'] ,
			Field2  = ['T2','tas']		):
		# Attributes
		self.proj    = LambertProjector(boundinglat=50,resolution=200.)
		self.Source  = Source
		ExpType0     = ExpType
		if ExpType == 'historical': ExpType0 = 'rcp85'
		self.ExpType = ExpType
		if Source != 'ERAInt':
		        self.d1 = cmipDataServer(Field=Field1[1],LevType='surface',Source=Source,ExpType=ExpType0,DataFreq='mon')
		        self.d2 = cmipDataServer(Field=Field2[1],LevType='surface',Source=Source,ExpType=ExpType0,DataFreq='mon')
		if Source == 'ERAInt':
		        self.d1 = MMDataServer(Field=Field1[0])
		        self.d2 = MMDataServer(Field=Field2[0])
			self.d3 = MMDataServer(Field='V',LevRange=(800,800))
			self.d4 = MMDataServer(Field='U',LevRange=(800,800))

	def trendFile(self,year0s=range(2007,2077+1,1),N=23):

		fname = 'trends/%s.%s.trends.p' % (self.Source,self.ExpType)
		print fname
		if not os.path.isfile(fname):
			T,P   = [],[]
			# First year0 trends
			years = range(year0s[0],year0s[0]+N,1)
			data  = self.proj(np.array([self.d2.getSeason(Year=year,Season='DJF') for year in years]),self.d2.lon,self.d2.lat)
			trend,pval = self.getTrend(years,data)
			T.append(trend)
			P.append(pval)
			print years[0],trend.shape
			for year0 in year0s[1:]:
				snap = self.proj(self.d2.getSeason(Year=years[-1]+1,Season='DJF'),self.d2.lon,self.d2.lat)
				data = np.append(data[1:,:,:],snap[np.newaxis,:,:],axis=0)
				years      = years[1:] + [years[-1] + 1]
				trend,pval = self.getTrend(years,data)
				T.append(trend)
				P.append(pval)
				print years[0],trend.shape
			T,P = np.array(T),np.array(P)
			toPick([year0s,T,P],fname)
		else:
			year0s,T,P = unpick(fname)
		return year0s,T,P

		"""
		#trendm,trends,pval = T.mean(axis=0),T.std(axis=0),P.mean(axis=0)
		for i in range(len(T)):
			trendm,pval = T[i],P[i]
			trends      = T.mean(axis=0)
			mx,my       = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
			pl.figure(1)
			cseq = np.arange(-3,3+0.5,0.5)
			cf   = pl.contourf(self.proj.x,self.proj.y,trendm,cseq,cmap=pl.cm.coolwarm,extend='both')
			cbar = pl.colorbar(cf)
			cl   = pl.contour(self.proj.x,self.proj.y,trends,7,colors='k',linewidths=0.6,alpha=0.3)
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
			cbar.set_label('%s [%s decade$^{-1}$]' % (self.d2.long_name,self.d2.units))
			#for x,y in bounds1:
			#        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
			#for x,y in bounds2:
			#       pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
			drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
			self.proj.m.drawparallels([70,80])
			pl.title('%s: %s' % (self.Source,year0s[i]))
			pl.show()
		"""
        def bivarPDF(self,x,y,xedges,yedges,norm=True):
                H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=None,normed=norm,weights=None)
                xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
                return xedges,yedges,H

	def boxBounds(self,m,LatRange,LonRange,n=35):
		lons = np.linspace(LonRange[0],LonRange[1],n)
		lats = np.linspace(LatRange[0],LatRange[1],n)
		mer1 = m([LonRange[0] for i in range(n)],lats)
		mer2 = m([LonRange[1] for i in range(n)],lats)
		zon1 = m(lons,[LatRange[0] for i in range(n)])
		zon2 = m(lons,[LatRange[1] for i in range(n)])
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

	def detrend(self,years,field1,field2,detrend=True):
		m1,c1 = np.polyfit(years,field1,1)
		m2,c2 = np.polyfit(years,field2,1)
		if detrend == True:
			line1  = np.array([m1*x + c1 for x in years])
			line2  = np.array([m2*x + c2 for x in years])
			field1 = field1 - line1
			field2 = field2 - line2
		r_pear = np.corrcoef(field1,field2)[0][1]
		return field1,field2,10*m1,10*m2,r_pear

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

	def modelRegressions(self,YearRange=(1980,2016),LatRange1=[65,80],LonRange1=[60,140],LatRange2=[55,75],LonRange2=[0,360],LatRange3=[80,90],LonRange3=[180,360],plot=True):
                # Input data
                years    = range(YearRange[0],YearRange[1]+1,1)
		bounds1  = self.boxBounds(self.proj.m,LatRange1,LonRange1,35)
		bounds2  = self.boxBounds(self.proj.m,LatRange2,LonRange2,35)
		bounds3  = self.boxBounds(self.proj.m,LatRange3,LonRange3,35)
		fftimed1 = self.proj(np.array([self.d1.getSeason(Year=year,Season='DJF') for year in years]),self.d1.lon,self.d1.lat)
                fftimed2 = self.proj(np.array([self.d2.getSeason(Year=year,Season='DJF') for year in years]),self.d2.lon,self.d2.lat)
		# Detrend data
		fftimed1_dt,trend1,pvaltrend1 = self.detrend2d(years,fftimed1)
		fftimed2_dt,trend2,pvaltrend2 = self.detrend2d(years,fftimed2)
		# Mask data
                fftime1_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange1[0],fftimed1_dt)
                fftime1_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange1[1],fftime1_dt)
                fftime1_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange1[0],fftime1_dt)
                fftime1_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange1[1],fftime1_dt)
		fftime1_dt = fftime1_dt.mean(axis=-1).mean(axis=-1)
                fftime1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange1[0],fftimed1)
                fftime1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange1[1],fftime1)
                fftime1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange1[0],fftime1)
                fftime1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange1[1],fftime1)
                fftime1 = fftime1.mean(axis=-1).mean(axis=-1)
		ttrend1,ptrend1 = np.polyfit(years,fftime1,1)

                fftime2_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange2[0],fftimed2_dt)
                fftime2_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange2[1],fftime2_dt)
                fftime2_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange2[0],fftime2_dt)
                fftime2_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange2[1],fftime2_dt)
                fftime2_dt = fftime2_dt.mean(axis=-1).mean(axis=-1)
                fftime2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange2[0],fftimed2)
                fftime2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange2[1],fftime2)
                fftime2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange2[0],fftime2)
                fftime2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange2[1],fftime2)
                fftime2 = fftime2.mean(axis=-1).mean(axis=-1)
                ttrend2,ptrend2 = np.polyfit(years,fftime2,1)

                fftime3_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange3[0],fftimed1_dt)
                fftime3_dt = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange3[1],fftime3_dt)
                fftime3_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange3[0],fftime3_dt)
                fftime3_dt = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange3[1],fftime3_dt)
                fftime3_dt = fftime3_dt.mean(axis=-1).mean(axis=-1)
                fftime3 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange3[0],fftimed1)
                fftime3 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange3[1],fftime3)
                fftime3 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange3[0],fftime3)
                fftime3 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange3[1],fftime3)
                fftime3 = fftime3.mean(axis=-1).mean(axis=-1)
                ttrend3,ptrend3 = np.polyfit(years,fftime3,1)

		grad = fftime1 - fftime3

		# Regress
		regcoeff,pval   = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
				slope, intercept, r_value, p_value, std_err = stats.linregress(fftime1_dt,fftimed2_dt[:,i,j])
				#slope, intercept, r_value, p_value, std_err = stats.linregress(grad,fftimed2_dt[:,i,j])
				regcoeff[i,j] = slope
				pval[i,j]     = p_value/2
		predtrend2 = ttrend1*regcoeff

                ptrend2 = np.ma.masked_where(self.proj.lat<LatRange2[0],predtrend2)
                ptrend2 = np.ma.masked_where(self.proj.lat>LatRange2[1],ptrend2)
                ptrend2 = np.ma.masked_where(self.proj.lon<LonRange2[0],ptrend2)
                ptrend2 = np.ma.masked_where(self.proj.lon>LonRange2[1],ptrend2)
		ptrend2 = ptrend2.mean()

		if plot == True:
			# Plot
			pl.figure(44)
			mx,my = np.ma.masked_where(pval>0.05,self.proj.x),np.ma.masked_where(pval>0.05,self.proj.y)
			cseq = 13#np.arange(-0.5,0.5+0.1,0.1)
			cf   = pl.contourf(self.proj.x,self.proj.y,regcoeff,cseq,cmap=pl.cm.coolwarm,extend='both')
			cbar = pl.colorbar(cf)
			cl   = pl.contour(self.proj.x,self.proj.y,fftimed1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
			cbar.set_label('Reg. coeff. [%s/%s]' % (self.d2.units,self.d1.units))
			for x,y in bounds1:
			        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
			drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
			self.proj.m.drawparallels([70,80])
			pl.title(self.Source)
			pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/tas.psl.reg.%s-%s.DJF.pdf' % (self.Source,YearRange[0],YearRange[1]),format='pdf')
			#pl.close()

                        pl.figure(45)
			mx,my = np.ma.masked_where(pvaltrend1>0.1,self.proj.x),np.ma.masked_where(pvaltrend1>0.1,self.proj.y)
                        cseq  = 13#np.arange(-480,480+60,60)
                        cf    = pl.contourf(self.proj.x,self.proj.y,10*trend1,cseq,cmap=pl.cm.coolwarm,extend='both')
                        cbar  = pl.colorbar(cf)
                        cl    = pl.contour(self.proj.x,self.proj.y,fftimed1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
                        pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
                        cbar.set_label('%s trend [%s decade$^{-1}$]' % (self.d1.long_name,self.d1.units))
                        for x,y in bounds1:
                                pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                        for x,y in bounds3:
                                pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                        pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
                        drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                        self.proj.m.drawparallels([70,80])
                        pl.title(self.Source)
                        pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/psl.trend.%s-%s.DJF.pdf' % (self.Source,YearRange[0],YearRange[1]),format='pdf')
			#pl.close()

                        pl.figure(46)
                        mx,my = np.ma.masked_where(pvaltrend2>0.1,self.proj.x),np.ma.masked_where(pvaltrend2>0.1,self.proj.y)
                        cseq  = 13#np.arange(-3,3+0.5,0.5)
                        cf    = pl.contourf(self.proj.x,self.proj.y,10*trend2,cseq,cmap=pl.cm.coolwarm,extend='both')
                        cbar  = pl.colorbar(cf)
                        cl    = pl.contour(self.proj.x,self.proj.y,fftimed1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
                        pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
                        cbar.set_label('%s trend [%s decade$^{-1}$]' % (self.d2.long_name,self.d2.units))
                        for x,y in bounds2:
                                pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                        pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
                        drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                        self.proj.m.drawparallels([70,80])
                        pl.title('%s: %s %s decade$^{-1}$' % (self.Source,round(10*ttrend2,2),self.d2.units))
                        pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/tas.trend.%s-%s.DJF.pdf' % (self.Source,YearRange[0],YearRange[1]),format='pdf')
			#pl.close()

			pl.figure(47)
                        cseq  = 13#np.arange(-2,2+0.25,0.25)
                        cf    = pl.contourf(self.proj.x,self.proj.y,10*predtrend2,cseq,cmap=pl.cm.coolwarm,extend='both')
                        cbar  = pl.colorbar(cf)
                        #cl    = pl.contour(self.proj.x,self.proj.y,10*trend2,7,colors='k',linewidths=0.6,alpha=0.3)
                        #pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
                        cbar.set_label('%s predicted trend [%s decade$^{-1}$]' % (self.d2.long_name,self.d2.units))
                        for x,y in bounds2:
                                pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                        #pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
                        drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                        self.proj.m.drawparallels([70,80])
			pl.title('%s: %s %s decade$^{-1}$' % (self.Source,round(10*ptrend2,2),self.d2.units))
                        pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/%s/tas.predtrend.%s-%s.DJF.pdf' % (self.Source,YearRange[0],YearRange[1]),format='pdf')
			#pl.close()

			return regcoeff,pval,trend1,trend2,10*ttrend2,10*ptrend2
		else:
			return regcoeff,pval,trend1,trend2,10*ttrend2,10*ptrend2

	def modelTrends(self,year0=1980,N=37,LatRange1=[65,80],LonRange1=[60,140],LatRange2=[75,83],LonRange2=[20,80],landmask=False,scalef=1,plot=True):
		# Input data
		years     = range(year0,year0+N,1)
		bounds1   = self.boxBounds(self.proj.m,LatRange1,LonRange1,35)
		bounds2   = self.boxBounds(self.proj.m,LatRange2,LonRange2,35)

		if self.Source == 'ERAInt':
			# Bandpass filtered slp
			slpstd,lon,lat = unpick('/mnt/climstorage/cian/scripts/synopfiles/ERAInt.msl.2-6days-pass.std.1980-2016.DJF.p')
			slpyears       = range(1980,2016+1,1)
			slpstd         = np.array([slpstd[i] for i in range(len(slpstd)) if slpyears[i] in years])
			slpyears       = [slpyears[i] for i in range(len(slpyears)) if slpyears[i] in years]
			slpstd         = self.proj(slpstd,lon,lat)
	
		# Time series in boxes
		fftimed1 = scalef*self.proj(np.array([self.d1.getSeason(Year=year,Season='DJF') for year in years]),self.d1.lon,self.d1.lat)
		fftimed2 = scalef*self.proj(np.array([self.d2.getSeason(Year=year,Season='DJF') for year in years]),self.d2.lon,self.d2.lat)
		if self.Source == 'ERAInt':
			fftimedv = scalef*self.proj(np.array([self.d3.getSeason(Year=year,Season='DJF') for year in years]),self.d3.lon,self.d3.lat).squeeze()
			fftimedu = scalef*self.proj(np.array([self.d4.getSeason(Year=year,Season='DJF') for year in years]),self.d4.lon,self.d4.lat).squeeze()
			fftimedU = np.sqrt(fftimedu**2 + fftimedv**2)
			fftimedd = np.sin(fftimedv/fftimedU)		

		# Data segments
		if landmask:
			lsm = proj(self.d1.lsm,self.d1.lon,self.d1.lat)
			lsm[np.where(lsm>=0.5)] = 1
			lsm[np.where(lsm<0.5)]  = 0
			fftimed1 = np.ma.masked_where(lsm==1,fftimed1)
			fftimed2 = np.ma.masked_where(lsm==1,fftimed2)
		
		fftime1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange1[0],fftimed1)
		fftime1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange1[1],fftime1)
		fftime1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange1[0],fftime1)
		fftime1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange1[1],fftime1)

		fftime2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange2[0],fftimed1)
		fftime2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange2[1],fftime2)
		fftime2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange2[0],fftime2)
		fftime2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange2[1],fftime2)

		fftime3 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange2[0],fftimed2)
		fftime3 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange2[1],fftime3)
		fftime3 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange2[0],fftime3)
		fftime3 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange2[1],fftime3)

                fftime4 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange1[0],fftimed2)
                fftime4 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange1[1],fftime4)
                fftime4 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange1[0],fftime4)
                fftime4 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange1[1],fftime4)

		if self.Source == 'ERAInt':
                	fftime5 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))<LatRange2[0],fftimedd)
                	fftime5 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(years),1,1))>LatRange2[1],fftime5)
                	fftime5 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))<LonRange2[0],fftime5)
                	fftime5 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(years),1,1))>LonRange2[1],fftime5)

			slpstd1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(slpyears),1,1))<LatRange1[0],slpstd)
			slpstd1 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(slpyears),1,1))>LatRange1[1],slpstd1)
			slpstd1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(slpyears),1,1))<LonRange1[0],slpstd1)
			slpstd1 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(slpyears),1,1))>LonRange1[1],slpstd1)

			slpstd2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(slpyears),1,1))<LatRange2[0],slpstd)
			slpstd2 = np.ma.masked_where(np.tile(self.proj.lat[np.newaxis,:,:],(len(slpyears),1,1))>LatRange2[1],slpstd2)
			slpstd2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(slpyears),1,1))<LonRange2[0],slpstd2)
			slpstd2 = np.ma.masked_where(np.tile(self.proj.lon[np.newaxis,:,:],(len(slpyears),1,1))>LonRange2[1],slpstd2)

		fftime1 = fftime1.mean(axis=1).mean(axis=1)
		fftime2 = fftime2.mean(axis=1).mean(axis=1)
		fftime3 = fftime3.mean(axis=1).mean(axis=1)
		fftime4 = fftime4.mean(axis=1).mean(axis=1)
		if self.Source == 'ERAInt':
			fftime5 = fftime5.mean(axis=1).mean(axis=1)
			slpstd1 = slpstd1.mean(axis=1).mean(axis=1)
			slpstd2 = slpstd2.mean(axis=1).mean(axis=1)

		"""
		a0,a1,a2 = np.polyfit(years,fftime3,2)
		line     = [a0*(x**2) + a1*x + a2 for x in years]
		pl.plot(years,fftime3,'k.',markersize=6,alpha=0.75)
		pl.plot(years,line,'k--',linewidth=1,alpha=0.35,label=r'$%s{x^{2}} + %sx + %s$' % (round(a0,4),round(a1,2),round(a2,2)))
		pl.xlabel('Year')
                pl.ylabel('%s [%s]' % (self.d2.long_name,self.d2.units))
		pl.legend(loc=0)
		pl.show()
		"""

		# Detrend two timeseries for plotting
		field1,field2,m1,m2,r_pear = self.detrend(years,fftime1,fftime3,detrend=True)

		"""
		slope, intercept, r_value, p_value, std_err = stats.linregress(fftime1,fftime3)
		print slope
                pred = slope*m1
                actu = m2
                print pred,actu
		"""
		"""
		pl.plot(field1,field2,'k+',markersize=9,alpha=0.75,mew=2)	
		pl.xlabel('%s [%s]' % (self.d1.long_name,self.d1.units))
		pl.ylabel('%s [%s]' % (self.d2.long_name,self.d2.units))
		pl.grid()
		pl.title(r'$\frac{dF}{dG}$ = %s; r = %s; m1 = %s; m2 = %s' % (round(slope,4),round(r_pear,3),round(m1,1),round(m2,1)))
		pl.show()
		"""

		if plot == False:
			return m1,m2,fftime1-fftime1.mean(),fftime3-fftime3.mean()
		else:
			# Standardise
			field1 = (field1-field1.mean())/field1.std()
			field2 = (field2-field2.mean())/field2.std()

			trend,pval = self.getTrend(years,fftimed1)
			mx,my      = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
			pl.figure(1)
			cseq = np.arange(-480,480+60,60)
			cf   = pl.contourf(self.proj.x,self.proj.y,trend,cseq,cmap=pl.cm.coolwarm,extend='both')
			cbar = pl.colorbar(cf)
			cl   = pl.contour(self.proj.x,self.proj.y,fftimed1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
			cbar.set_label('%s [%s decade$^{-1}$]' % (self.d1.long_name,self.d1.units))
			for x,y in bounds1:
				pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
			#for x,y in bounds2:
			#	pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
			drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
			self.proj.m.drawparallels([70,80])
			pl.savefig('/mnt/climstorage/cian/scripts/figs/ERAInt/trends/slp.%s.%syears.pdf' % (year0,N),format='pdf')
			pl.close()

			trend,pval = self.getTrend(slpyears,slpstd)
			mx,my      = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
			pl.figure(2)
			cseq = np.arange(-40,40+5,5)
			cf   = pl.contourf(self.proj.x,self.proj.y,trend,cseq,cmap=pl.cm.coolwarm,extend='both')
			cbar = pl.colorbar(cf)
			cl   = pl.contour(self.proj.x,self.proj.y,slpstd.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
			cbar.set_label('%s [%s decade$^{-1}$]' % ('Standard deviation of\nmean sea level pressure','Pa'))
			#for x,y in bounds1:
			#	pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
			for x,y in bounds2:
				pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
			drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
			self.proj.m.drawparallels([70,80])
			pl.savefig('/mnt/climstorage/cian/scripts/figs/ERAInt/trends/slpvar.%s.%syears.pdf' % (year0,N),format='pdf')
			pl.close()

			trend,pval = self.getTrend(years,fftimed2)
			mx,my      = np.ma.masked_where(pval>0.1,self.proj.x),np.ma.masked_where(pval>0.1,self.proj.y)
			pl.figure(3)
			if self.d2.Field != 'ci': cseq = np.arange(-3,3+0.5,0.5)
			if self.d2.Field == 'ci': cseq = np.arange(-0.3,0.3+0.05,0.05)
			cf   = pl.contourf(self.proj.x,self.proj.y,trend,cseq,cmap=pl.cm.coolwarm,extend='both')
			cbar = pl.colorbar(cf)
			cl   = pl.contour(self.proj.x,self.proj.y,fftimed2.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
			pl.clabel(cl,fmt='%0.0f',colors='k',fontsize=9)
			cbar.set_label('%s [%s decade$^{-1}$]' % (self.d2.long_name,self.d2.units))
			#for x,y in bounds1:
			#	pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
			for x,y in bounds2:
				pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.3,markersize=4)
			drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
			self.proj.m.drawparallels([70,80])
			pl.savefig('/mnt/climstorage/cian/scripts/figs/ERAInt/trends/T2.%s.%syears.pdf' % (year0,N),format='pdf')
			pl.close()

                        trendv,pvalv = self.getTrend(years,fftimedv)
			trendu,pvalu = self.getTrend(years,fftimedu)
			trendd,pvald = self.getTrend(years,fftimedd)
			trendU,pvalU = self.getTrend(years,fftimedU)
			urot,vrot    = self.proj.m.rotate_vector(trendu,trendv,self.proj.lon,self.proj.lat,returnxy=False)
			urot,vrot    = np.ma.masked_where(pvald>0.1,urot),np.ma.masked_where(pvald>0.1,vrot)
			mx,my        = np.ma.masked_where(pvalU>0.1,self.proj.x),np.ma.masked_where(pvalU>0.1,self.proj.y)
                        pl.figure(4)
			cseq = np.arange(-1,1+0.2,0.2)
                        cf   = pl.contourf(self.proj.x,self.proj.y,trendU,cseq,cmap=pl.cm.coolwarm,extend='both')
                        cbar = pl.colorbar(cf)  
                        cbar.set_label('%s [%s decade$^{-1}$]' % ('Wind speed',self.d3.units))
                	Q  = pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],units='inches',scale=3,\
                       			scale_units='inches',headwidth=3.5,headlength=3.5,headaxislength=3.5,pivot='tail',alpha=0.6)
                	qk = pl.quiverkey(Q, 0.2, 1.02, 1, '%s %s decade$^{-1}$' % (1,self.d3.units), labelpos='W',fontproperties={'weight': 'bold'})
                        #for x,y in bounds1:
                        #        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
                        for x,y in bounds2:
                                pl.plot(x,y,'r',linewidth=1.5,alpha=0.6)
			pl.plot(mx,my,'k.',alpha=0.8,markersize=4)
                        drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                        self.proj.m.drawparallels([70,80])
			pl.savefig('/mnt/climstorage/cian/scripts/figs/ERAInt/trends/u.%s.%syears.pdf' % (year0,N),format='pdf')
			pl.close()

			pl.figure(5)
			pl.plot(years,field1,'b',linewidth=1.5,alpha=0.7)
			pl.plot(years,field2,'r',linewidth=1.5,alpha=0.7)
			pl.ylabel('Standard deviations')
			pl.xlabel('Year')
			pl.title('r = %s; m1 = %s; m2 = %s' % (round(r_pear,3),round(m1,1),round(m2,1)))
			pl.grid()	
			pl.savefig('/mnt/climstorage/cian/scripts/figs/ERAInt/trends/time_series.%s.%syears.pdf' % (year0,N),format='pdf')
			pl.close()
	
			return m1,m2,fftime1-fftime1.mean(),fftime3-fftime3.mean()

	def testReg(self):
		years   = range(1850,2075+15,15)
		S       = {}
		Models_ = [g[9:] for g in glob.glob('../rcp85/*')]
		for Model in Models_:
			S[Model] = np.zeros(len(years))
		for jj in range(len(years)):
			year         = years[jj]
			Models,m1,m2 = unpick('trends/histo/65-80N_60-140E_75-83N_20-80E/model_trends_rcp85_23years_%d-%d.p' % (year,year+15))
			slopes       = [stats.linregress(m1[i],m2[i])[0] for i in range(len(m1)) if m1[i] != []]
			models       = [Models[i] for i in range(len(m1)) if m1[i] != []]
			for ii in range(len(models)):
				S[models[ii]][jj] = slopes[ii]
		ss = np.array([S[i] for i in S])
		ss = np.ma.masked_where(ss==0,ss)
		for i in range(len(ss)):
			pl.plot(years,ss[i],linewidth=0.75,alpha=0.5,label=S.keys()[i])
			pl.grid()
			pl.ylabel('SLP to T2 trend regression (23 years)')
			pl.xlabel('Year')
			pl.legend()
			pl.xlim(1850,2100)
			pl.ylim(-0.02,0.02)
			pl.show()

if __name__ == "__main__":

	from toPick import *
	from UnPickle import *
	import glob

	"""
	LatRange1,LonRange1 = [55,75],[30,100]
        LatRange2,LonRange2 = [75,83],[20,80]
	Dir                 = '%s-%sN_%s-%sE_%s-%sN_%s-%sE' % (LatRange1[0],LatRange1[1],LonRange1[0],LonRange1[1],LatRange2[0],LatRange2[1],LonRange2[0],LonRange2[1])
	fname               = '/mnt/climstorage/cian/scripts/trends/histo/%s/ERAInt_predtrends.%syears_%s-%s.p' % (Dir,23,1980,1994)
	if not os.path.isfile(fname):
		print 'Making %s ...' % (fname)
		fc    = FieldChange(Source='ERAInt',ExpType='rcp85')
		TT,PT = [],[]
		for i in range(1980,1994+1,1):
			YearRange = (i,i+23-1)
			regcoeffE,pvalE,trend1,trend2,ttrend2,ptrend2 = fc.modelRegressions(YearRange=YearRange,LatRange1=LatRange1,LonRange1=LonRange1,LatRange2=LatRange2,LonRange2=LonRange2,plot=False)
			TT.append(ttrend2)
			PT.append(ptrend2)
		toPick([TT,PT],fname)
	else:
		print 'File %s already exists ...' % (fname)
	"""

	YearRange           = (1980,2016)
        LatRange2,LonRange2 = [55,75],[45,90]
        LatRange1,LonRange1 = [75,83],[20,80]
	LatRange3,LonRange3 = [80,90],[-180,0]
	Field2,Field1       = ['slp','psl'],['T2','tas']
	fc                  = FieldChange(Source='ERAInt',ExpType='rcp85',Field1=Field1,Field2=Field2)
	regcoeffE,pvalE,trend1,trend2,ttrend2,ptrend2 = fc.modelRegressions(YearRange=YearRange,LatRange1=LatRange1,LonRange1=LonRange1,LatRange2=LatRange2,LonRange2=LonRange2,LatRange3=LatRange3,LonRange3=LonRange3,plot=True)
	pl.show()
	Models              = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*')]
	R,P,T1,T2,models    = [],[],[],[],[]
	for Model in Models[0:3]:
		try:
			fc = FieldChange(Source=Model,ExpType='rcp85',Field1=Field1,Field2=Field2)
			regcoeff,pval,trend1,trend2,ttrend2,ptrend2 = fc.modelRegressions(YearRange=YearRange,LatRange1=LatRange1,LonRange1=LonRange1,LatRange2=LatRange2,LonRange2=LonRange3,LatRange3=LatRange3,LonRange3=LonRange3,plot=False)
			R.append(regcoeff)
			P.append(pval)
			T1.append(trend1)
			T2.append(trend2)
			models.append(Model)
		except:
			pass
	R,P,T1,T2  = np.array(R),np.array(P),np.array(T1),np.array(T2)
	Ra         = R - regcoeffE[np.newaxis,:,:]
	Ram        = Ra.mean(axis=0)
	Rm,T1m,T2m = R.mean(axis=0),T1.mean(axis=0),T2.mean(axis=0)
	sx,sy      = stipling(R ,xx=fc.proj.x,yy=fc.proj.y,x=None,y=None,thresh=0.8)
	sxa,sya    = stipling(Ra,xx=fc.proj.x,yy=fc.proj.y,x=None,y=None,thresh=0.8)

	# Plot
	pl.figure(1)
	cseq = 13#np.arange(-0.5,0.5+0.1,0.1)
	cf   = pl.contourf(fc.proj.x,fc.proj.y,Rm,cseq,cmap=pl.cm.coolwarm,extend='both')
	cbar = pl.colorbar(cf)
	#cl   = pl.contour(self.proj.x,self.proj.y,fftimed1.mean(axis=0),7,colors='k',linewidths=0.6,alpha=0.3)
	#pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=9)
	cbar.set_label('Reg. coeff. [%s/%s]' % (fc.d2.units,fc.d1.units))
	#for x,y in bounds:
	#        pl.plot(x,y,'b',linewidth=1.5,alpha=0.6)
	pl.plot(sx,sy,'k.',alpha=0.3,markersize=4)
	drawCoastlinesNoRivers(fc.proj.m,linewidth=0.6,color='0.2')
	fc.proj.m.drawparallels([70,80])
	pl.title('all: %s models' % (len(models)))
	pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.%s.%s.reg.full.%s-%s.DJF.pdf' % (Field2[0],Field1[0],YearRange[0],YearRange[1]),format='pdf')

        pl.figure(2)
        cseq = 13#np.arange(-0.3,0.3+0.05,0.05)
        cf   = pl.contourf(fc.proj.x,fc.proj.y,Ram,cseq,cmap=pl.cm.coolwarm,extend='both')
        cbar = pl.colorbar(cf)
        cbar.set_label('Reg. coeff. [%s/%s]' % (fc.d2.units,fc.d1.units))
        pl.plot(sxa,sya,'k.',alpha=0.3,markersize=4)
        drawCoastlinesNoRivers(fc.proj.m,linewidth=0.6,color='0.2')
        fc.proj.m.drawparallels([70,80])
        pl.title('all: %s models' % (len(models)))
        pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.%s.%s.reg.bias.%s-%s.DJF.pdf' % (Field2[0],Field1[0],YearRange[0],YearRange[1]),format='pdf')

	pl.figure(3)
	cseq  = 13#np.arange(-3,3+0.5,0.5)
	cf    = pl.contourf(fc.proj.x,fc.proj.y,10*T2m,cseq,cmap=pl.cm.coolwarm,extend='both')
	cbar  = pl.colorbar(cf)
	cbar.set_label('%s trend [%s decade$^{-1}$]' % (fc.d2.long_name,fc.d2.units))
	drawCoastlinesNoRivers(fc.proj.m,linewidth=0.6,color='0.2')
	fc.proj.m.drawparallels([70,80])
	pl.title('all: %s models' % (len(models)))
	pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/all.%s.trend.%s-%s.DJF.pdf' % (Field2[0],YearRange[0],YearRange[1]),format='pdf')

	pl.show()
 
	"""
        fc            = FieldChange(Source='ERAInt',ExpType='rcp85')
	LatRange1,LonRange1 = [55,75],[30,100]
        LatRange2,LonRange2 = [75,83],[20,80]
	xedges_m,yedges_m   = np.linspace(-450,450,12),np.linspace(-2.5,5.5,12)
	xedges_f,yedges_f   = np.linspace(-600,600,15),np.linspace(-5,5,15)
	#xedges_f,yedges_f   = np.linspace(100000,103000,15),np.linspace(230,280,15)
	Dir                 = '%s-%sN_%s-%sE_%s-%sN_%s-%sE' % (LatRange1[0],LatRange1[1],LonRange1[0],LonRange1[1],LatRange2[0],LatRange2[1],LonRange2[0],LonRange2[1])
	N                   = 23
	# Open Files
	TT,PT              = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/ERAInt_predtrends.%syears_%s-%s.p' 		% (Dir,N,1980,1994))
	modhst,m1hst,m2hst = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_trends_historical_%syears_%s-%s.p' 	% (Dir,N,1850,1983))
	modrcp,m1rcp,m2rcp = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_trends_rcp85_%syears_%s-%s.p' 		% (Dir,N,2005,2078))
	modpic,m1pic,m2pic = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_trends_piControl_%syears_%s-%s.p' 	% (Dir,N,1,53))
	mod,m1con,m2con    = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_trends_rcp85_%syears_%s-%s.p' 		% (Dir,N,1980,1994))
        modhst,f1hst,f2hst = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_clims_historical_%syears_%s-%s.p' 	% (Dir,N,1850,1983))
        modrcp,f1rcp,f2rcp = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_clims_rcp85_%syears_%s-%s.p' 		% (Dir,N,2005,2078))
        modpic,f1pic,f2pic = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_clims_piControl_%syears_%s-%s.p' 	% (Dir,N,1,53))
        mod,f1con,f2con    = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_clims_rcp85_%syears_%s-%s.p' 		% (Dir,N,1980,1994))
	modera,m1ERA,m2ERA = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_trends_ERAInt_%syears_%s-%s.p' 	% (Dir,N,1980,1994))
	modera,f1ERA,f2ERA = unpick('/mnt/climstorage/cian/scripts/trends/histo/%s/model_clims_ERAInt_%syears_%s-%s.p' 		% (Dir,N,1980,1994))
	# Sort arrays
	m1hst,m2hst = [j for i in m1hst for j in i],[j for i in m2hst for j in i]
	m1rcp,m2rcp = [j for i in m1rcp for j in i],[j for i in m2rcp for j in i]
	m1pic,m2pic = [j for i in m1pic for j in i],[j for i in m2pic for j in i]
	m1ERA,m2ERA = [j for i in m1ERA for j in i],[j for i in m2ERA for j in i]
	m1con,m2con = [j for i in m1con for j in i],[j for i in m2con for j in i]

        f1hst,f2hst = np.array([j[0] for i in f1hst for j in i]),np.array([j[0] for i in f2hst for j in i])
        f1rcp,f2rcp = np.array([j[0] for i in f1rcp for j in i]),np.array([j[0] for i in f2rcp for j in i])
        f1pic,f2pic = np.array([j[0] for i in f1pic for j in i]),np.array([j[0] for i in f2pic for j in i])
        f1ERA,f2ERA = np.array([j[0] for i in f1ERA for j in i]),np.array([j[0] for i in f2ERA for j in i])
	f1con,f2con = np.array([j[0] for i in f1con for j in i]),np.array([j[0] for i in f2con for j in i])

        #f1hst,f2hst = [k for i in f1hst for j in i for k in j],[k for i in f2hst for j in i for k in j]
        #f1rcp,f2rcp = [k for i in f1rcp for j in i for k in j],[k for i in f2rcp for j in i for k in j]
        #f1pic,f2pic = [k for i in f1pic for j in i for k in j],[k for i in f2pic for j in i for k in j]
        #f1ERA,f2ERA = [k for i in f1ERA for j in i for k in j],[k for i in f2ERA for j in i for k in j]
        #f1con,f2con = [k for i in f1con for j in i for k in j],[k for i in f2con for j in i for k in j]

	#fzhst = list(set(zip(f1hst,f2hst)))
	#f1hst,f2hst = np.array([i[0] for i in fzhst]),np.array([i[1] for i in fzhst])
        #fzrcp = list(set(zip(f1rcp,f2rcp)))
        #f1rcp,f2rcp = np.array([i[0] for i in fzrcp]),np.array([i[1] for i in fzrcp])
        #fzpic = list(set(zip(f1pic,f2pic)))
        #f1pic,f2pic = np.array([i[0] for i in fzpic]),np.array([i[1] for i in fzpic])
        #fzERA = list(set(zip(f1ERA,f2ERA)))
        #f1ERA,f2ERA = np.array([i[0] for i in fzERA]),np.array([i[1] for i in fzERA])
        #fzcon = list(set(zip(f1con,f2con)))
        #f1con,f2con = np.array([i[0] for i in fzcon]),np.array([i[1] for i in fzcon])

	print len(f1hst),len(f1rcp),len(f1pic),len(f1ERA),len(f1con)

	# Extract means and add ERAInt means
	#f1hst,f2hst = f1hst-f1hst.mean(),f2hst-f2hst.mean()
        #f1rcp,f2rcp = f1rcp-f1rcp.mean(),f2rcp-f2rcp.mean()
        #f1pic,f2pic = f1pic-f1pic.mean(),f2pic-f2pic.mean()
        #f1ERA,f2ERA = f1ERA-f1ERA.mean(),f2ERA-f2ERA.mean()
        #f1con,f2con = f1con-f1con.mean(),f2con-f2con.mean()
	# Compute bivar PDFs
	xem,yem,H1m = fc.bivarPDF(m1hst,m2hst,xedges_m,yedges_m,norm=True)
	xem,yem,H2m = fc.bivarPDF(m1rcp,m2rcp,xedges_m,yedges_m,norm=True)
	xem,yem,H4m = fc.bivarPDF(m1pic,m2pic,xedges_m,yedges_m,norm=True)
	H1m,H2m,H4m  = np.rollaxis(H1m,1,0),np.rollaxis(H2m,1,0),np.rollaxis(H4m,1,0)
        xef,yef,H1f = fc.bivarPDF(f1hst,f2hst,xedges_f,yedges_f,norm=True)
        xef,yef,H2f = fc.bivarPDF(f1rcp,f2rcp,xedges_f,yedges_f,norm=True)
        xef,yef,H4f = fc.bivarPDF(f1pic,f2pic,xedges_f,yedges_f,norm=True)
        H1f,H2f,H4f  = np.rollaxis(H1f,1,0),np.rollaxis(H2f,1,0),np.rollaxis(H4f,1,0)

	# Plot
	pl.figure(1)
	cf   = pl.contourf(xem,yem,H1m,np.arange(0.0002,0.0009+0.00015,0.00015),cmap=pl.cm.OrRd,extend='both')
	cf.cmap.set_under('w')
	cl1  = pl.contour(xem,yem,H2m,5,colors='k',linewidths=0.6,alpha=0.3)
	cl2  = pl.contour(xem,yem,H4m,4,colors='g',linewidths=0.6,alpha=0.75,linestyles='--')
	cbar = pl.colorbar(cf)
	cbar.set_label('Frequency')
	#manual_locations1 = [(-227,3.3),(-165,3.3),(-125,2.3),(-45,2.5),(40,1.3)]
	pl.clabel(cl1,fmt='%0.4f',colors='k',fontsize=9)#,manual=manual_locations1)
	#manual_locations2 = [(0,-2),(-50,-1.4),(-100,-1),(-50,-0.8)]
	pl.clabel(cl2,fmt='%0.4f',colors='g',fontsize=9)#,manual=manual_locations2)
	pl.plot(m1ERA,m2ERA,'k.-',markersize=7,alpha=0.65)
	pl.plot(m1ERA,PT,'c.-',markersize=7,alpha=0.65)

	pl.plot(np.mean(m1ERA),np.mean(m2ERA)  ,'k+',markersize=10,alpha=0.65,mew=2)
	pl.plot(np.mean(m1hst),np.mean(m2hst)  ,'r+',markersize=10,alpha=0.95,mew=2)
	pl.plot(np.mean(m1rcp),np.mean(m2rcp)  ,'k+',markersize=10,alpha=0.40,mew=2)
	pl.plot(np.mean(m1pic),np.mean(m2pic)  ,'g+',markersize=10,alpha=1.00,mew=2)
	pl.plot(np.mean(m1con),np.mean(m2con)  ,'b+',markersize=10,alpha=0.85,mew=2)
	pl.plot(np.mean(m1ERA),np.mean(PT)     ,'c+',markersize=10,alpha=0.85,mew=2)

	pl.plot(m1con,m2con,'b.',markersize=4,alpha=0.3)
	slopePIC, interceptPIC, r_value, p_value, std_err = stats.linregress(m1pic,m2pic)
	slopeERA, interceptERA, r_value, p_value, std_err = stats.linregress(m1ERA,m2ERA)

	xs = np.linspace(-600,600,10)
	linePIC = [slopePIC*x + interceptPIC for x in xs]
	lineERA = [slopeERA*x + interceptERA for x in xs]

	pl.plot(xs,linePIC,'g--',linewidth=1,alpha=0.35)
	pl.plot(xs,lineERA,'k--',linewidth=1,alpha=0.35)

	years = range(1980,2016-N+2,1)
	for label, x, y in zip([years[i] for i in [0,-1]], [m1ERA[i] for i in [0,-1]], [m2ERA[i] for i in [0,-1]]):
    		plt.annotate(   label,xy=(x, y),xytext=(0, 0),color='k',
        			textcoords='offset points', ha='center', va='bottom',alpha=0.65)
	pl.grid(color='b')
	pl.xlim(-475,475)
	pl.ylim(-2.75,5.25)
	pl.xlabel('Sea level pressure [Pa decade$^{-1}$]')
	pl.ylabel('Surface temperature [K decade$^{-1}$]')
	pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/bivar.trends.%syears.pdf' % (N),format='pdf')

        pl.figure(2)
        cf   = pl.contourf(xef,yef,H1f,5,cmap=pl.cm.OrRd,extend='both')
        cf.cmap.set_under('w')
        cl1  = pl.contour(xef,yef,H2f,5,colors='k',linewidths=0.6,alpha=0.3)
        cl2  = pl.contour(xef,yef,H4f,4,colors='g',linewidths=0.6,alpha=0.75,linestyles='--')
        cbar = pl.colorbar(cf)
        cbar.set_label('Frequency')
        #manual_locations1 = [(-227,3.3),(-165,3.3),(-125,2.3),(-45,2.5),(40,1.3)]
        pl.clabel(cl1,fmt='%0.5f',colors='k',fontsize=9)#,manual=manual_locations1)
        #manual_locations2 = [(0,-2),(-50,-1.4),(-100,-1),(-50,-0.8)]
        pl.clabel(cl2,fmt='%0.5f',colors='g',fontsize=9)#,manual=manual_locations2)
        pl.plot(f1ERA,f2ERA,'k.',markersize=7,alpha=0.65)

        pl.plot(np.mean(f1ERA),np.mean(f2ERA)  ,'k+',markersize=10,alpha=0.65,mew=2)
        pl.plot(np.mean(f1hst),np.mean(f2hst)  ,'r+',markersize=10,alpha=0.95,mew=2)
        pl.plot(np.mean(f1rcp),np.mean(f2rcp)  ,'k+',markersize=10,alpha=0.40,mew=2)
        pl.plot(np.mean(f1pic),np.mean(f2pic)  ,'g+',markersize=10,alpha=1.00,mew=2)
        pl.plot(np.mean(f1con),np.mean(f2con)  ,'b+',markersize=10,alpha=0.85,mew=2)

        pl.plot(f1con,f2con,'b.',markersize=4,alpha=0.3)
        slopePIC, interceptPIC, r_value, p_value, std_err = stats.linregress(f1pic,f2pic)
        slopeERA, interceptERA, r_value, p_value, std_err = stats.linregress(f1ERA,f2ERA)

        xs = xedges_f
        linePIC = [slopePIC*x + interceptPIC for x in xs]
        lineERA = [slopeERA*x + interceptERA for x in xs]
        pl.plot(xs,linePIC,'g--',linewidth=1,alpha=0.35)
        pl.plot(xs,lineERA,'k--',linewidth=1,alpha=0.35)

        #years = range(1980,2016-N+2,1)
        #for label, x, y in zip([years[i] for i in [0,-1]], [f1ERA[i] for i in [0,-1]], [f2ERA[i] for i in [0,-1]]):
        #        plt.annotate(   label,xy=(x, y),xytext=(0, 0),color='k',
        #                        textcoords='offset points', ha='center', va='bottom',alpha=0.65)

        pl.grid(color='b')
        #pl.xlim(-475,475)
        #pl.ylim(-2.75,5.25)
        pl.xlabel('Sea level pressure [Pa]')
        pl.ylabel('Surface temperature [K]')
        pl.savefig('/mnt/climstorage/cian/scripts/figs/trends/bivar.clims.%syears.pdf' % (N),format='pdf')

	pl.show()
	"""

	"""
	Exp = str(sys.argv[1])
	N   = int(sys.argv[2])
        #LatRange1,LonRange1 = [67,80],[-50,-30]
        #LatRange2,LonRange2 = [65,80],[-80,-50]
	#LatRange1,LonRange1 = [65,80],[60,140]
	LatRange1,LonRange1 = [55,75],[30,100]
	LatRange2,LonRange2 = [75,83],[20,80]
	Field2              = ['pw','prw']
	if Exp == 'rcp85':
		Models_no = ['EC-Earth','MIROC4h','MPI-ESM-P']
		Models    = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*') if x[28:] not in Models_no]
		#year0s   = range(2005,2100-N+2,1)
		year0s    = range(1980,2016-N+2,1)
	elif Exp == 'historical':
		Models_no = ['MPI-ESM-P','CMCC-CM','EC-Earth','MIROC4h']
		Models    = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*') if x[28:] not in Models_no]
		year0s    = range(1850,2005-N+2,1)
        elif Exp == 'piControl':
                Models_no = ['GFDL-ESM2G','EC-EARTH','GFDL-CM3','GFDL-ESM2M','CMCC-CESM']
                Models    = [x[32:] for x in glob.glob('/mnt/climstorage/cian/piControl/*') if x[32:] not in Models_no]
                year0s    = range(1,75-N+2,1)
	elif Exp == 'ERAInt':
		Models = ['ERAInt']
		year0s = range(1980,2016-N+2,1)
	else:
		Exp    = 'ERA20c'
		Models = [Exp]
		year0s = range(1900,2010-N+2,1)
	M1,M2,F1,F2 = [],[],[],[]
	print Models
	fname_m = '/mnt/climstorage/cian/scripts/trends/histo/%s-%sN_%s-%sE_%s-%sN_%s-%sE/model_trends_%s_%syears_%s-%s.p' % \
			(LatRange1[0],LatRange1[1],LonRange1[0],LonRange1[1],LatRange2[0],LatRange2[1],LonRange2[0],LonRange2[1],Exp,N,year0s[0],year0s[-1])
        fname_f = '/mnt/climstorage/cian/scripts/trends/histo/%s-%sN_%s-%sE_%s-%sN_%s-%sE/model_clims_%s_%syears_%s-%s.p' % \
                        (LatRange1[0],LatRange1[1],LonRange1[0],LonRange1[1],LatRange2[0],LatRange2[1],LonRange2[0],LonRange2[1],Exp,N,year0s[0],year0s[-1])
	os.system('mkdir -p /mnt/climstorage/cian/scripts/trends/histo/%s-%sN_%s-%sE_%s-%sN_%s-%sE' % (LatRange1[0],LatRange1[1],LonRange1[0],LonRange1[1],LatRange2[0],LatRange2[1],LonRange2[0],LonRange2[1]))
	if (not os.path.isfile(fname_m)) and (not os.path.isfile(fname_f)):
		print 'Making file %s ...' % (fname_m)
		print 'Making file %s ...' % (fname_f)
		for i in range(len(Models)):
			Model = Models[i]
			M1.append([])
			M2.append([])
			F1.append([])
			F2.append([])
			fc = FieldChange(Source=Model,ExpType=Exp,Field2=Field2)
			for year0 in year0s:
				if Exp == 'piControl': year0 = year0 + fc.d1.Yearmin
				try:
					print 'base year is %s (%s) ...' % (year0,Model)
					m1,m2,f1,f2 = fc.modelTrends(year0,N,LatRange1,LonRange1,LatRange2,LonRange2,plot=False)
					M1[i].append(m1)
					M2[i].append(m2)
                                        F1[i].append(f1)
                                        F2[i].append(f2)
				except:
					print 'skipped base year %s (%s) ...' % (year0,Model)
		toPick([Models,M1,M2],fname_m)
		toPick([Models,F1,F2],fname_f)
	else:
		print 'File %s already exists ... did nothing.' % (fname_m)
		print 'File %s already exists ... did nothing.' % (fname_f)
	"""
