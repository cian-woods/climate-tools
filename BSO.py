import scipy.io
import matplotlib.pyplot as pl
import numpy as np
import sys,os

from drawCoastlinesNoRivers import *
from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from scipy import interpolate,stats
from scipy import ndimage as nd
from toPick import *

class BSO:

	def __init__(	self	):
	
		# Attributes
		self.mat  = scipy.io.loadmat('/mnt/climstorage/cian/BSO/obs/fluks_varmetransport1997_2014.mat')
		self.keys = self.mat.keys()
		print self.keys
		# Data structures
		self.structs = {}
		for key in self.keys:
			self.structs[key] = self.mat[key]
		# DataServers
		self.ds   = reDataServer(Field='T2',LevType='surface_analysis',Source='ERAInt',Year0=None)
		self.cids = MMDataServer(Field='ci')
		self.pwds = MMDataServer(Field='pw')
		self.spds = MMDataServer(Field='slp')
		self.shds = MMDataServer(Field='sshf')
		self.lwds = MMDataServer(Field='fls')
		self.lwdds = MMDataServer(Field='flds')
		self.swds = MMDataServer(Field='fsns')
		self.lhds = MMDataServer(Field='slhf')
		#self.vsds = MMDataServer(Field='V10')
		#self.usds = MMDataServer(Field='U10')
                self.vsds = MMDataServer(Field='V',LevRange=(1000,1000))
                self.usds = MMDataServer(Field='U',LevRange=(1000,1000))
		# Datestart and end
		self.datestart = (1997,8,22,0)
		self.dateend   = (2014,5,28,0)
		print 'Data from %s to %s' % (str(self.datestart),str(self.dateend))
		self.hourlist  = np.arange(self.ds.getHours(*self.datestart),self.ds.getHours(*self.dateend)+24,24)
		self.datelist  = [self.ds.getDate(hour) for hour in self.hourlist]
		# Data
		datelist,self.hfTOT30d   = self.extendAndMask(self.data('heatfluxTOT30d')*1e12, (1997,1,1,0), (2014,12,31,0))
		datelist,self.hfAW30d    = self.extendAndMask(self.data('heatfluxAW30d')*1e12, (1997,7,1,0), (2014,12,31,0))#(2014,6,30,0))
                datelist,self.hfTOT3mnd  = self.extendAndMask(self.data('heatfluxTOT3mnd')*1e12, (1997,7,1,0),(2014,12,31,0)) #(2014,6,30,0))
                datelist,self.hfAW3mnd   = self.extendAndMask(self.data('heatfluxAW3mnd')*1e12, (1997,7,1,0),(2014,12,31,0)) #(2014,6,30,0))
                datelist,self.hfTOT12mnd = self.extendAndMask(self.data('heatfluxTOT12mnd')*1e12, (1997,7,1,0), (2014,12,31,0))#(2014,6,30,0))
                datelist,self.hfAW12mnd  = self.extendAndMask(self.data('heatfluxAW12mnd')*1e12, (1997,7,1,0), (2014,12,31,0))#(2014,6,30,0))
		datelist,self.fTOT30d    = self.extendAndMask(self.data('fluxTOT30d'), (1997,1,1,0), (2014,12,31,0))
		datelist,self.fAW30d     = self.extendAndMask(self.data('fluxAW30d'), (1997,1,1,0), (2014,12,31,0))
		self.datelist            = datelist
		print self.datelist[0],self.datelist[-1]
		# Lambert
		self.proj = LambertProjector(boundinglat=45,resolution=150.)

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

	def regress(self,YearRange,Season,Field,lag=0):
		# Data
		years  = np.arange(YearRange[0],YearRange[1]+1,1)
		CI     = self.proj(np.array([self.cids.getSeason(Year=year+lag,Season=Season) for year in years]),self.cids.lon,self.cids.lat)
		SP     = self.proj(np.array([self.spds.getSeason(Year=year+lag,Season=Season) for year in years]),self.spds.lon,self.spds.lat)
		PW     = self.proj(np.array([self.pwds.getSeason(Year=year+lag,Season=Season) for year in years]),self.pwds.lon,self.pwds.lat)
		SH     = self.proj(np.array([self.shds.getSeason(Year=year+lag,Season=Season) for year in years]),self.shds.lon,self.shds.lat)
		LH     = self.proj(np.array([self.lhds.getSeason(Year=year+lag,Season=Season) for year in years]),self.lhds.lon,self.lhds.lat)
		LW     = self.proj(np.array([self.lwds.getSeason(Year=year+lag,Season=Season) for year in years]),self.lwds.lon,self.lwds.lat)
		LWD    = self.proj(np.array([self.lwdds.getSeason(Year=year+lag,Season=Season) for year in years]),self.lwdds.lon,self.lwdds.lat)
		SW     = self.proj(np.array([self.swds.getSeason(Year=year+lag,Season=Season) for year in years]),self.swds.lon,self.swds.lat)
		VS     = self.proj(np.array([self.vsds.getSeason(Year=year+lag,Season=Season).squeeze() for year in years]),self.vsds.lon,self.vsds.lat)
		US     = self.proj(np.array([self.usds.getSeason(Year=year+lag,Season=Season).squeeze() for year in years]),self.usds.lon,self.usds.lat)
		UVS    = np.sqrt(US**2 + VS**2)
		HF     = SH+LH+LW+SW
		HT     = np.array([self.getSeason(Year=year,Season=Season,Field=Field) for year in years])/1e12
		# Detrend
		ci ,m_ci ,p_ci  = self.detrend2d(years, CI)
		sp ,m_sp ,p_sp  = self.detrend2d(years, SP)
		pw ,m_pw ,p_pw  = self.detrend2d(years, PW)
		hf ,m_hf ,p_hf  = self.detrend2d(years, HF)
		us ,m_us ,p_us  = self.detrend2d(years, US)
		vs ,m_vs ,p_vs  = self.detrend2d(years, VS)
		uvs,m_uvs,p_uvs = self.detrend2d(years,UVS)
		sh ,m_sh ,p_sh  = self.detrend2d(years, SH)
		lh ,m_lh ,p_lh  = self.detrend2d(years, LH)
		lw ,m_lw ,p_lw  = self.detrend2d(years, LW)
		mht,cht         = np.polyfit(years,HT,1)
		ht              = HT - np.array([mht*x + cht for x in years])

		# Regress
		mci ,pci  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		msp ,psp  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mpw ,ppw  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mhf ,phf  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mvs ,pvs  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mus ,pus  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		muvs,puvs = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		msh ,psh  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mlh ,plh  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		mlw ,plw  = np.zeros((self.proj.nx,self.proj.ny)),np.zeros((self.proj.nx,self.proj.ny))
		for i in range(self.proj.nx):
			for j in range(self.proj.ny):
				slope, intercept, r_value, p_value, std_err = stats.linregress(ht,ci[:,i,j])
				mci[i,j],pci[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,sp[:,i,j])
                                msp[i,j],psp[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,pw[:,i,j])
                                mpw[i,j],ppw[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,hf[:,i,j])
                                mhf[i,j],phf[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,us[:,i,j])
                                mus[i,j],pus[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,vs[:,i,j])
                                mvs[i,j],pvs[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,uvs[:,i,j])
                                muvs[i,j],puvs[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,sh[:,i,j])
                                msh[i,j],psh[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,lh[:,i,j])
                                mlh[i,j],plh[i,j] = slope,p_value/2.
                                slope, intercept, r_value, p_value, std_err = stats.linregress(ht,lw[:,i,j])
                                mlw[i,j],plw[i,j] = slope,p_value/2.
				
		# Plot
		cseq  = 13
		pl.figure(1)
		sx,sy = np.ma.masked_where(pci>0.05,self.proj.x),np.ma.masked_where(pci>0.05,self.proj.y)
		cf    = pl.contourf(self.proj.x,self.proj.y,mci,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cbar  = pl.colorbar(cf)
		cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (self.cids.units))
		pl.plot(sx,sy,'k.',alpha=0.5)
		drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		self.proj.m.drawparallels([70,80],latmax=90)
		pl.title(self.cids.long_name)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('ci',YearRange[0],YearRange[1],Season),format='pdf')

		pl.figure(2)
                sx,sy = np.ma.masked_where(psp>0.05,self.proj.x),np.ma.masked_where(psp>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,msp,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (self.spds.units))
                pl.plot(sx,sy,'k.',alpha=0.5)
		urot,vrot = self.proj.m.rotate_vector(mus,mvs,self.proj.lon,self.proj.lat,returnxy=False)
		rotmask   = (pus>0.05)&(pvs>0.05)
		urot,vrot = np.ma.masked_array(urot,mask=rotmask),np.ma.masked_array(vrot,mask=rotmask)
		pl.quiver(self.proj.x[::1,::1],self.proj.y[::1,::1],urot[::1,::1],vrot[::1,::1],pivot='tail',alpha=0.85)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title(self.spds.long_name)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('slp',YearRange[0],YearRange[1],Season),format='pdf')

		pl.figure(3)
                sx,sy = np.ma.masked_where(ppw>0.05,self.proj.x),np.ma.masked_where(ppw>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mpw,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (self.pwds.units))
                pl.plot(sx,sy,'k.',alpha=0.5)
		pl.quiver(self.proj.x,self.proj.y,urot,vrot,pivot='tail',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title(self.pwds.long_name)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('pw',YearRange[0],YearRange[1],Season),format='pdf')

		pl.figure(4)
                sx,sy = np.ma.masked_where(phf>0.05,self.proj.x),np.ma.masked_where(phf>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mhf,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [W m$^{-2}$ TW$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title('Total heat flux')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('hf',YearRange[0],YearRange[1],Season),format='pdf')

                pl.figure(5)
		cseq0 = np.arange(-0.6,0.6+0.1,0.1)
                sx,sy = np.ma.masked_where(psh>0.05,self.proj.x),np.ma.masked_where(psh>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,msh,cseq0,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [W m$^{-2}$ TW$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title('Sensible heat flux')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('sh',YearRange[0],YearRange[1],Season),format='pdf')

                pl.figure(6)
                sx,sy = np.ma.masked_where(plh>0.05,self.proj.x),np.ma.masked_where(plh>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mlh,cseq0,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [W m$^{-2}$ TW$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title('Latent heat flux')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('lh',YearRange[0],YearRange[1],Season),format='pdf')

                pl.figure(7)
                sx,sy = np.ma.masked_where(plw>0.05,self.proj.x),np.ma.masked_where(plw>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mlw,cseq0,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [W m$^{-2}$ TW$^{-1}$]')
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title('Longwave flux')
                pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('lw',YearRange[0],YearRange[1],Season),format='pdf')

		pl.figure(8)
                sx,sy = np.ma.masked_where(pus>0.05,self.proj.x),np.ma.masked_where(pus>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mus,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (self.usds.units))
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title(self.usds.long_name)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('US',YearRange[0],YearRange[1],Season),format='pdf')

                pl.figure(9)
                sx,sy = np.ma.masked_where(pvs>0.05,self.proj.x),np.ma.masked_where(pvs>0.05,self.proj.y)
                cf    = pl.contourf(self.proj.x,self.proj.y,mvs,cseq,cmap=pl.cm.RdBu_r,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (self.vsds.units))
                pl.plot(sx,sy,'k.',alpha=0.5)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                self.proj.m.drawparallels([70,80],latmax=90)
                pl.title(self.vsds.long_name)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % ('VS',YearRange[0],YearRange[1],Season),format='pdf')

		pl.show()

	def extendAndMask(self,data,datestart,dateend):
		hourlist = np.arange(self.ds.getHours(*datestart),self.ds.getHours(*dateend)+24,24)
		datelist = [self.ds.getDate(hour) for hour in hourlist]
		mask     = [date not in self.datelist for date in datelist]
		ix0,ix1  = datelist.index(self.datestart),datelist.index(self.dateend)
		data     = np.append(np.zeros(ix0),data)
		data     = np.append(data,np.zeros(len(mask)-ix1-1))
		data     = np.ma.masked_array(data,mask=mask)
		return datelist,data

	def data(self,key):
		f = np.array(self.structs[key]).squeeze()
		return f

	def moving_average(self,x,y,n=3):
		if n%2 == 0: return 0,0
        	ret = np.cumsum(y, dtype=float)
        	ret[n:] = ret[n:] - ret[:-n]
        	return x[n/2:-n/2 + 1],ret[n - 1:]/n

	def convertDatesToYear(self,datelist):
		hoursinyear             = 365*24
		year0,month0,day0,hour0 = datelist[0]
		year1,month1,day1,hour1 = datelist[-1]
		h0,h1 = self.ds.getHours(year0,1,1,0),self.ds.getHours(year1,1,1,0)
		hs,he = self.ds.getHours(*datelist[0]),self.ds.getHours(*datelist[-1])
		years = year0 + 1.*(hs-h0)/hoursinyear
		yeare = year1 + 1.*(he-h1)/hoursinyear
		return np.linspace(years,yeare,len(datelist))

	def getTimeSlice(self,date0=None,date1=None,Field='hfTOT30d'):
		if (date0 == None) and (date1 == None): date0,date1 = self.datestart,self.dateend
		xs0,xs1  = self.datelist.index(date0),self.datelist.index(date1)
		yearlist = self.convertDatesToYear(self.datelist[xs0:xs1+1])
		if Field ==   'hfTOT30d' : return yearlist,self.hfTOT30d[xs0:xs1+1]
		if Field ==    'hfAW30d' : return yearlist,self.hfAW30d[xs0:xs1+1]
                if Field == 'hfTOT12mnd' : return yearlist,self.hfTOT12mnd[xs0:xs1+1]
                if Field ==  'hfAW3mnd' : return yearlist,self.hfAW3mnd[xs0:xs1+1]
		if Field == 'hfTOT3mnd' : return yearlist,self.hfTOT3mnd[xs0:xs1+1]
		if Field ==  'hfAW12mnd' : return yearlist,self.hfAW12mnd[xs0:xs1+1]
		if Field ==    'fTOT30d' : return yearlist,self.fTOT30d[xs0:xs1+1]

	def getSeason(self,Year,Season,Field='hfTOT30d'):
		if Season == 'Annual':
                        date0 = (Year,  1,  1, 0)
                        date1 = (Year, 12, 31, 0)
		if Season == 'Annual_wc':	# Winter centered annual mean (July-June)
                        date0 = (Year-1, 7,  1, 0)
                        date1 = (Year  , 6, 30, 0)
                if Season == 'Annual_MF':       # Feb ending annual mean (March-Feb)
                        date0 = (Year-1, 3,  1, 0)
                        date1 = (Year  , 2, 28, 0)
		if Season == 'MAMJJASO_l1':
			date0 = (Year-1 ,  3,  1, 0)
			date1 = (Year-1 , 10, 31, 0)
		if Season == 'DJF':
			date0 = (Year-1, 12,  1, 0)
			date1 = (Year  ,  2, 28, 0)
                if Season == 'NDJFM':
                        date0 = (Year-1, 11,  1, 0)
                        date1 = (Year  ,  3, 31, 0)
                if Season == 'ONDJFM':
                        date0 = (Year-1, 10,  1, 0)
                        date1 = (Year  ,  3, 31, 0)
                if Season == 'NDJ':
                        date0 = (Year-1, 11,  1, 0)
                        date1 = (Year  ,  1, 31, 0)
                if Season == 'NDJF':
                        date0 = (Year-1, 11,  1, 0)
                        date1 = (Year  ,  2, 28, 0)
                if Season == 'ONDJF':
                        date0 = (Year-1, 10,  1, 0)
                        date1 = (Year  ,  2, 28, 0)
		if Season == 'SON':
                        date0 = (Year,  9,  1, 0)
                        date1 = (Year, 11, 30, 0)
                if Season == 'OND':
                        date0 = (Year, 10,  1, 0)
                        date1 = (Year, 12, 31, 0)
                if Season == 'ONDJ':
                        date0 = (Year-1, 10,  1, 0)
                        date1 = (Year  ,  1, 31, 0)
                if Season == 'JFM':
                        date0 = (Year, 1,  1, 0)
                        date1 = (Year, 3, 31, 0)
                if Season == 'JJA':
                        date0 = (Year,  6,  1, 0)
                        date1 = (Year,  8, 31, 0)
                if Season == 'MAM':
                        date0 = (Year,  3,  1, 0)
                        date1 = (Year,  5, 31, 0)
                if Season == 'DJ':
                        date0 = (Year-1, 12,  1, 0)
                        date1 = (Year,    1, 31, 0)
                if Season == 'J':
                        date0 = (Year, 1,  1, 0)
                        date1 = (Year, 1, 31, 0)
                if Season == 'D':
                        date0 = (Year, 12,  1, 0)
                        date1 = (Year, 12, 31, 0)
		xs0,xs1 = self.datelist.index(date0),self.datelist.index(date1)
		if Field == 'hfTOT30d':
			data = self.hfTOT30d[xs0:xs1+1].mean(axis=0)
		if Field == 'hfAW30d':
			data = self.hfAW30d[xs0:xs1+1].mean(axis=0)
		if Field == 'hfTOT12mnd':
			data = self.hfTOT12mnd[xs0:xs1+1].mean(axis=0)
		if Field == 'hfAW12mnd':
			data = self.hfAW12mnd[xs0:xs1+1].mean(axis=0)
                if Field == 'hfTOT3mnd':
			data = self.hfTOT3mnd[xs0:xs1+1].mean(axis=0)
                if Field == 'hfAW3mnd':
			data = self.hfAW3mnd[xs0:xs1+1].mean(axis=0)
		if Field == 'fTOT30d':
			data = self.fTOT30d[xs0:xs1+1].mean(axis=0)
		if Field == 'fAW30d':
			data = self.fAW30d[xs0:xs1+1].mean(axis=0)
		print Year,data.data
		if data.data == 0: data = -999
		return data

	def standardize(self,data):
		return (data-data.mean())/data.std()

if __name__ == "__main__":

	b = BSO()
 
	YearRange = (int(sys.argv[1]),int(sys.argv[2]))
	Season    = str(sys.argv[3])
	Field1    = str(sys.argv[4])
	Field2    = str(sys.argv[5])

	#b.regress(YearRange,Season,Field1)
	#sys.exit()

	#yr,ht_d = b.getTimeSlice(date0=(1997,1,1,0),date1=(2014,12,31,0),Field=Field1)
	#yr0,ht0 = b.moving_average(yr,ht_d,365)
	#sis     = len(b.ds.getDateList(1997,Season=Season))*(6*60*60)
	years   = range(YearRange[0],YearRange[1]+1,1)
	F1      = np.array([b.getSeason(year,Season=Season,Field=Field1) for year in years])
	F2      = np.array([b.getSeason(year,Season=Season,Field=Field2) for year in years])
	F1,F2   = np.ma.masked_where(F1==-999,F1),np.ma.masked_where(F2==-999,F2)
	F1mean  = F1.mean()
	F2mean  = F2.mean()

	slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(years,F1)
	line1 = np.array([slope1*x + intercept1 for x in years])
	f1    = F1 - line1
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(years,F2)
        line2 = np.array([slope2*x + intercept2 for x in years])
        f2    = F2 - line2

	# Save timeseries
	fname1 = '/mnt/climstorage/cian/BSO/%s.bso.full.%s-%s.%s.p'   % (Field1,years[0],years[-1],Season)
	fname2 = '/mnt/climstorage/cian/BSO/%s.bso.dtrend.%s-%s.%s.p' % (Field1,years[0],years[-1],Season)
	fname3 = '/mnt/climstorage/cian/BSO/%s.bso.full.%s-%s.%s.p'   % (Field2,years[0],years[-1],Season)
	fname4 = '/mnt/climstorage/cian/BSO/%s.bso.dtrend.%s-%s.%s.p' % (Field2,years[0],years[-1],Season)

	if not os.path.isfile(fname1):
		print fname1
	        toPick([years,F1],fname1)
	else:
		os.system('rm %s' % (fname1))
		toPick([years,F1],fname1)

	if not os.path.isfile(fname2):
		print fname2
	        toPick([years,f1],fname2)
        else:
                os.system('rm %s' % (fname2))
                toPick([years,f1],fname2)

        if not os.path.isfile(fname3):
                print fname3
                toPick([years,F2],fname3)
        else:
                os.system('rm %s' % (fname3))
                toPick([years,F2],fname3)

        if not os.path.isfile(fname4):
                print fname4
                toPick([years,f2],fname4)
        else:
                os.system('rm %s' % (fname4))
                toPick([years,f2],fname4)

	pl.figure(1)
	pl.plot(years,b.standardize(F1),'r')
	pl.plot(years,b.standardize(F2),'b')

        pl.figure(2)
        pl.plot(years,b.standardize(f1),'r')
        pl.plot(years,b.standardize(f2),'b')

	pl.show()

	sys.exit()

	print slope1*sis/1e18
	
	fig,ax1 = pl.subplots(num=1)
	ax2     = ax1.twinx()
	ax1.plot(years,F1/1e12,'k.-',linewidth=0.6,alpha=0.5,label='HT$_{BSO}$')
	ax2.plot(years,F2,'b.-',linewidth=1,alpha=0.35,label='VT$_{BSO}$')
	ax1.set_ylabel('Heat transport [TW]')
	ax2.set_ylabel('Volume transport [Sv]')
	ax1.set_xlabel('Year')
	ax1.grid()
	pl.title('%s; trend = %s TW year$^{-1}$; p = %s; r = %s' % (Season,round(slope1/1e12,3),round(p_value1/2.,4),round(np.corrcoef(f1,f2)[0][1],3)))
	ax1.set_xlim(years[0],years[-1])
	ax1.legend(loc=2)
	ax2.legend(loc=1)
	pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_moor.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

	pl.figure(2)
	slope, intercept, r_value, p_value, std_err = stats.linregress(yr,ht_d)
	line = np.array([slope*x + intercept for x in yr])
	pl.plot(yr,ht_d/1e12,'k-',linewidth=0.5,alpha=0.6)
	pl.plot(yr0,ht0/1e12,'r-',linewidth=1.0,alpha=0.7)
	#pl.plot(yr,line/1e12,'k--',linewidth=0.7,alpha=0.5)
	pl.ylabel('Heat transport [TW]')
	pl.xlabel('Year')
	pl.xticks(range(1998,2014+2,2))
	pl.xlim(yr[0],yr[-1])
	pl.grid()
	pl.title('trend = %s TW year$^{-1}$; p = %s' % (round(slope/1e12,3),round(p_value/2.,4)))
	pl.show()

