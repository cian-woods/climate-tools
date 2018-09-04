import sys
sys.path.insert(0, '/home/cian/scripts')
from UnPickle import *
from toPick import *
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from Trajectories import *
from scipy import stats,interpolate
from drawCoastlinesNoRivers import drawCoastlinesNoRivers

import matplotlib.pyplot as pl
import numpy as np
import scipy.io
import random

class Intrusions:

	def __init__(	self,
			OpenTraj = True,
			Source   = 'ERAInt'	):

		self.Source = Source
		# Trajectories class
		self.blat = 55
		self.res  = 350
		self.tr   = Trajectories(Source=self.Source,dT=0.25,ff=1,blat=self.blat,res=self.res)
		# Intrusions; ERAInt is in 4 files (1979 JFM, 1980-2015 NDJFM, 2016 NDJFM, 2016 ND)
		G1,Q1,D1 = unpick('/mnt/climstorage/cian/intrusions/%s_intrusions.1979-1979.FM.6x6hr.9deg.200.6dt.20.5.filtered.p'   % (self.Source))
		G2,Q2,D2 = unpick('/mnt/climstorage/cian/intrusions/%s_intrusions.1980-2015.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p' % (self.Source))
		G3,Q3,D3 = unpick('/mnt/climstorage/cian/intrusions/%s_intrusions.2016-2016.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p' % (self.Source))
		G4,Q4,D4 = unpick('/mnt/climstorage/cian/intrusions/%s_intrusions.2016-2016.ND.6x6hr.9deg.200.6dt.20.5.filtered.p'    % (self.Source))
		self.G,self.Q,self.D = G1+G2+G3+G4,Q1+Q2+Q3+Q4,D1+D2+D3+D4
		print len(self.G)
		self.maxdates,self.maxlons = self.datesAtMaxInten(self.G,self.Q,self.D)
		# Forward and backward trajectories
		if OpenTraj == True:
			print 'Opening forward trajectories for NDJFM 1979 to 2016 ...'
			LONf1,LATf1,Pf1 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_fwrd_tracks.1979-1979.FM.6hour.filtered.p' % (self.Source))
			print len(LONf1[0][0][0])
			LONf2,LATf2,Pf2 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_fwrd_tracks.1980-2015.NDJFM.6hour.filtered.p' % (self.Source))
			print len(LONf2[0][0][0])
			LONf3,LATf3,Pf3 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_fwrd_tracks.2016-2016.NDJFM.6hour.filtered.p' % (self.Source))
			print len(LONf3[0][0][0])
			LONf4,LATf4,Pf4 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_fwrd_tracks.2016-2016.ND.6hour.filtered.p' % (self.Source))
			print len(LONf4[0][0][0])
			LONf,LATf,Pf    = LONf1+LONf2+LONf3+LONf4,LATf1+LATf2+LATf3+LATf4,Pf1+Pf2+Pf3+Pf4
			print 'Opening backward trajectories for NDJFM 1979 to 2016 ...'
                	LONb1,LATb1,Pb1 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_back_tracks.1979-1979.FM.6hour.filtered.p' % (self.Source))
			print len(LONb1[0][0][0])
                	LONb2,LATb2,Pb2 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_back_tracks.1980-2015.NDJFM.6hour.filtered.p' % (self.Source))
			print len(LONb2[0][0][0])
			LONb3,LATb3,Pb3 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_back_tracks.2016-2016.NDJFM.6hour.filtered.p' % (self.Source))
			print len(LONb3[0][0][0])
			LONb4,LATb4,Pb4 = unpick('/mnt/climstorage/cian/scripts/intrusiontracks/%s_back_tracks.2016-2016.ND.6hour.filtered.p' % (self.Source))
			print len(LONb4[0][0][0])
			LONb,LATb,Pb    = LONb1+LONb2+LONb3+LONb4,LATb1+LATb2+LATb3+LATb4,Pb1+Pb2+Pb3+Pb4
			# Stitch intrusion trajectories together
			for i in range(len(LONf)):
				for j in range(len(LONf[i])):
                			LONf[i][j] = np.append(LONb[i][j][:,::-1][:,20:],LONf[i][j][:,1:20+1],axis=1)
                			LATf[i][j] = np.append(LATb[i][j][:,::-1][:,20:],LATf[i][j][:,1:20+1],axis=1)
                			Pf[i][j]   = np.append(np.array(Pb[i][j])[:,::-1][:,0:],np.array(Pf[i][j])[:,1:40+1],axis=1)
			self.LON,self.LAT,self.P = LONf,LATf,Pf
			self.nsteps = len(self.LON[0][0][0])
		dur = np.array([len(self.D[ii]) for ii in range(len(self.D))]).mean()
		print 'Mean duration of injection events is %s days' % (dur/4.)
		# DataServers
		self.Field = 'T2'
		self.reds  = reDataServer(Field=self.Field,Source=self.Source,LevType='surface_analysis',LatRange=(-15,15),LonRange=(80,160))
		self.cmds  = cmipDataServer(Field='psl',LevType='surface',Source='CanESM2',DataFreq='day',ExpType='rcp85')

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

	def makeDataset1(self,case):
		dates,datelist = self.topEvents(50,case=case)
		xs,ls          = self.associateEvents(dates)
#		G,Q,D,LON,LAT,P          = self.extractIntrusions(xs)
		File = open('AssociatedIntrusions_%s.txt' % (case),'w')
		File.write('Rank|Intrusion|Timestep|Date|LonStart|LonEnd\n')
		for i in range(50):
			date = str(dates[i])
			if ls[i] == []:
				File.write('%s|%s|%s|%s|%s|%s\n' % (i+1,-999,-999,-999,-999,-999))
				continue
			for j in range(len(ls[i])):
				x = ls[i][j]
				for k in range(len(self.D[x])):
					date = str(self.D[x][k])
					File.write('%s|%s|%s|%s|%s|%s\n' % (i+1,j+1,k,date,self.G[x][k][0],self.G[x][k][-1]))
		File.close()

        def makeDataset2(self,case):
                dates,datelist = self.topEvents(50,case=case)
                xs,ls          = self.associateEvents(dates)
#               G,Q,D,LON,LAT,P          = self.extractIntrusions(xs)
                File = open('SectorIntrusions_%s.txt' % (case),'w')
                File.write('Rank|Date|ALL|ATL|PAC\n')
                for i in range(50):
			inds = ls[i]
			lons = [self.maxlons[j] for j in inds]
                        date = str(dates[i])
			Ntot = len(inds)
			NAtl = len(self.filterLons(lons,Sector=(290,110)))
			NPac = len(self.filterLons(lons,Sector=(110,290)))
			File.write('%s|%s|%s|%s|%s\n' % (i+1,date,Ntot,NAtl,NPac))
                File.close()

	def datesAtMaxInten(self,G,Q,D):
		dates,lons = [],[]
		for i in range(len(Q)):
			ymax    = [max(Q[i][j]) for j in range(len(Q[i]))]
			yargmax = [np.argmax(Q[i][j]) for j in range(len(Q[i]))]
			k       = np.argmax(ymax)
			l       = yargmax[k]
			dates.append(D[i][k])
			lons.append(G[i][k][l])
		return dates,lons

	def intrusionVar(self,years,mio=(11,3)):
		INT = {}
		for year in years:
			INT[year] = []
		for ii in range(len(self.D)):
			date = self.D[ii][0]
			year,month,day,hour = date
			if (month >=mio[0]) or (month <= mio[1]):
				if month >= 10: year = year + 1
				try:
					INT[year].append(self.Q[ii])
				except:
					pass
		# Mean number, flux and intensity per event
		N = np.array([len(INT[key]) for key in INT.keys()])
		E = np.array([sum([sum([sum(q1) for q1 in q]) for q in INT[year]])/4. for year in INT.keys()])
		I = E/N
		# Total duration
		dur   = np.array([sum([len(q) for q in INT[year]]) for year in INT.keys()])

		durtrend, durc, r_value, p_value, std_err = stats.linregress(years,dur)
		#durtrend,durc = np.polyfit(years,dur,1)
		linedur       = [durtrend*year + durc for year in years]
		pl.figure(1)
		pl.plot(years,dur,'k')
		pl.plot(years,linedur,'k--')
		pl.xlabel('Year')
		pl.ylabel('Duration')
		pl.title('trend = %s decade$^{-1}$; p = %s' % (round(durtrend*10,2),round(p_value/2,3)))
		pl.grid()
		pl.xlim(1980,2016)
		pl.savefig('figs/durtrend.pdf',format='pdf')

		Ntrend, Nc, r_value, p_value, std_err = stats.linregress(years,N)
                #Ntrend,Nc = np.polyfit(years,N,1)
                lineN     = [Ntrend*year + Nc for year in years]
		pl.figure(2)
                pl.plot(years,N,'k')
		pl.plot(years,lineN,'k--')
                pl.xlabel('Year')
                pl.ylabel('Number of injection events')
		pl.title('trend = %s decade$^{-1}$; p = %s' % (round(Ntrend*10,2),round(p_value/2,3)))
                pl.grid()
		pl.xlim(1980,2016)
		pl.savefig('figs/Ntrend.pdf',format='pdf')
		pl.show()

		inten,x = [],0
		width   = []
		for year in INT.keys():
			inten.append([])
			width.append([])
			for q in INT[year]:
				for q1 in q:
					width[x].append(len(q1))
					for q2 in q1:
						inten[x].append(q2)
			x = x + 1
		inten = np.array([np.mean(i) for i in inten])
		width = np.array([np.mean(i) for i in width])
                # Terms
                durc   = dur.mean()
                duri   = dur - dur.mean()
                intenc = inten.mean()
                inteni = inten - intenc
                widthc = width.mean()
                widthi = width-widthc
                t0 =  (E-E.mean())/E.mean()
                t1 =  (inteni/intenc)
                t2 =  (duri/durc)
                t3 =  (widthi/widthc)
                t4 =  (duri*inteni)/(durc*intenc)
                t5 =  (inteni*widthi)/(intenc*widthc)
                t6 =  (duri*widthi)/(durc*widthc)
                t7 =  (duri*inteni*widthi)/(durc*intenc*widthc)
		return N,E,I,t1,t2,t3

	def seasonalVar(self,years,Season='NDJFM',blat=80.):
		# Seasonal mean field
		proj = LambertProjector(boundinglat=blat,resolution=80.)
                #reds = reDataServer(Field='T2',LevType='surface_analysis')
                #s    = np.array([reds.getDataSnaps(year,Season=Season).mean(axis=0) for year in years])
		mmds = MMDataServer(Field='T2')
		s    = np.array([mmds.getSeason(year,Season=Season) for year in years])
		s    = proj(s,mmds.lon,mmds.lat).mean(axis=-1).mean(axis=-1)
		return s

	def moistureSeries(self,years,mio=(11,3)):
		vq,d = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/70/%s/%s.moist.1979-2016.30-1000hPa.70N.Annual.p' % (self.Source,self.Source))
		vq   = self.interpolate(vq)
		F    = {}
                for year in years:
                        F[year] = []
		for ii in range(len(d)):
			date = d[ii]
			year,month,day,hour = date
			if (month >= mio[0]) or (month <= mio[1]):
				if month >= 10: year = year + 1
				try:
					flux = np.ma.masked_where(vq[ii]<0,vq[ii]).sum()
					if type(flux) != np.ma.core.MaskedConstant: F[year].append(flux)
				except:
					pass
		f = np.array([np.sum(F[key]) for key in F.keys()])
		return f

	def plotInjectionSeries(self,years,detrend=True,Slice=['NDJFM',(11,3)],blat=80.):
		# Time series
		N,E,I,i,d,w = self.intrusionVar(years,Slice[1])
		dw,di,iw = d + w, d + i, i + w
		F        = self.moistureSeries(years,Slice[1])
		R        = F - E
		s        = self.seasonalVar(years,Slice[0],blat=blat)
		# Detrend
		if detrend:
			N  = N - np.array([np.polyfit(years,N,1)[0]*x + np.polyfit(years,N,1)[1] for x in years])
			E  = E - np.array([np.polyfit(years,E,1)[0]*x + np.polyfit(years,E,1)[1] for x in years])
			I  = I - np.array([np.polyfit(years,I,1)[0]*x + np.polyfit(years,I,1)[1] for x in years])
			F  = F - np.array([np.polyfit(years,F,1)[0]*x + np.polyfit(years,F,1)[1] for x in years])
			R  = R - np.array([np.polyfit(years,R,1)[0]*x + np.polyfit(years,R,1)[1] for x in years])
			i  = i - np.array([np.polyfit(years,i,1)[0]*x + np.polyfit(years,i,1)[1] for x in years])
			d  = d - np.array([np.polyfit(years,d,1)[0]*x + np.polyfit(years,d,1)[1] for x in years])
			w  = w - np.array([np.polyfit(years,w,1)[0]*x + np.polyfit(years,w,1)[1] for x in years])
			dw = dw - np.array([np.polyfit(years,dw,1)[0]*x + np.polyfit(years,dw,1)[1] for x in years])
			di = di - np.array([np.polyfit(years,di,1)[0]*x + np.polyfit(years,di,1)[1] for x in years])
			iw = iw - np.array([np.polyfit(years,iw,1)[0]*x + np.polyfit(years,iw,1)[1] for x in years])
			s  = s - np.array([np.polyfit(years,s,1)[0]*x + np.polyfit(years,s,1)[1] for x in years])
		# Standard deviations
		sigN,sigE,sigI    = N.std(),E.std(),I.std()
		sigi,sigd,sigw    = i.std(),d.std(),w.std()
		sigdw,sigdi,sigiw = dw.std(),di.std(),iw.std()
		sigF,sigR         = F.std(),R.std()
		sigs              = s.std()
		# Standardise
		N,E,I    = (N-N.mean())/N.std(),(E-E.mean())/E.std(),(I-I.mean())/I.std()
		i,d,w    = (i-i.mean())/i.std(),(d-d.mean())/d.std(),(w-w.mean())/w.std()
		dw,di,iw = (dw-dw.mean())/dw.std(),(di-di.mean())/di.std(),(iw-iw.mean())/iw.std()
		F,R      = (F-F.mean())/F.std(),(R-R.mean())/R.std()
		s        = (s-s.mean())/s.std()
		# Plots
		self.plotSeries(years,N,E,sigN,sigE,'Number of injections','Injection transport','figs/N_E_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,E,I,sigE,sigI,'Injection transport','Injection intensity','figs/E_I_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,N,I,sigN,sigI,'Number of injections','Injection intensity','figs/N_I_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,E,F,sigE,sigF,'Injection transport','Total transport','figs/E_F_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,E,R,sigE,sigR,'Injection transport','Residual transport','figs/E_R_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,E,s,sigE,sigs,'Injection transport','%s %s %s' % (self.Field,blat,Slice[0]),'figs/E_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,F,s,sigF,sigs,'Total transport','%s %s %s' % (self.Field,blat,Slice[0]),'figs/F_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,R,s,sigR,sigs,'Residual transport','%s %s %s' % (self.Field,blat,Slice[0]),'figs/R_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,N,s,sigN,sigs,'Number of injections','%s %s %s' % (self.Field,blat,Slice[0]),'figs/N_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,I,s,sigI,sigs,'Injection intensity','%s %s %s' % (self.Field,blat,Slice[0]),'figs/I_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,i,s,sigi,sigs,'Mean intensity','%s %s %s' % (self.Field,blat,Slice[0]),'figs/i_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,d,s,sigd,sigs,'Duration','%s %s %s' % (self.Field,blat,Slice[0]),'figs/d_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,w,s,sigw,sigs,'Width','%s %s %s' % (self.Field,blat,Slice[0]),'figs/w_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,dw,s,sigdw,sigs,'Duration + Width','%s %s %s' % (self.Field,blat,Slice[0]),'figs/dw_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,di,s,sigdi,sigs,'Duration + Intensity','%s %s %s' % (self.Field,blat,Slice[0]),'figs/di_s_%s.pdf' % (Slice[0]),Slice[0])
		self.plotSeries(years,iw,s,sigiw,sigs,'Intensity + Width','%s %s %s' % (self.Field,blat,Slice[0]),'figs/iw_s_%s.pdf' % (Slice[0]),Slice[0])

	def plotSeries(self,x,y1,y2,sig1,sig2,label1,label2,savename,Season):
		r = np.corrcoef(y1,y2)[0][1]
		PP,MM = [],[]
                pp, = pl.plot(x,y1,'k',linewidth=1.5)
                PP.append(pp)
                MM.append(label1)
                pp, = pl.plot(x,y2,'r',linewidth=1.5)
                PP.append(pp)
                MM.append(label2)
                pl.xlim(x[0],x[-1])
                pl.ylabel('Standard deviation')
                pl.xlabel('Years')
                pl.title('%s, $\sigma1$ = %s, $\sigma2$ = %s, n = %s, r = %s' % (Season,round(sig1,4),round(sig2,4),len(x),round(r,4)))
                lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
                pl.setp(lg.get_title(),fontsize=10)
                pl.grid()
                pl.savefig('%s' % (savename),format='pdf')
                pl.close()

	def topEvents(self,n=50,case='warm'):
		if case == 'warm': var_rank,I = 'var_rankw','Iw'
		if case == 'cold': var_rank,I = 'var_rankc','Ic'
		# Warm Events (extract top n)
		years = range(1979,2017+1,1)
		mat = scipy.io.loadmat('/mnt/climstorage/gabriele/Data/Warm_Arctic/temp_Artext_nf_ERAiAll_y_y_7_NM.mat')
		Iw    = mat[I].squeeze()[:n]
		# Make dateList
		datelist = []
		for year in years:
			datelist = datelist + self.cmds.getDateList(year,Season='NDJFM')[1:]
		datelist = datelist[60:-91+1]	
		dates    = [datelist[iw-1] for iw in Iw]
		# Save dates
		fname = 'dates.%s.%s.%s.p' % (case,n,self.Source)
		if not os.path.isfile(fname): toPick(dates,fname)
		return dates,datelist	

	def associateEvents(self,dates,mind=-8,maxd=-2):
		# Associate intrusions to warm arctic events
		days0b = np.array([self.reds.getDays(*self.D[i][0]) for i in range(len(self.D))])	# Beginning days of injections
		days0e = np.array([self.reds.getDays(*self.D[i][-1]) for i in range(len(self.D))])	# Ending days of injections
		days   = np.array([self.reds.getDays(*date) for date in dates])				# Warm event days
		xs,ls  = [],[]
		for i in range(len(days)):
			day0  = days[i]							# Day of warm event
			dtEnd = days0e - day0						# Day difference between intrusion start dates and warm event
			dtBeg = days0b - day0						# Day difference between intrusion end dates and warm event
			xxs   = list(np.where(((dtBeg<mind)&(dtEnd>maxd))|((dtBeg>=mind)&(dtBeg<=maxd))|((dtEnd>=mind)&(dtEnd<=maxd))==True)[0])	
			xs   = xs + xxs
			ls.append(xxs)
		xs = list(set(xs))							# Exclude intrusions associated multiple times
		return xs,ls

	def lagAssociation(self,Sector=(0,360)):
		YearRange        = (1979,2000)
		years            = range(YearRange[0],YearRange[1]+1,1)
		dates1,datelist1 = self.topEvents(n=50,case='warm')
		dates2,datelist2 = self.topEvents(n=50,case='cold')
		dates1,dates2    = [date for date in dates1 if date[0] in years],[date for date in dates2 if date[0] in years]
		dates1n,dates2n  = len(dates1),len(dates2)	
		lags,r1n,r1d,r2n,r2d,c = np.arange(-20,10+1,1),[],[],[],[],2
		for i in lags:
			# Warm
			xs1,ls1  = self.associateEvents(dates1,mind=i-c,maxd=i+c)
			mxl1n    = [self.maxlons[xs0i] for xs0i in xs1]
			ll1n     = self.filterLons(mxl1n,Sector=Sector)
			xs1n     = len(ll1n)
			r1n.append(xs1n)
                        xs1d     = np.mean([ np.sum([len(self.D[xxsi]) for xxsi in xxs]) for xxs in ls1])
			r1d.append(xs1d)
			# Cold
                        xs2,ls2 = self.associateEvents(dates2,mind=i-c,maxd=i+c)
                        mxl2n    = [self.maxlons[xs0i] for xs0i in xs2]
                        ll2n     = self.filterLons(mxl2n,Sector=Sector)
                        xs2n     = len(ll2n)
                        r2n.append(xs2n)
                        xs2d     = np.mean([ np.sum([len(self.D[xxsi]) for xxsi in xxs]) for xxs in ls2])
                        r2d.append(xs2d)
		r1d,r2d = np.array(r1d)/4,np.array(r2d)/4
		r1n,r2n = 1.*np.array(r1n)/dates1n,1.*np.array(r2n)/dates2n
		pl.figure(1)
		pl.plot(lags,r1n,'r',linewidth=2,label='Warm extremes')
		pl.plot(lags,r2n,'b',linewidth=2,label='Cold extremes')
		#pl.plot([-20,10],[1.28,1.28],'k--')	# -8 to -2 1979-2016
		#pl.plot([-20,10],[0.96,0.96],'k--')
		#pl.plot([-20,10],[1.62,1.62],'k--')
                #pl.plot([-20,10],[1.20,1.20],'k--')	# -8 to -2 1979-2000
                #pl.plot([-20,10],[0.76,0.76],'k--')
                #pl.plot([-20,10],[1.63,1.63],'k--')
                #pl.plot([-20,10],[1.00,1.00],'k--')	# -7 to -3 1979-2016
                #pl.plot([-20,10],[0.72,0.72],'k--')
                #pl.plot([-20,10],[1.28,1.28],'k--')
                pl.plot([-20,10],[0.92,0.92],'k--')     # -7 to -3 1979-2010
                pl.plot([-20,10],[0.66,0.66],'k--')
                pl.plot([-20,10],[1.20,1.20],'k--')
		pl.xlabel('Days after extreme')
		pl.ylabel('Mean number of injection events')
		lg = pl.legend(loc=0,frameon=False,prop={'size':11},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		#pl.ylim(0.4,2.2)
		pl.xlim(-20,10)
		pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/lagged_association_number.pdf',format='pdf')	
                pl.figure(2)
                pl.plot(lags,r1d,'r',linewidth=2,label='Warm extremes')
                pl.plot(lags,r2d,'b',linewidth=2,label='Cold extremes')
                #pl.plot([-20,10],[3.99,3.99],'k--')	# -8 to -2 1979-2016
		#pl.plot([-20,10],[2.86,2.86],'k--')
		#pl.plot([-20,10],[5.28,5.28],'k--')
                #pl.plot([-20,10],[3.69,3.69],'k--')	# -8 to -2 1979-2000
                #pl.plot([-20,10],[2.25,2.25],'k--')
                #pl.plot([-20,10],[5.40,5.40],'k--')
                #pl.plot([-20,10],[3.16,3.16],'k--')	# -7 to -3 1979-2016
                #pl.plot([-20,10],[2.15,2.15],'k--')
                #pl.plot([-20,10],[4.37,4.37],'k--')
                pl.plot([-20,10],[2.93,2.93],'k--')     # -7 to -3 1979-2010
                pl.plot([-20,10],[1.93,1.93],'k--')
                pl.plot([-20,10],[4.06,4.06],'k--')
                pl.xlabel('Days after extreme')
                pl.ylabel('Mean duration of injection events [days]')
                lg = pl.legend(loc=0,frameon=False,prop={'size':11},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
                pl.ylim(0,7)
                pl.xlim(-20,10)
                pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/lagged_association_duration.pdf',format='pdf')
		pl.show()

	def extractIntrusions(self,xs):
		G   = [self.G[i]   for i in xs]
		Q   = [self.Q[i]   for i in xs]
		D   = [self.D[i]   for i in xs]
		LON = [self.LON[i] for i in xs]
		LAT = [self.LAT[i] for i in xs]
		P   = [self.P[i]   for i in xs]
		return G,Q,D,LON,LAT,P

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

	def plotTrackDensity(self,llon,llat,LLON,LLAT,DD,case,plot=True):
		# Compute densities and plot
		nn,ntot,ttot = self.tr.density(llon,llat,nend=None)
		if plot:
			# Compute climatological density
			NN,Ntot,Ttot = self.tr.density(LLON,LLAT,nend=None)
			dtot         = ttot - Ttot
			# Compute trend
			dates          = [d[0] for d in DD]
			years,ns       = self.tr.densityToSeason(NN,dates,YearRange=(1980,2016),Season='NDJFM')
			field,trend,pv = self.tr.detrend2d(years,ns)
			# Plot
			n,x,y  = self.tr.interp2d(trend,self.tr.x,self.tr.y,6,kind='cubic')
	                self.tr.plotMap(x,y,n,20,cmap=pl.cm.RdBu_r,extend='both',cbarlabel='Number density trend {400x400 km$^{2}$}$^{-1}$',\
	                             title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_all_trend.pdf' % (self.Source))
			n,x,y  = self.tr.interp2d(Ntot,self.tr.x,self.tr.y,6,kind='cubic')
			self.tr.plotMap(x,y,n,20,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$}$^{-1}$',\
			             title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_all.pdf' % (self.Source),Tang=Ttot)
			n,x,y  = self.tr.interp2d(ntot,self.tr.x,self.tr.y,6,kind='cubic')
			self.tr.plotMap(x,y,n,20,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$}$^{-1}$',\
			             title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s.pdf' % (self.Source,case),Tang=ttot)
	                self.tr.plotMap(x,y,n,20,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density {400x400 km$^{2}$}$^{-1}$',\
	                             title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s_dtot.pdf' % (self.Source,case),Tang=dtot)
		else:
			return ntot,ttot

	def plotTrackDensity2(self,llon,llat):

		# List of intrusion llon; (N,date,width,tsteps)
		# Need to compute tangents for each trajecotry point
		LON,LAT,XT,YT = [],[],[],[]
		for i in range(len(llon)):
			for t in range(len(llon[i])):
				xs,ys     = self.tr.proj.m(llon[i][t],llat[i][t])
				xs,ys     = xs.mean(axis=0),ys.mean(axis=0)
				yt,xt     = self.tr.getTangents(xs,ys,mag=1)
				lons,lats = self.tr.proj.m(xs,ys,inverse=True)
				LON,LAT   = LON+list(lons),LAT+list(lats)
				XT,YT     = XT+list(xt),YT+list(yt)
		LON,LAT = np.array(LON),np.array(LAT)
		XT,YT   = np.array(XT),np.array(YT)
                # Set up lat/lon grid
                res           = 3.
                nlon          = 360/res
                nlat          = 180/res
                lon           = np.arange(nlon)*res - 180.
                lat           = np.arange(nlat+1)*res - 90.
                lon,lat       = np.meshgrid(lon,lat)	
		Count0  = np.zeros(lon.shape)
		XC,YC  = np.zeros(lon.shape),np.zeros(lon.shape)
		N      = np.zeros(lon.shape)
		XC1,YC1 = [],[]
                # Lagged composite
                for ii in range(len(LON)):
			ix,iy     = np.where(self.haversine(LON[ii],LAT[ii],lon,lat)<564)
			a         = np.zeros(lon.shape)
			a[ix,iy]  = XT[ii]
			b         = np.zeros(lon.shape)
			b[ix,iy]  = YT[ii]
			XC1.append(a)
			YC1.append(b)
			#XC[ix,iy] += XT[ii]
			#YC[ix,iy] += YT[ii]
			#N[ix,iy]  += 1
                	Count0    += np.where(self.haversine(LON[ii],LAT[ii],lon,lat) < 564, 1, 0)
		XC1,YC1 = np.array(XC1).mean(axis=0),np.array(YC1).mean(axis=0)
		Count0  = self.tr.proj(Count0,lon[0,:],lat[:,0])
		XC1     = self.tr.proj(XC1,lon[0,:],lat[:,0])
		YC1     = self.tr.proj(YC1,lon[0,:],lat[:,0])
		#N       = self.tr.proj(N,lon[0,:],lat[:,0])
		XC,YC  = 1.*XC/N,1.*YC/N
		Tang   = np.array([XC1,YC1])
		Tang   = np.rollaxis(Tang,-1,0)
		Tang   = np.rollaxis(Tang,2,1)
		Tang = np.zeros((self.tr.proj.nx,self.tr.proj.ny,2))
		#return Count0,Tang
		self.tr.plotMap(self.tr.x,self.tr.y,Count0,14,cmap=pl.cm.OrRd,extend='max',cbarlabel='Number density',title='test',savename=None,Tang=Tang)

	def plotLagged(self,Dates0=None):
		# Given intrusions, plot lagged figures of moisture etc.
		N,i0      = np.arange(-30,10+1,1),-1
		days      = [self.reds.getDays(*Date0) for Date0 in Dates0]
		datelists = [ [self.reds.getDate(day0+day,nd=24) for day in N] for day0 in days ]
		# DataServers and LambertProjector
		mmt2ds = MMDataServer(Field='ttr', Source='ERAInt', LevType='surface_forecast')
		t2ds   = reDataServer(Field='ttr', Source='ERAInt', LevType='surface_forecast')
		spds   = reDataServer(Field='slp', Source='ERAInt', LevType='surface_analysis')
		proj   = LambertProjector(boundinglat=0.01,resolution=100.)
		# Clim
		t2clim  = np.array([mmt2ds.getSeason(Year=year,Season='NDJFM') for year in range(1979,2012+1,1)]).mean(axis=0)
		#t2clim = proj(np.array([mmt2ds.getMonth(Year=year,Month=datelist[0][1]) for year in range(1979,2012+1,1)]).mean(axis=0), mmt2ds.lon, mmt2ds.lat)
		T,SLP = [],[]
		for datelist in datelists:
			print datelist[i0]
			T.append([t2ds.snapshot(*date) for date in datelist])
			SLP.append([spds.snapshot(*date) for date in datelist])
#			T.append(proj(np.array([t2ds.snapshot(*date) for date in datelist]),t2ds.lon,t2ds.lat))
#			SLP.append(proj(np.array([spds.snapshot(*date) for date in datelist]),spds.lon,spds.lat)/100.)
		T,SLP = np.array(T).mean(axis=0),np.array(SLP).mean(axis=0)
		T     = T - t2clim[np.newaxis,:,:]
		# Plot
		#cseq1 = np.arange(-16,16+2,2)
		cseq1 = np.arange(-25,25+2.5,2.5)
		cseq2 = np.arange(960,1050+10,10)
		m     = Basemap(projection='ortho',lon_0=110,lat_0=25)
		tj,tk = m(*np.meshgrid(t2ds.lon,t2ds.lat))
		sj,sk = m(*np.meshgrid(spds.lon,spds.lat))
		for ii in range(len(T)):
			cf   = m.contourf(tj,tk,T[ii],cseq1,cmap=pl.cm.RdBu,extend='both')
			cl   = m.contour(sj,sk,SLP[ii],cseq2,colors='k',linewidths=1,alpha=0.3)
			cbar = pl.colorbar(cf)
			cbar.set_label('OLR anomaly [W m$^{-2}$]')
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=10)
			m.drawcoastlines()
			m.drawparallels([70,80],latmax=90)
			pl.title(N[ii])
			pl.savefig('figs/lagged/comps/%s.pdf' % (ii+1),format='pdf')
			pl.close()

	def makeRandDateList(self,N,daylist,space):
		# Add first random day
		daysrand = np.array([random.choice(daylist)])
		while len(daysrand) < N:
		        day0r = random.choice(daylist)
		        if np.abs(daysrand-day0r).min() > space:
		        	daysrand = np.append(daysrand,np.array([day0r]),axis=0)
		daysrand = np.array(daysrand)
		daysrand.sort()
		dt        = np.diff(daysrand)
		avg,std   = dt.mean(),dt.std()
		datesrand = [self.reds.getDate(day,nd=24) for day in daysrand]
		return datesrand,avg,std,stats.skew(dt)

	def imposeStats(self,N,daylist,space,avg0,std0,skw0):
		avg,std,skw = -999,-999,-999
		#from statsmodels.sandbox.distributions.extras import pdf_mvsk
		#pf = pdf.mvsk([avg0,std0,skw0,kur0])
#		while (np.abs((avg-avg0)/avg0) > 0.01) or (np.abs((std-std0)/std0) > 0.01) or (np.abs((skw-skw0)/skw0) > 0.01):
		datesrand,avg,std,skw = self.makeRandDateList(N,daylist,space)
		return datesrand,avg,std,skw

	def spacing(self,Sector=(0,360),case='cold'):
		# Average spacing between intrusions and warm events
		# Using random dates, same spacing or larger? Is the presence of an injection 
		# preceeding warm even significant? Number injections existing within in certain time?
		dates,datelist = self.topEvents(50,case=case)
		# Extract YearRange
		mind,maxd      = -7,-3
		years          = range(1979,2016+1,1)
		dates,datelist = [date for date in dates if date[0] in years],[date for date in datelist if date[0] in years]
		print len(dates)
		daylist        = [self.reds.getDays(*date) for date in datelist]
		# Associate events
		xs0,ls0        = self.associateEvents(dates,mind=mind,maxd=maxd)
		#mxl            = [self.maxlons[xs0i] for xs0i in xs0]
		#ll0            = self.filterLons(mxl,Sector=Sector)	
		#ll0            = np.array([len(ls0[i]) for i in range(len(ls0))])
		#xs0            = ll0.mean()
		#xv0            = ll0.std()
		dr0           = [ np.sum([len(self.D[xxsi]) for xxsi in xxs]) for xxs in ls0 ]
		xs0           = np.mean(dr0)/4.
		print '\nxs0   = %s' % (xs0)
		# Random series of dates
		r,Avg,Std,Skw = [],[],[],[]
		while len(r) < 10000:	
			#datesrand,avg,std,skw = self.makeRandDateList(50,daylist,7)
			datesrand,avg,std,skw = self.imposeStats(50,daylist,7,avg0=None,std0=None,skw0=None)
			xs,ls		  = self.associateEvents(datesrand,mind=mind,maxd=maxd)
			#mxl		  = [self.maxlons[xsi] for xsi in xs]
			#ll		  = self.filterLons(mxl,Sector=Sector)
                	#ll                = np.array([len(ls[i]) for i in range(len(ls))])
                	#xs                = ll.mean()
                	#xv                = ll.std()
			dr		  = [ np.sum([len(self.D[xxsi]) for xxsi in xxs]) for xxs in ls]
			xs		  = np.mean(dr)/4.
			r.append(xs)
			Avg.append(avg)
			Std.append(std)
			Skw.append(skw)
		print 'mean   = %s' % (np.mean(r))
		print 'std    = %s' % (np.std(r))
		print 'min    = %s' % (np.min(r))
		print 'argmin = %s' % (np.argmin(r))
		print 'max    = %s' % (np.max(r))
		print '1st    = %s' % stats.scoreatpercentile(r,1)
		print '50th   = %s' % stats.scoreatpercentile(r,50)
		print '99th   = %s' % stats.scoreatpercentile(r,99)
		print 'P(xs0) = %s' % stats.percentileofscore(r,xs0)
		print 'Average spacing = %s' % (np.mean(Avg))
		print 'Average std of spacing = %s' % (np.mean(Std))
		print 'Average skew of spacing = %s' % (np.mean(Skw))

	def filterLons(self,lons,Sector=(330,105)):
		ll = []
                for i in range(len(lons)):
                        if Sector[0] > Sector[1]:
                                if (Sector[0] <= lons[i] <= 360) or (0 <= lons[i] < Sector[1]):
					ll.append(lons[i])
                        elif Sector[0] < Sector[1]:
                                if (Sector[0] <= lons[i] < Sector[1]):
					ll.append(lons[i])
                return ll

	def pdf(self):
		i0 = 859+1
		# Events
                Iw,vr,vtn,dates,datelist = self.topEvents(50,case='warm')
                xs,ls = self.associateEvents(dates)
		# Durations
		t1 = np.array([len(d) for d in self.D[0:i0]])/4.
		t2 = np.array([len(d) for d in self.D[i0:]])/4.
		t3 = np.array([len(self.D[xi]) for xi in xs])/4.
		# Widths
		w1 = [ np.mean([len(self.G[i][j]) for j in range(len(self.G[i]))]) for i in range(0,i0)]
		w2 = [ np.mean([len(self.G[i][j]) for j in range(len(self.G[i]))]) for i in range(i0,len(self.G)-1)]
		w3 = [ np.mean([len(self.G[i][j]) for j in range(len(self.G[i]))]) for i in xs]
		# Flux
                q1 = np.array([ np.sum([np.sum(self.Q[i][j]) for j in range(len(self.Q[i]))]) for i in range(0,i0)])/4.
		q2 = np.array([ np.sum([np.sum(self.Q[i][j]) for j in range(len(self.Q[i]))]) for i in range(i0,len(self.Q)-1)])/4.
                q3 = np.array([ np.sum([np.sum(self.Q[i][j]) for j in range(len(self.Q[i]))]) for i in xs])/4.
		# PDF
		binst   = np.arange(1.5,12+0.5,0.5)
		binsw   = np.arange(9,60+4,4)
		binsq   = 20
		ht1,et1 = np.histogram(t1,bins=binst,normed=True)
		ht2,et2 = np.histogram(t2,bins=binst,normed=True)
		ht3,et3 = np.histogram(t3,bins=binst,normed=True)
                hw1,ew1 = np.histogram(w1,bins=binsw,normed=True)
                hw2,ew2 = np.histogram(w2,bins=binsw,normed=True)
		hw3,ew3 = np.histogram(w3,bins=binsw,normed=True)
                hq1,eq1 = np.histogram(q1,bins=binsq,normed=True)
                hq2,eq2 = np.histogram(q2,bins=binsq,normed=True)
		hq3,eq3 = np.histogram(q3,bins=binsq,normed=True)
		# Plot
		pl.figure(1)
		pl.plot(et1[0:-1],ht1,'k',linewidth=1.5,label='1950-1983')
		pl.plot(et2[0:-1],ht2,'r',linewidth=1.5,label='1984-2016')
		pl.plot(et3[0:-1],ht3,'b--',linewidth=1.5,label='Top 50')
		pl.xlabel('Duration [days]')
		pl.ylabel('Normalised frequency')
		pl.grid()
		lg = pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('figs/dur.norm.pdf',format='pdf')
                pl.figure(2)
                pl.plot(ew1[0:-1],hw1,'k',linewidth=1.5,label='1950-1983')
                pl.plot(ew2[0:-1],hw2,'r',linewidth=1.5,label='1984-2016')
		pl.plot(ew3[0:-1],hw3,'b--',linewidth=1.5,label='Top 50')
                pl.xlabel('Width [degrees longitude]')
                pl.ylabel('Normalised frequency')
                pl.grid()
		lg = pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('figs/width.norm.pdf',format='pdf')
                pl.figure(3)
                pl.plot(eq1[0:-1],hq1,'k',linewidth=1.5,label='1950-1983')
                pl.plot(eq2[0:-1],hq2,'r',linewidth=1.5,label='1984-2016')
		pl.plot(eq3[0:-1],hq3,'b--',linewidth=1.5,label='Top 50')
                pl.xlabel('Flux [Tg]')
                pl.ylabel('Normalised frequency')
                pl.grid()
		lg = pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('figs/F.norm.pdf',format='pdf')

        def interpolate(self,vq):
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
		return vq

	def randDensity(self,N=100,Sector=(0,360)):
		mind,maxd = -7,-3
		YearRange = [1979,2000]
		Nrand     = 27
		years     = range(YearRange[0],YearRange[1]+1,1)
		if YearRange == [1979,2016]:
			if [mind,maxd] == [-8,-2]: clim_factor = 1.28	#  Full period: 1979 - 2016
			if [mind,maxd] == [-7,-3]: clim_factor = 1.00   #  Full period: 1979 - 2016
			if [mind,maxd] == [-6,-2]: clim_factor = 1.00   #  Full period: 1979 - 2016
                if YearRange == [1979,2000]:
                        if [mind,maxd] == [-8,-2]: clim_factor = 1.20   #  Early period: 1979 - 2000
                        if [mind,maxd] == [-7,-3]: clim_factor = 0.92   #  Early period: 1979 - 2000
			if [mind,maxd] == [-6,-2]: clim_factor = 0.92   #  Early period: 1979 - 2000
		# Cold and warm events density
		datesc,datelist          = self.topEvents(50,case='cold')
		datesw,datelist          = self.topEvents(50,case='warm')
		datesc                   = [date for date in datesc if date[0] in years]
		datesw                   = [date for date in datesw if date[0] in years]
		datesNc                  = len(datesc)
		datesNw                  = len(datesw)
		#Nrand                   = np.array([datesNc,datesNw]).min()
		#datesc,datesw           = datesc[0:Nrand],datesw[0:Nrand]
		#datesNc                 = len(datesc)
                #datesNw                 = len(datesw)
		xs,ls                    = self.associateEvents(datesc,mind=mind,maxd=maxd)
		days        		 = np.array([self.reds.getDays(*date) for date in datesc])
		days.sort()
		dtcold          	 = np.diff(days)
		avgcold,stdcold   	 = dtcold.mean(),dtcold.std()
                G,Q,D,LON,LAT,P          = self.extractIntrusions(xs)
                gg,qq,dd,llon,llat,pp    = self.tr.filterInjections(G,Q,D,LON,LAT,P,Sector=Sector) 
                ntotcold,ttotcold        = self.plotTrackDensity(llon,llat,LLON=None,LLAT=None,DD=None,case='cold',plot=False)
		# Warm events density
		xs,ls                    = self.associateEvents(datesw,mind=mind,maxd=maxd)
                days                     = np.array([self.reds.getDays(*date) for date in datesw])
                days.sort()
                dtwarm                   = np.diff(days)
                avgwarm,stdwarm,skwwarm  = dtwarm.mean(),dtwarm.std(),stats.skew(dtwarm)
                G,Q,D,LON,LAT,P          = self.extractIntrusions(xs)
                gg,qq,dd,llon,llat,pp    = self.tr.filterInjections(G,Q,D,LON,LAT,P,Sector=Sector)
                ntotwarm,ttotwarm        = self.plotTrackDensity(llon,llat,LLON=None,LLAT=None,DD=None,case='warm',plot=False)
		# Climatological density
		Nclim,Tclim = self.plotTrackDensity(self.LON,self.LAT,LLON=None,LLAT=None,DD=None,case='rand',plot=False)
		print 'warm: %s' % (datesNw)
		print 'cold: %s' % (datesNc)
		# Random sampling
                r,daylist,Narr,Tarr      = [],[self.reds.getDays(*date) for date in datelist],[],[]
                for t in range(N):
			print '%s of %s' % (t+1,N)
                        # Add first random day
			datesrand,avg,std,skw = self.imposeStats(Nrand,daylist,7,avgwarm,stdwarm,skwwarm)
                        xs,ls 		      = self.associateEvents(datesrand,mind=mind,maxd=maxd)
			G,Q,D,LON,LAT,P       = self.extractIntrusions(xs)
			gg,qq,dd,llon,llat,pp = self.tr.filterInjections(G,Q,D,LON,LAT,P,Sector=Sector)
			ntot,ttot             = self.plotTrackDensity(llon,llat,LLON=None,LLAT=None,DD=None,case='rand',plot=False)
			Narr.append(ntot)
			Tarr.append(ttot)
                print 'warm: %s' % (datesNw)
                print 'cold: %s' % (datesNc)
		Narr,Tarr   = np.array(Narr),np.array(Tarr)
		Narrm,Tarrm = Narr.mean(axis=0),Tarr.mean(axis=0)
		nwarm,ncold = np.zeros(ntotwarm.shape),np.zeros(ntotcold.shape)
		for i in range(len(nwarm)):
			for j in range(len(nwarm[i])):
				nwarm[i,j] = stats.percentileofscore(Narr[:,i,j],ntotwarm[i,j])
				ncold[i,j] = stats.percentileofscore(Narr[:,i,j],ntotcold[i,j])				

		Nclim0    = clim_factor*Nclim/len(self.G)	# clim_factor x clim per injection event
                dntotwarm = ntotwarm - datesNw*Nclim0       	# Anomaly wrt clim
                dntotcold = ntotcold - datesNc*Nclim0       	# Anomaly wrt clim
                dTwarm    = ttotwarm - Tclim
                dTcold    = ttotcold - Tclim

                nwarm_i,x,y   = self.tr.interp2d(nwarm,self.tr.x,self.tr.y,6,kind='cubic')
                ncold_i,x,y   = self.tr.interp2d(ncold,self.tr.x,self.tr.y,6,kind='cubic')
                nwarm_t       = np.tile(nwarm[:,:,np.newaxis],(1,1,2))
                ncold_t       = np.tile(ncold[:,:,np.newaxis],(1,1,2))
		# Plot
		cseqr0 = np.arange(0,1+0.1,0.1)
		cseqr  = np.arange(0,6+1,1)
		cseqa  = np.arange(-6,6+1,1)

                n,x,y  = self.tr.interp2d(Nclim0,self.tr.x,self.tr.y,6,kind='cubic')
                self.tr.plotMap(x,y,n,cseqr,cmap=pl.cm.OrRd,extend='max',cbarlabel=r'Number density $%s\times${$350 \times 350$ km$^{2}$ injection}$^{-1}$' % (clim_factor),\
                                title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s_clim.pdf' % (self.Source,'all'),Tang=Tclim)

		n,x,y  = self.tr.interp2d(Narrm,self.tr.x,self.tr.y,6,kind='cubic')
               	self.tr.plotMap(x,y,n/Nrand,cseqr,cmap=pl.cm.OrRd,extend='max',cbarlabel=r'Number density {$350 \times 350$ km$^{2}$}$^{-1}$',\
				title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s_rand.pdf' % (self.Source,'all'),Tang=Tarrm)

                n,x,y  = self.tr.interp2d(dntotcold,self.tr.x,self.tr.y,6,kind='cubic')
		n      = np.ma.masked_where(ncold_i>3,n)
		dTcold = np.ma.masked_where(ncold_t>3,dTcold)
                self.tr.plotMap(x,y,n/datesNc,cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=r'Number density {$350 \times 350$ km$^{2}$}$^{-1}$',\
                                title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s_rand.pdf' % (self.Source,'cold'),\
				Tang=dTcold,stip=None)

                n,x,y  = self.tr.interp2d(dntotwarm,self.tr.x,self.tr.y,6,kind='cubic')
		n      = np.ma.masked_where(nwarm_i<97,n)
		dTwarm = np.ma.masked_where(nwarm_t<97,dTwarm)
                self.tr.plotMap(x,y,n/datesNw,cseqa,cmap=pl.cm.RdBu_r,extend='both',cbarlabel=r'Number density {$350 \times 350$ km$^{2}$}$^{-1}$',\
                                title=self.Source,savename='/mnt/climstorage/cian/WarmArctic/figs/%s_%s_rand.pdf' % (self.Source,'warm'),\
				Tang=dTwarm,stip=None)

        def intrusionsTimesteps(self):
                # Returns (time)x(coordinate) array of all trajectory points
                # Find trajectory points
                time,points = [],[]
		nsteps      = np.linspace(-20,20,41)
                for i in range(len(self.LON)): 
                        for j in range(len(self.LON[i])):
					lons,lats = self.LON[i][j],self.LAT[i][j]
					xs,ys     = self.tr.proj.m(lons,lats)
					xs,ys     = xs.mean(axis=0),ys.mean(axis=0)	
					#lon,lat   = self.tr.proj.m(xs,ys,inverse=True)  
                                        for k in range(len(nsteps)): 
                                                hour  = self.reds.getHours(*self.D[i][j]) + nsteps[k]*6
						date  = self.reds.getDate(hour)
                                                time.append(date)
                                                points.append([xs[k],ys[k]])
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

        def interpTrack(self,xs,n):
                xold = np.linspace(0,1,len(xs))
                xnew = np.linspace(0,1,len(xs)*n)
                f  = interpolate.interp1d(xold,xs)
                xs = f(xnew)
                return xs

	def density0(self,xs,ys):
                # Data holders
                N   = np.zeros((self.tr.nx,self.tr.ny))
                for i in range(len(xs)):
                        xi,yi    = np.argmin(abs(self.tr.x-xs[i])),np.argmin(abs(self.tr.y-ys[i])) 
                        N[yi,xi] = N[yi,xi] + 1
		return N

	def intrusionLaggedDates(self,lags=[-10,0],case='warm',rankN=50,plot=True,YearRange=(1979,2016)):
		# computes density of trajectory particles are each lagged timestep
		# Time and points of all trajectory particles
	        t,p = self.intrusionsTimesteps()
	        t,p = self.sortTimePoint(t,p)
		# All points list
		pall = np.array([p[i][j] for i in range(len(p)) for j in range(len(p[i])) if t[i][0] in range(YearRange[0],YearRange[1]+1,1)])
		tall = [t[i] for i in range(len(t)) if t[i][0] in range(YearRange[0],YearRange[1]+1,1)]
		X,Y  = pall[:,0],pall[:,1]
		# Dates to composite over
		days  = np.arange(lags[0],lags[1]+0.25,0.25)
		Xs,Ys = [[] for k in days],[[] for k in days]
		dates,datelist = self.topEvents(rankN,case)
		dates  = [date for date in dates if date[0] in range(YearRange[0],YearRange[1]+1,1)]
		datesN = len(dates)
		for date0 in dates:
			# Lagged datelist
			dlist = [self.reds.getDate(self.reds.getHours(*date0) + k*24) for k in days]
			for ii in range(len(dlist)):
				try:
					points = p[t.index(dlist[ii])]
					for jj in range(len(points)):
						Xs[ii].append(points[jj][0])
						Ys[ii].append(points[jj][1])	
				except:
					pass
		# Lagged composite
		Nl = []
		for ii in range(len(Xs)):
			Count0 = self.density0(Xs[ii],Ys[ii])
			Nl.append(Count0)
		Nl  = 1.*np.array(Nl)/datesN
                # Climatology
		climname = '/mnt/climstorage/cian/WarmArctic/intrusion_clim_%s-%s.p' % (YearRange[0],YearRange[1])
		if not os.path.isfile(climname):
			print '\nMaking climatological intrusion density file %s ...\n' % (climname)
			Nf = self.density0(X,Y)
			Nf = 1.*Nf/len(tall)
			toPick(Nf,climname)
		else:
			Nf = unpick(climname)	
		Nl  = Nl - Nf
		# Plot
		if plot:
			cseqa = np.arange(-0.25,0.25+0.025,0.025)
			cseqf = np.arange(0,1+0.1,0.1)
			for ii in range(len(Nl)):
				cf   = pl.contourf(self.tr.x,self.tr.y,Nl[ii],cseqa,cmap=pl.cm.RdBu_r,extend='both')
				cbar = pl.colorbar(cf)
				cbar.set_label('Number density {350x350 km$^{2}$}$^{-1}$')
				#self.tr.proj.m.drawcoastlines(color='0.5',linewidth=0.7)
				drawCoastlinesNoRivers(self.tr.proj.m,color='0.5',linewidth=0.7)
				self.tr.proj.m.drawparallels([70,80,85],latmax=90)
				pl.title('day %s' % (days[ii]))
				pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/lagged/intrusions/%s/%s.pdf' % (case,ii+1))
				pl.close()
			Nf,x,y = self.tr.interp2d(Nf,self.tr.x,self.tr.y,6,kind='cubic')
			cf     = pl.contourf(x,y,4*Nf,cseqf,cmap=pl.cm.OrRd,extend='max')
			cbar   = pl.colorbar(cf)
			cbar.set_label(r'Number density {$350 \times 350$ km$^{2}$ day}$^{-1}$')
			#self.tr.proj.m.drawcoastlines(color='0.5',linewidth=0.7)
			drawCoastlinesNoRivers(self.tr.proj.m,color='0.5',linewidth=0.7)
			self.tr.proj.m.drawparallels([70,80],latmax=90)
			pl.title('Climatology')
			pl.savefig('/mnt/climstorage/cian/WarmArctic/figs/lagged/intrusions/%s/clim.pdf' % (case))
			pl.close()
			return Nl
		else:
			return Nl

	def particleDensity(self,lons,lats):
		xs,ys = self.tr.proj.m(lons,lats)
                # Data holders
                N   = np.zeros((self.tr.nx,self.tr.ny)) 
                for i in range(len(xs)):
                        xi,yi    = np.argmin(abs(self.tr.x-xs[i])),np.argmin(abs(self.tr.y-ys[i])) 
                        N[yi,xi] = N[yi,xi] + 1 
                return N

	def lagged1d(self,case='warm'):
		# Dates
		Iw,vr,vtn,dates,datelist = self.topEvents(50,case=case)
		days = np.arange(-8,8+0.25,0.25)
		# DataServers and LambertProjector
		proj = LambertProjector(boundinglat=80,resolution=80.)
		ds1  = reDataServer(Field='nlwrs',Source=self.Source,LevType='surface')
		ds2  = reDataServer(Field='t2m',Source=self.Source,LevType='surface')
		# Composite
		s1,s2 = [],[]
		for date in dates:
			print date
			hour = self.reds.getHours(*date)
			s1.append([ds1.snapshot(*self.reds.getDate(hour+day*24)) for day in days])
			s2.append([ds2.snapshot(*self.reds.getDate(hour+day*24)) for day in days])
		s1,s2 = np.array(s1).mean(axis=0),np.array(s2).mean(axis=0)
		s1,s2 = proj(s1,ds1.lon,ds1.lat).mean(axis=-1).mean(axis=-1),proj(s2,ds2.lon,ds2.lat).mean(axis=-1).mean(axis=-1)
		fig,ax1 = pl.subplots()
		ax2     = ax1.twinx()

#		VQH,VQR = np.array(VQH).mean(axis=0),np.array(VQR).mean(axis=0)
#		ax1.plot(edgesvq[0:-1],VQH,'b',linewidth=2)
#		ax1.plot(edgesvq[0:-1],VQR,'r',linewidth=2)
#		ax2.plot(edgesvq[0:-1],VQR/VQH,'k--',linewidth=1)
#		pl.xlim(edgesvq[0],edgesvq[-2])
#		ax2.set_ylim(0,200)
#		ax1.set_xlabel('Standard deviation')
#		pl.title('Moisture flux: %s-%sE' % (lon1,lon2))
#		pl.savefig('figs/pdf/vq.pdf.%s-%sE.pdf' % (lon1,lon2))
#		pl.show()

		ax1.plot(days,s1,'b')
		ax2.plot(days,s2,'r')
		pl.show()

	def rankIntrusionDays(self,metric='number'):	
		# Rank association for all days
		dates,datelist = self.topEvents(50,case='warm')
		R,index = [],range(1,len(datelist)+1,1)
		if metric == 'number':
			for date in datelist:
				xs,ls = self.associateEvents([date],-8,-2)
				R.append(len(xs))
		elif metric == 'duration':
			for date in datelist:
				xs,ls = self.associateEvents([date],-8,-2)
				xs    = np.mean([ np.sum([len(self.D[xxsi]) for xxsi in xxs]) for xxs in ls])
				R.append(xs)
		# Rank events
		rank_index = [i for (r,i) in sorted(zip(R,index))][::-1]
		rank_dur   = [r for (r,i) in sorted(zip(R,index))][::-1]
		# Make File
		f=open('ranked_%s.txt' % (metric),'w')
		f.write('index|date|%s\n' % (metric))
		for i in range(len(rank_index)):
			date = str(datelist[rank_index[i]-1])
			f.write('%s|%s|%s\n' % (rank_index[i],date,rank_dur[i]))
		f.close()

	def makeSatelliteFile(self):
		FileName = '../IntrusionTrajectoryPoints_%sx%skm_%sN.nc' % (self.res,self.res,self.blat)
		if not os.path.isfile(FileName):
			# Time ane points of centroid trajectories
			t,p     = self.intrusionsTimesteps()
                	t,p     = self.sortTimePoint(t,p)
			# Datelist for NDJFM
			dates0 = self.reds.getDateList(Year=1979,Season='JFM')
			dates1 = self.reds.getDateList(Year=2016,Season='ND')
			dates  = [self.reds.getDateList(Year=year,Season='NDJFM') for year in range(1980,2016+1,1)]
			dates  = dates0 + [dates[i][j] for i in range(len(dates)) for j in range(len(dates[i]))] + dates1
			times  = [self.reds.getHours(*date) for date in dates]
			# Density for each timestep
			N  = []
			n0 = np.zeros((self.tr.proj.nx,self.tr.proj.ny))
			for i in range(len(dates)):
				if dates[i] in t:
					j     = t.index(dates[i])
					p0    = np.array(p[j])
					xs,ys = p0[:,0],p0[:,1]
					n     = self.density0(xs,ys)
					n[np.where(n>0)] = 1
					N.append(n)
				else:
					N.append(n0)
			N = np.array(N)
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
                        File.createDimension('X',self.tr.proj.nx)
                        File.createDimension('Y',self.tr.proj.ny)
                        var           = File.createVariable('lat','f',('X','Y',))
                        var.long_name = 'Latitude'
                        var.units     = 'degrees_north'
                        var[:]        = self.tr.proj.lat.astype('f')
                        var           = File.createVariable('lon','f',('X','Y',))
                        var.long_name = 'Longitude'
                        var.units     = 'degrees_east'
                        var[:]        = self.tr.proj.lon.astype('f')
                        # Create Variables
                        var           = File.createVariable('index','f',('time','X','Y',))
                        var.long_name = 'Intrusion centroid trajectory at girdpoint'
                        var.units     = '0/1'
                        var[:]        = N
                        # Close file
                        File.close()
		else:
			N = Dataset(FileName,'r').variables['index'][:]
                # Plot
		self.tr.plotMap(self.tr.proj.x,self.tr.proj.y,N.mean(axis=0),14,cmap=pl.cm.OrRd,extend='max',\
	        		cbarlabel=r'Number density {80 $\times$ 80 km$^{2}$}$^{-1}$',title=self.Source)

	def sectorDensity(self):

		LON = [self.LON[i] for i in range(len(self.maxlons)) if (self.maxlons[i]<60) or (self.maxlons[i]>340)]
		LAT = [self.LAT[i] for i in range(len(self.maxlons)) if (self.maxlons[i]<60) or (self.maxlons[i]>340)]

		inds = []
		XS,YS = np.zeros((0,28)),np.zeros((0,28))
		for i in range(len(LON)):
			for t in range(len(LON[i])):
				if LAT[i][t].min()>0:
					xs,ys   = self.tr.proj.m(LON[i][t],LAT[i][t])
					xs,ys   = xs.mean(axis=0)[:28],ys.mean(axis=0)[:28]
					lon,lat = self.tr.proj.m(xs,ys,inverse=True)
					try:
						i0      = np.where(lat>80)[0][0]
						i1      = np.where(lat[i0:]<80)[0][0]
						dd      = np.sqrt((xs[i0]-xs[i1])**2 + (ys[i0]-ys[i1])**2)
						if dd > 1200e03:
							XS,YS = np.append(XS,xs[np.newaxis,:],axis=0),np.append(YS,ys[np.newaxis,:],axis=0)
					except:
						pass
					#XS,YS = np.append(XS,xs[np.newaxis,:],axis=0),np.append(YS,ys[np.newaxis,:],axis=0)
		N               = len(XS)
		xmean,ymean     = np.mean(XS,axis=0),np.mean(YS,axis=0)
		lonmean,latmean = self.tr.proj.m(xmean,ymean,inverse=True)
		xmed,ymed       = np.median(XS,axis=0),np.median(YS,axis=0)
		lonmed,latmed   = self.tr.proj.m(xmed,ymed,inverse=True)
		lon,lat         = lonmed,latmed

		basename = 'median_trajectory_cross80N_20W-60E'

		"""
		# Intial points
		x0,y0 = self.tr.proj.m(0,70)
		x1,y1 = self.tr.proj.m(-132,70)
		XS = np.array([x0/1e06,4.47,4.65,4.57,4.33,4.00,3.43,2.83,x1/1e06])*1e06
		YS = np.array([y0/1e06,2.20,3.03,3.66,3.97,4.27,4.65,5.00,y1/1e06])*1e06
		"""

		XS,YS = xmed,ymed
		f     = interpolate.interp1d(YS,XS,kind='linear')
		XNEW  = np.linspace(YS.min(),YS.max(),10000)
		YNEW  = f(XNEW)

		dist = 400e03
		c    = 0
		ind  = [0]
		for i in range(len(XNEW)):
			xc,yc = XNEW[c],YNEW[c]
			xi,yi = XNEW[i],YNEW[i]	
			if np.sqrt((xi - xc)**2 + (yi-yc)**2) >= dist:
				print np.sqrt((xi - xc)**2 + (yi-yc)**2)
				ind.append(i)
				c = i
		x,y     = [YNEW[i] for i in ind],[XNEW[i] for i in ind]
		lon,lat = self.tr.proj.m(x,y,inverse=True)

		f = open('%s.txt' % (basename),'w')
		f.write('Longitude | Latitude\n')
		for i in range(len(lon)):
			line = '%s | %s\n' % (round(lon[i],2),round(lat[i],2))
			f.write(line)
		f.close()

		pl.plot(x,y,'k',linewidth=1.75,alpha=0.85)
		pl.plot(x,y,'k.',linewidth=1.75,alpha=0.85)
                #pl.plot(XS,YS,'b',linewidth=1.5,alpha=0.65)
                #pl.plot(XS,YS,'b.',linewidth=1.5,alpha=0.65)
		drawCoastlinesNoRivers(self.tr.proj.m,color='0.5',linewidth=0.7)
		self.tr.proj.m.drawparallels([70,80],latmax=90)
		pl.savefig('figs/%s.pdf' % (basename),format='pdf')
		pl.show()


if __name__ == "__main__":

	OpenTraj = str(sys.argv[1])=='True'
	IN       = Intrusions(Source='ERAInt',OpenTraj=OpenTraj)
 
	#IN.sectorDensity()

#	IN.makeSatelliteFile()
#	IN.plotTrackDensity(IN.LON[436:437+1],IN.LAT[436:437+1],IN.LON,IN.LAT,IN.D,case='warm',plot=True)
#	IN.plotTrackDensity2(IN.LON[:],IN.LAT[:])

#	IN.intrusionLaggedDates(lags=[-10,0],case='warm',rankN=50,plot=True,YearRange=(1989,2000))

#	dates,datelist = IN.topEvents(50,'warm')
#	dates,datelist = IN.topEvents(1000,'cold')
#	print dates

#	IN.rankIntrusionDays(metric='duration')

#	IN.makeDataset1(case='warm')
#	IN.makeDataset2(case='warm')
	IN.lagAssociation()
#	IN.lagged1d(case='cold')

#	IN.pdf()
#	IN.randDensity(N=100)

#	case           = 'cold'
#	dates,datelist = IN.topEvents(50,case=case)
#	xs,ls          = IN.associateEvents(datelist)
#	ls0            = np.array([len(i) for i in ls])#.mean()
#	xs0            = np.array([ np.sum([len(IN.D[xxsi]) for xxsi in xxs]) for xxs in ls])#.mean()
#	print ls0.mean(),xs0.mean()/4.
#	print ls0.std(),xs0.std()/4.

#	G,Q,D,LON,LAT,P          = IN.extractIntrusions(xs)
#	gg,qq,dd,llon,llat,pp    = IN.tr.filterInjections(G,Q,D,LON,LAT,P,Sector=(0,360),YearRange=(1950,2017))
#	GG,QQ,DD,LLON,LLAT,PP    = IN.tr.filterInjections(IN.G,IN.Q,IN.D,IN.LON,IN.LAT,IN.P,Sector=(0,360),YearRange=(1950,2017))
#	dates = [dd[i][0] for i in range(len(dd))]
#	dates = [date for date in dates if 1980 <= date[0] <= 2012]
#	IN.plotLagged(Dates0=dates)
#	IN.plotTrackDensity(llon,llat,LLON,LLAT,DD,case=case,plot=True)

#	IN.spacing(Sector=(0,360))
#	IN.plotInjectionSeries(years=range(1980,2015+1,1),detrend=False,Slice=['NDJFM',(11,3)],blat=80)

