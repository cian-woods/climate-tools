from TrackDataServer import *
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from ReanalysisDataServer import DataServer as reDataServer
from LambertProjector import *
from drawCoastlinesNoRivers import *

import numpy as np
import matplotlib.pyplot as pl

class SelectTracks:

	def __init__(	self	):
		# Attributes
		self.proj           = LambertProjector(boundinglat=80,resolution=200.)	# Projection of initial points of trajectories
		self.x,self.y       = self.proj.x[0,:],self.proj.y[:,0]
		self.ix_            = np.where(self.proj.lat.reshape(-1)>=85)[0]	# Only points north of 80N
		#self.mask          = (self.proj.lat.reshape(-1)>83)&(self.proj.lon.reshape(-1)>-140)&(self.proj.lon.reshape(-1)<-60)
		#self.ix_           = np.where(self.mask==True)[0]
		self.proj_          = LambertProjector(boundinglat=60,resolution=400.)	# Projection for density figures
		self.nsq            = self.proj_.nx**2
		self.x_,self.y_     = self.proj_.x[0,:],self.proj_.y[:,0]
		self.d              = DataServer(Type='back2')
		self.hours          = np.arange(0,240+6,6)
		self.weights        = self.hours[::-1]**np.e
		self.weights        = np.ones(len(self.hours))
		self.mmds           = MMDataServer(Field='pw')
		self.reds           = reDataServer(Field='pw',LevType='surface_analysis',LatRange=(50,90))
		self.s1             = self.pw()
		self.interpfunc     = interpolate.interp2d(self.x_,self.y_,self.s1, kind='linear')

	def pdf(self,Years1,Years2):
		pw1    = [self.proj(self.reds.getDataSnaps(Year=year,Season='DJF'),self.reds.lon,self.reds.lat)[:,6,3].reshape(-1) for year in Years1]
		pw1    = [pw1[i][j] for i in range(len(pw1)) for j in range(len(pw1[i]))]
		"""
		pw2    = [self.proj(self.reds.getDataSnaps(Year=year,Season='DJF'),self.reds.lon,self.reds.lat)[:,6,3].reshape(-1) for year in Years2]
		pw2    = [pw2[i][j] for i in range(len(pw2)) for j in range(len(pw2[i]))]
		h1,e1  = np.histogram(pw1,50,normed=True)
		h2,e2  = np.histogram(pw2,50,normed=True)

		pl.plot(e1[0:-1],h1,'b',linewidth=1.25,alpha=0.65)
		pl.plot(e2[0:-1],h2,'r',linewidth=1.25,alpha=0.65)
		pl.grid()
		pl.ylabel('Frequency')
		pl.xlabel('%s [%s]' % (self.reds.long_name,self.reds.units))
		pl.show()
		"""
		return pw1

	def pw(self):
		Years1,Years2 = range(1981,1998+1,1),range(1999,2012+1,1)
		s1 = self.proj_(np.array([self.mmds.getSeason(Year=year,Season='DJF') for year in Years1]).mean(axis=0),self.mmds.lon,self.mmds.lat)
		#s2 = self.proj_(np.array([self.mmds.getSeason(Year=year,Season='DJF') for year in Years2]).mean(axis=0),self.mmds.lon,self.mmds.lat)
		#ds = s2-s1

		"""
		cseqa,cseqf = np.arange(-1,1+0.2,0.2),np.arange(0,10+1,1)
		pl.figure(1)
		pl.contourf(self.proj_.x,self.proj_.y,ds,cseqa,cmap=pl.cm.RdBu_r,extend='both')
		pl.colorbar()
		drawCoastlinesNoRivers(self.proj_.m)
		self.proj_.m.drawparallels([70,80],latmax=90)
		"""

		"""
		pl.figure(2)
                pl.contourf(self.proj_.x,self.proj_.y,s1,14,cmap=pl.cm.OrRd,extend='max')
                pl.colorbar()
                drawCoastlinesNoRivers(self.proj_.m)
                self.proj_.m.drawparallels([70,80],latmax=90)
                pl.show()
		pl.show()	
		"""
		return s1

	def getTracks(self,Years,LatRange=(70,83),LonRange=(-10,80)):
		nn = 25
		# Get trajectories in Years
		LATS,LONS = np.ones((0,len(self.ix_),nn)),np.ones((0,len(self.ix_),nn))
		dates     = []
		for year in Years:
			print year
			s         = self.d.getDataSnaps(Year=year,Season='DJF')
			dates     = dates + self.d.getDateList(Year=year,Season='DJF')
			lat,lon,p = np.rollaxis(s[:,self.ix_,:nn,0],1,0),np.rollaxis(s[:,self.ix_,:nn,1],1,0),np.rollaxis(s[:,self.ix_,:nn,2],1,0)
			LATS = np.append(LATS,lat[:,:,:],axis=0)
			LONS = np.append(LONS,lon[:,:,:],axis=0)
		# Timesteps of trajectory coordinates
		hours       = np.array([self.d.getHours(*date) for date in dates])
		hours       = np.tile(hours[:,np.newaxis,np.newaxis],(1,len(self.ix_),nn))
		dhours      = (np.arange(nn)*6)[np.newaxis,np.newaxis,:]
		hours       = hours - dhours
		# Reshape to single axes
		hours       = hours.reshape((-1,nn))
		lons,lats   = LONS.reshape((-1,nn)),LATS.reshape((-1,nn))
		print hours.shape,lons.shape,lats.shape
                # Filter out trajectories for sector
		latargmin       = np.argmin(lats,axis=1)
		indices         = np.where(lats.min(axis=1)<=LatRange[1])[0]
		latargmin       = latargmin[indices]
		hours,lons,lats = hours[indices,:],lons[indices,:],lats[indices,:]
		print hours.shape,lons.shape,lats.shape
		blatarg         = [np.where(lats[i,:]<LatRange[1])[0][0] for i in range(len(lats))]
                indices         = []
                for t in range(len(lats)):
			i0 = blatarg[t]
			i1 = latargmin[t]
			indices.append(1.*len(np.where((lons[t,i0:i1]>LonRange[0])&(lons[t,i0:i1]<LonRange[1])==True)[0])/(i1-i0+1))
		indices = np.where(np.array(indices)>0.8)[0]
		hours,lons,lats = hours[indices,:],lons[indices,:],lats[indices,:]
		print hours.shape,lons.shape,lats.shape
		hours,lons,lats = hours.reshape(-1),lons.reshape(-1),lats.reshape(-1)
		# Sort points by timestep
		indices         = np.where(lats>60)[0]
		hours,lons,lats = hours[indices],lons[indices],lats[indices]
		x,y             = self.proj_.m(lons,lats)
		Lats,points     = {},{}
		for ii in range(len(hours)):
			try:
				points[hours[ii]].append([x[ii],y[ii]])
				Lats[hours[ii]].append(lats[ii])
			except:
				points[hours[ii]] = []
				points[hours[ii]].append([x[ii],y[ii]])
                                Lats[hours[ii]] = []
                                Lats[hours[ii]].append(lats[ii])
		hours,points,lats = points.keys(),[points[i] for i in points],[Lats[i] for i in Lats]

		points_   = np.array([points[ii][jj] for ii in range(len(points)) for jj in range(len(points[ii]))])
		XS,YS     = points_[:,0],points_[:,1]
		#Nd,x0,y0  = self.density(XS,YS,xy=True)
		#self.plot(x0,y0,Nd)	

		# Interpolate trajectory coordinates to field
		PW,LT = {},{}
		for year in range(1979,2016+1,1):
			PW[year] = []
			LT[year] = []
		for ii in range(len(hours)):
			hour  = hours[ii]
			yr,mn = self.reds.getDate(hour)[0:2]
			yr    = yr + mn/11
			snap  = self.proj_(self.reds.snapshot(*self.reds.getDate(hour)),self.reds.lon,self.reds.lat)
			func  = interpolate.interp2d(self.x_,self.y_,snap, kind='linear')
			for jj in range(len(points[ii])):
				x,y = points[ii][jj]
                        	#xi,yi = np.argmin(abs(self.x_-x)),np.argmin(abs(self.y_-y))
				PW[yr].append(func(x,y)[0])
				#PW[yr].append(snap[yi,xi])
				LT[yr].append(lats[ii][jj])
		for key in PW.keys():
			if PW[key] == []:
				del PW[key]
				del LT[key]
		#PW = [PW[i][j] for i in PW for j in range(len(PW[i]))]
		#LT = [LT[i][j] for i in LT for j in range(len(LT[i]))]
		return PW,LT

	def getPoints(self,Years,LatRange,LonRange):
		if LonRange[0]<0: LonRange[0] = LonRange[0] + 360
		if LonRange[1]<0: LonRange[1] = LonRange[1] + 360
		if LonRange[0] > LonRange[1]:
			xs  = np.where((self.reds.lon<LonRange[1])|(self.reds.lon>LonRange[0])==True)[0]
		else:
			xs  = np.where((self.reds.lon<LonRange[1])&(self.reds.lon>LonRange[0])==True)[0]
		data,lats = {},{}
		for year in Years:
			print year
			n          = len(self.reds.getDateList(Year=year,Season='DJF'))
			data[year] = self.reds.getDataSnaps(Year=year,Season='DJF')[:,:,xs].reshape(-1)
			lats[year] = np.tile(self.reds.lat[np.newaxis,:,np.newaxis],(n,1,len(xs))).reshape(-1)
		return data,lats
	
	def sectors(self,Years=range(1981,2016+1,1)):
                # Bounds
                normed   = 'absolute' 
		norm     = normed=='normed'
		case     = 'tracks'

		xe,ye,H0,years,H01m,H02m = self.pdf(Years, [-180,180],[0,80],norm=norm,case=case)
		xe,ye,H1,years,H11m,H12m = self.pdf(Years,   [-30,80],[0,80],norm=norm,case=case)
		xe,ye,H2,years,H21m,H22m = self.pdf(Years,   [80,180],[0,80],norm=norm,case=case)
		xe,ye,H3,years,H31m,H32m = self.pdf(Years,[-180,-120],[0,80],norm=norm,case=case)
		xe,ye,H4,years,H41m,H42m = self.pdf(Years,  [-120,-30],[0,80],norm=norm,case=case)
		#Hr             = H0 - (H1+H2+H3+H4)
		#Hr1m = ((Hr[0:18].sum(axis=0))*ye[:,np.newaxis]).sum(axis=0)/Hr[0:18,:,:].sum(axis=0).sum(axis=0)
                #Hr2m = (( Hr[18:].sum(axis=0))*ye[:,np.newaxis]).sum(axis=0)/Hr[18:,:,:].sum(axis=0).sum(axis=0)

		trend0,pval0 = self.getTrend(years,H0)
		trend1,pval1 = self.getTrend(years,H1)
		trend2,pval2 = self.getTrend(years,H2)
		trend3,pval3 = self.getTrend(years,H3)
		trend4,pval4 = self.getTrend(years,H4)
		#trendr,pvalr = self.getTrend(years,Hr)

		self.plotPDF(xe,ye,trend0,pval0,H01m,H02m, [-180,180],norm,normed,case)
		self.plotPDF(xe,ye,trend1,pval1,H11m,H12m,   [-30,80],norm,normed,case)
		self.plotPDF(xe,ye,trend2,pval2,H21m,H22m,   [80,180],norm,normed,case)
		self.plotPDF(xe,ye,trend3,pval3,H31m,H32m,[-180,-120],norm,normed,case)
		self.plotPDF(xe,ye,trend4,pval4,H41m,H42m, [-120,-30],norm,normed,case)
		#self.plotPDF(xe,ye,trendr,pvalr,Hr1m,Hr2m,'RES',      (0,0],norm,normed)

	def pdf(self,Years,LonRange,LatRange,norm,case='tracks'):
		#xedges = np.linspace(60,90,17)
		xedges  = np.arange(60,90+1.5,1.5)
		yedges = np.linspace(0,10,24)
		if case == 'tracks': PW1,lats1 = self.getTracks(Years,LatRange,LonRange)
		if case == 'points': PW1,lats1 = self.getPoints(Years,LatRange,LonRange)
		HH1           = []
		for key in PW1.keys():
			xe,ye,H1 = self.bivarPDF(lats1[key],PW1[key],xedges,yedges,norm=norm)
			HH1.append(H1)
		HH1 = np.array(HH1)
		HH1 = np.rollaxis(HH1,2,1)
                # Means
		i0  = 18
                H1m = ((HH1[0:i0].sum(axis=0))*ye[:,np.newaxis]).sum(axis=0)/(HH1[0:i0,:,:].sum(axis=0).sum(axis=0))
                H2m = (( HH1[i0:].sum(axis=0))*ye[:,np.newaxis]).sum(axis=0)/(HH1[i0:,:,:].sum(axis=0).sum(axis=0))
		return xe,ye,HH1,PW1.keys(),H1m,H2m

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

	def plotPDF(self,xe,ye,trend,pval,H1m,H2m,LonRange,norm,normed,case):
		pl.figure(1)
		#cseq = np.arange(-105,105+15,15)
		cseq  = 15
                xxe,yye = np.meshgrid(xe,ye)
                mx,my   = np.ma.masked_where(pval>0.1,xxe),np.ma.masked_where(pval>0.1,yye)
		cf      = pl.contourf(xe,ye,trend,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cbar    = pl.colorbar(cf)
		cbar.set_label('Frequency [decade$^{-1}$]')
		pl.plot(mx,my,'k.',markersize=5,alpha=0.35)
		pl.plot(xe,H1m,'k--',linewidth=1.2,alpha=0.8)
		pl.plot(xe,H2m,'r--',linewidth=1.2,alpha=0.8)
		pl.xlabel('Latitude')
		pl.ylabel('%s [%s]' % (self.reds.long_name,self.reds.units))	
		pl.savefig('figs/trackfield/dH.%s-%sE.%s.%s.pdf' % (LonRange[0],LonRange[1],normed,case),format='pdf')
		pl.close()

        def bivarPDF(self,x,y,xedges,yedges,norm=True):
                H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=[[xedges[0],xedges[-1]],[yedges[0],yedges[-1]]],normed=norm,weights=None)
                xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
                return xedges,yedges,H

	def tracksInBox(self,Years,frac=0.6):
		PW1,PW2 = [],[]
		PW3,PW4 = [],[]
		PW5,PW6 = [],[]
		PW0     = np.ones((0,len(self.ix_)))
		for Year in Years:
			print Year
			pw        = self.proj(self.reds.getDataSnaps(Year=Year,Season='DJF'),self.reds.lon,self.reds.lat).reshape((-1,100))
			pw        = pw[:,self.ix_]
			PW0       = np.append(PW0,pw,axis=0)
			s         = self.d.getDataSnaps(Year=Year,Season='DJF')
                	lat,lon,p = s[:,self.ix_,:28,0],s[:,self.ix_,:28,1],s[:,self.ix_,:28,2]
			lat,lon,p = np.rollaxis(lat,1,0),np.rollaxis(lon,1,0),np.rollaxis(p,1,0)	# Axes rotated when extracting slice for some reason...
                	for t in range(len(lat)):
                	        for i in range(len(lat[t])):
                	                lats,lons = lat[t,i,:],lon[t,i,:]
                                        if lats.min()<=70:                              # Only tracks which reached 70N
						i0 = np.where(lats<=80)[0][0]
						i1 = np.where(lats<=70)[0][0]+1
						if ((lons[i0:i1]>-30)&(lons[i0:i1]<110)==True).all():
							PW1.append(pw[t,i])
                                                elif ((lons[i0:i1]>80)&(lons[i0:i1]<180)==True).all():
                                                        PW2.append(pw[t,i])
                                                elif ((lons[i0:i1]>-180)&(lons[i0:i1]<-30)==True).all():
                                                        PW3.append(pw[t,i])
                                                elif ((lons[i0:i1]>-30)&(lons[i0:i1]<20)==True).all():
                                                        PW4.append(pw[t,i])
						else:
							PW5.append(pw[t,i])
						""""
						if 1.*np.sum((lons[i0:]>20)&(lons[i0:]<80).astype(int))/len(lons[i0:])      > frac:
							PW1.append(pw[t,i])
                                                elif 1.*np.sum((lons[i0:]>80)&(lons[i0:]<180).astype(int))/len(lons[i0:])   > frac:
                                                        PW2.append(pw[t,i])
                                                elif 1.*np.sum((lons[i0:]>-180)&(lons[i0:]<-30).astype(int))/len(lons[i0:]) > frac:
                                                        PW3.append(pw[t,i])
                                                elif 1.*np.sum((lons[i0:]>-30)&(lons[i0:]<20).astype(int))/len(lons[i0:])   > frac:
                                                        PW4.append(pw[t,i])
						else:
							PW5.append(pw[t,i])
						"""
					else:
						PW6.append(pw[t,i])
		return PW0.reshape(-1),PW1,PW2,PW3,PW4,PW5,PW6

	def plotTest(self):
		Full,Bar,Eur,Pac,Atl,Res0,Res1 = self.tracksInBox(Years=range(1981,1998+1,1))

		edges = np.arange(0,9+0.2,0.2)
		h0,e0 = np.histogram(Full,edges,normed=False)
		h1,e1 = np.histogram(Bar ,edges,normed=False)
		h2,e2 = np.histogram(Eur ,edges,normed=False)
		h3,e3 = np.histogram(Pac ,edges,normed=False)
		h4,e4 = np.histogram(Atl ,edges,normed=False)
		h5,e5 = np.histogram(Res0,edges,normed=False)
		h6,e6 = np.histogram(Res1,edges,normed=False)

		pl.plot(e0[0:-1],h0,'k'  ,linewidth=1.75,alpha=0.85,label='All points')
		pl.plot(e1[0:-1],h1,'r'  ,linewidth=1.25,alpha=0.65,label='Barents and Kara')
		pl.plot(e2[0:-1],h2,'b'  ,linewidth=1.25,alpha=0.65,label='Eurasia')
                pl.plot(e3[0:-1],h3,'g'  ,linewidth=1.25,alpha=0.65,label='Pacific and North America')
		pl.plot(e4[0:-1],h4,'y'  ,linewidth=1.25,alpha=0.65,label='Atlantic')
		pl.plot(e5[0:-1],h5,'b--',linewidth=1.00,alpha=0.45,label='Residual0')
		pl.plot(e6[0:-1],h6,'k--',linewidth=1.00,alpha=0.45,label='Tracks north of 70N')

		lg = pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)

		pl.xlabel('%s [%s]' % (self.reds.long_name,self.reds.units))
		pl.ylabel('Frequency')
		pl.yscale('linear')
		pl.show()

	def fieldTracks(self,LATS,LONS):
		PW = []
		xs,ys = self.proj_.m(LONS[:,:],LATS[:,:])
		for i in range(len(xs)):
			#weights = self.trackWeights(xs[i,:],ys[i,:])
			pws0    = []
			for j in range(len(xs[i])):
				pws0.append(self.interpfunc(xs[i,j],ys[i,j]))
			pws0 = np.array(pws0).squeeze()
			PW.append((pws0*self.weights).sum()/(self.weights.sum()))
			#PW.append(pws0[0:20].mean())
		return PW

	def trackWeights(self,xs,ys):
		d = np.sqrt((xs - xs[0])**2 + (ys - ys[0])**2)
		#d = np.diff(d)

	def density(self,LONS,LATS,xy=False):
		# if xy == True then LONS and LATS are treated as XS and YS, respectively (using proj_)
		if len(LATS.shape) > 1: LATS,LONS = LATS.reshape(-1),LONS.reshape(-1)
		if not xy:
			xs,ys = self.proj_.m(LONS,LATS)
		else:
			xs,ys = LONS,LATS
		N = np.zeros((self.proj_.nx,self.proj_.ny))
		for i in range(len(xs)):
		        xi,yi    = np.argmin(abs(self.x_-xs[i])),np.argmin(abs(self.y_-ys[i]))
		        N[yi,xi] = N[yi,xi] + 1
		N,x,y    = N[1:-1,1:-1],self.x_[1:-1],self.y_[1:-1]
		#N,xx,yy = d.interp2d(N,x,y,6,kind='cubic')
		return N,x,y

	def plot(self,x,y,N):
		pl.figure(1)
		cf    = pl.contourf(x,y,N,15,cmap=pl.cm.OrRd,extend='max')
		cbar  = pl.colorbar(cf)
		cbar.set_label(r'Number density {400$\times$400 km$^{2}$}$^{-1}$')
		drawCoastlinesNoRivers(self.proj_.m)
		self.proj_.m.drawparallels([70,80],latmax=90)
		pl.xlim(self.x_[1],self.x_[-2])
		pl.ylim(self.y_[1],self.y_[-2])
		pl.show()

if __name__ == "__main__":

	st = SelectTracks()
	st.sectors()

	"""
	Years1        = range(1981,2016+1,1)
	Years2        = range(1999,2012+1,1)
	LATS1,LONS1   = st.getTracks(Years1)
	pw = st.pdf(Years1,Years2)
	PW = st.fieldTracks(LATS1,LONS1)

        h1,e1  = np.histogram(pw,50,normed=False)
	h2,e2  = np.histogram(PW,50,normed=False)

	pl.plot(e1[0:-1],h1,'b',linewidth=1.25,alpha=0.65)
	pl.plot(e2[0:-1],h2,'r',linewidth=1.25,alpha=0.65)
	pl.grid()
	pl.ylabel('Frequency')
	pl.xlabel('%s [%s]' % (st.reds.long_name,st.reds.units))
	pl.show()

	pl.show()
	"""

