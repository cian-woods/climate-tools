from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from mpl_toolkits import basemap
from toPick import *
from UnPickle import *
from scipy import interpolate
from scipy import ndimage as nd
from scipy import stats
from datetime import datetime
from matplotlib.colors import LinearSegmentedColormap

import numpy as np
import matplotlib.pyplot as pl
import sys
import heapq
import subprocess

class FluxFile:

	def __init__(	self,
			Source    = 'CMCC-CESM',
			blat      =  	     80,
			YearRange = (1981,2005),
			LevRange  =    (0,1000),
			Season    =       'DJF',
			makefile  =       True		):

		# Attributes
		self.Source    = Source
		self.nlevs     = 100
		self.plevs     = np.linspace(LevRange[0],LevRange[1],self.nlevs)
		self.blat      = blat
		self.LatRange  = (blat-2,blat+3)
		self.YearRange = YearRange
		self.LevRange  = LevRange
		self.years     = range(YearRange[0],YearRange[1]+1,1)
		self.Season    = Season
		self.Re        = 6371000
		self.xres      = 2*np.pi*np.cos(blat*np.pi/180)*self.Re/360.

		# DataServers
		if makefile == True:
			if Source[0:3] == 'ERA':
				self.vds = reDataServer(Field='V',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source=Source)
				self.qds = reDataServer(Field='q',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source=Source)
				self.dailymean,self.step = False,6
			elif Source[0:4] == 'NCEP':
				self.vds = reDataServer(Field='vplev',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source=Source)
				self.qds = reDataServer(Field='shum',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source=Source)
				self.dailymean,self.step = False,6
                        elif Source == 'MERRA2':
                                self.vds = reDataServer(Field='V',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source='MERRA2')
                                self.qds = reDataServer(Field='QV',LevType='plev',LevRange=(0,1000),LatRange=self.LatRange,Source='MERRA2')
				self.dailymean,self.step = False,24
			elif Source[0:3] == 'CAM':
				self.vds = reDataServer(Field='V',Source=Source,LevRange=(0,1000),LatRange=self.LatRange)
				self.qds = reDataServer(Field='q',Source=Source,LevRange=(0,1000),LatRange=self.LatRange)
				self.dailymean,self.step = False,6
                        else:
				ExpType = 'rcp85'
				#if YearRange[1]  > 2005: ExpType='rcp85'
				#if YearRange[1] <= 2005: ExpType='historical'
                                self.vds = cmipDataServer(Field='va' ,LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=self.LatRange,ExpType=ExpType)
                                self.qds = cmipDataServer(Field='hus',LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=self.LatRange,ExpType=ExpType)
				self.dailymean,self.step = True,24

                # Mass of girdbox on each level
		self.dP = []
		dP = np.diff(self.plevs)/2.
		self.dP.append(dP[0])
		for i in range(len(dP)-1):
			self.dP.append(dP[i]+dP[i+1])
		self.dP.append(dP[-1])
		self.dP = np.array(self.dP)
		self.dM = 100.*self.dP/9.80665                             # Mass per unit area between pressure levels of data (kg/m^2)

	def fill(self,data,invalid=None):
		if invalid is None: invalid = np.isnan(data)
		ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
		return data[tuple(ind)]

	def interpolateND(self,field,xold,xnew,axis):
		# field in (time)x(lev)x(lon)
		#f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
		f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=False,fill_value='extrapolate')
		field = f(xnew)	
		return field

	def pwError(self):

		vqERA,datesERA = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/ERAInt.vapor.1981-2005.30-1000hPa.%sN.DJF.p' % (self.blat))
		vqERA          = np.array(vqERA).mean(axis=0)
		vqlonsERA      = np.arange(0,360,360./len(vqERA))
		pwERA          = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/ERAInt.pw.1981-2005.%sN.DJF.p' % (self.blat))
		pwlonsERA      = np.arange(0,360,360./len(pwERA))

		vq,dates = self.flux(case='vapor')
		vq       = vq.mean(axis=0)

		pw = np.array([self.interpolateND(self.pds.getDataSnaps(Year=year,Season=self.Season),self.pds.lat,self.blat,axis=1).mean(axis=0)\
				for year in self.years]).mean(axis=0)	

		pl.plot(self.vds.lon,vq,'r',linewidth=1,alpha=0.5)
		pl.plot(self.pds.lon,pw,'b',linewidth=1.5)
		pl.plot(vqlonsERA,vqERA,'grey',linewidth=1,alpha=0.5)
		pl.plot(pwlonsERA,pwERA,'k',linewidth=1.5)

		pl.grid()
		pl.ylim(0,9)
		pl.xlim(0,360)
		pl.title(Source)
		pl.xlabel('Longitude')
		pl.ylabel('Precipitable water [kg m$^{-2}$]')
		pl.show()

	def flux(self,case='moist'):

		fname = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.%s-%s.%s-%shPa.%sN.%s.p' % (self.blat,self.Source,self.Source,case,self.YearRange[0],\
									      self.YearRange[1],self.LevRange[0],self.LevRange[1],self.blat,self.Season)
		makedir = "mkdir -p /mnt/climstorage/cian/scripts/newfluxfiles/%s/%s" % (self.blat,self.Source)
		print makedir
		print fname
		os.system(makedir)
		if os.path.isfile(fname) == False:
			print 'Calculating fluxes for year %s ...' % (self.years[0])

			# First year (fill missing with nearest and interpolate to blat and self.plevs)
			if case == 'moist':
				v = self.vds.getDataSnaps(Year=self.years[0],Season=self.Season,dailymean=self.dailymean,step=self.step)
				q = self.qds.getDataSnaps(Year=self.years[0],Season=self.Season,dailymean=self.dailymean,step=self.step)	
				if (type(q)==np.ma.core.MaskedArray) and (q.mask.shape == q.data.shape):
					v  = self.fill(v.data,invalid=v.mask)
					q  = self.fill(q.data,invalid=q.mask)
				v    = self.interpolateND(v,self.vds.lat,self.blat,axis=2)
				v    = self.interpolateND(v,self.vds.lev,self.plevs,axis=1)
				q    = self.interpolateND(q,self.qds.lat,self.blat,axis=2)
				q    = self.interpolateND(q,self.qds.lev,self.plevs,axis=1)
			elif case == 'mass':
				v = self.vds.getDataSnaps(Year=self.years[0],Season=self.Season,dailymean=self.dailymean,step=self.step)
				if (type(v)==np.ma.core.MaskedArray) and (v.mask.shape == v.data.shape):
					v  = self.fill(v.data,invalid=v.mask)
				v    = self.interpolateND(v,self.vds.lat,self.blat,axis=2)
				v    = self.interpolateND(v,self.vds.lev,self.plevs,axis=1)
				q = 1.
			elif case == 'vapor':
				v    = 1.
				q    = self.qds.getDataSnaps(Year=self.years[0],Season=self.Season,dailymean=self.dailymean,step=self.step)
				if (type(q)==np.ma.core.MaskedArray) and (q.mask.shape == q.data.shape):
					q  = self.fill(q.data,invalid=q.mask)
				q    = self.interpolateND(q,self.qds.lat,self.blat,axis=2)
				q    = self.interpolateND(q,self.qds.lev,self.plevs,axis=1)

	                # Fix overflow problem
	                v,q = np.float64(v),np.float64(q)

			# Vertically integrate (time)x(lon) now
			vq0    = ((v*q)*self.dM[np.newaxis,:,np.newaxis]).sum(axis=1)
			dates0 = self.vds.getDateList(Year=self.years[0],Season=self.Season,step=self.step)

			for year in self.years[1:]:
				print 'Calculating fluxes for year %s ...' % (year)

				# v and q data
				if case == 'moist':
					v = self.vds.getDataSnaps(Year=year,Season=self.Season,dailymean=self.dailymean,step=self.step)
					q = self.qds.getDataSnaps(Year=year,Season=self.Season,dailymean=self.dailymean,step=self.step)
					if (type(q)==np.ma.core.MaskedArray) and (q.mask.shape == q.data.shape):
						v  = self.fill(v.data,invalid=v.mask)
						q  = self.fill(q.data,invalid=q.mask)
					v    = self.interpolateND(v,self.vds.lat,self.blat,axis=2)
					v    = self.interpolateND(v,self.vds.lev,self.plevs,axis=1)
					q    = self.interpolateND(q,self.qds.lat,self.blat,axis=2)
					q    = self.interpolateND(q,self.qds.lev,self.plevs,axis=1)
				elif case == 'mass':
					v = self.vds.getDataSnaps(Year=year,Season=self.Season,dailymean=self.dailymean,step=self.step)
					if (type(v)==np.ma.core.MaskedArray) and (v.mask.shape == v.data.shape):
						v  = self.fill(v.data,invalid=v.mask)
					v    = self.interpolateND(v,self.vds.lat,self.blat,axis=2)
					v    = self.interpolateND(v,self.vds.lev,self.plevs,axis=1)
					q = 1.
				elif case == 'vapor':
					v = 1.
					q = self.qds.getDataSnaps(Year=year,Season=self.Season,dailymean=self.dailymean,step=self.step)
					if (type(q)==np.ma.core.MaskedArray) and (q.mask.shape == q.data.shape):
                                        	q  = self.fill(q.data,invalid=q.mask)
                                	q    = self.interpolateND(q,self.qds.lat,self.blat,axis=2)
					q    = self.interpolateND(q,self.qds.lev,self.plevs,axis=1)

				# Datelist
				dates = self.vds.getDateList(Year=year,Season=self.Season,step=self.step)

		                # Fix overflow problem
		                v,q = np.float64(v),np.float64(q)

				# Vertically integrate (time)x(lon) now [kg s**-1 m**-1]
				vq = ((v*q)*self.dM[np.newaxis,:,np.newaxis]).sum(axis=1)

				# Append to vq0
				vq0    = np.append(vq0,vq,axis=0)
				dates0 = dates0 + dates

			# Convert to Tg day**-1 deg**-1
			if (case == 'mass') or (case == 'moist'):
				vq0 = vq0*24*60*60*self.xres/1e09
			toPick([vq0,dates0],fname)
		else:
			vq0,dates0 = unpick(fname)
			vq0 = np.array(vq0)
		return vq0,dates0

	def pdf(self,vq,edges,integrate=True,Field0=[]):

		if Field0 == []: Field0 = vq[:,:]
		n0,n1 = vq.shape
		n     = len(edges)
		H     = np.zeros((n-1,n1))
		N     = np.zeros((n-1,n1))
		for i in range(n0):
			for j in range(n1):
				for k in range(n-1):
					if edges[k] < vq[i][j] <= edges[k+1]:
						H[k,j] = H[k,j] + Field0[i][j]
						N[k,j] = N[k,j] + 1
						break
		# Shiftgrid
#		lons   = np.arange(0,359+1,1)
#		H,lons = basemap.shiftgrid(330,H,lons,start=True)
		# Integrate
		if integrate == True:
			n0,n1 = H.shape
			H     = H.reshape((n0,n1/10,-1)).sum(axis=2)
			N     = N.reshape((n0,n1/10,-1)).sum(axis=2)
		n0,n1 = H.shape
		# Rotate H to -180 --> 180
		lons  = np.arange(0,360,360./n1)
		x     = np.argmin((lons-(360-80))**2)+1
		H     = np.append(H[:,x:],H[:,0:x],axis=1)
		N     = np.append(N[:,x:],N[:,0:x],axis=1)
		Hmean = (np.ma.masked_where(N==0,H)/np.ma.masked_where(N==0,N)).data	
		return H[::-1,:],Hmean[::-1,:]

	def pdf_new(self,vq,edges,integrate=True,Field0=[]):

		if Field0 == []: Field0 = vq[:,:]
		n0,n1 = vq.shape

		H = np.zeros((len(edges)-1,n1))
		N = np.zeros((len(edges)-1,n1))
		for i in range(n1):
			for j in range(len(edges)-1):
				xs     = np.where((edges[j]<=vq[:,i])&(vq[:,i]<edges[j+1])==True)[0]
				H[j,i] = vq[xs,i].sum()
				N[j,i] = N[j,i] + 1
		if integrate == True:
			n0,n1 = H.shape
			H     = H.reshape((n0,n1/10,-1)).sum(axis=2)*self.xres*(60*60*24)/1e09
			N     = N.reshape((n0,n1/10,-1)).sum(axis=2)
		# Rotate H to -180 --> 180
		n0,n1 = H.shape
		lons  = np.arange(0,360,360./n1)
		x     = np.argmin((lons-(360-80))**2)+1
		H     = np.append(H[:,x:],H[:,0:x],axis=1)
		N     = np.append(N[:,x:],N[:,0:x],axis=1)
		Hmean = (np.ma.masked_where(N==0,H)/np.ma.masked_where(N==0,N)).data
		return H[::-1,:],Hmean[::-1,:]

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

	def plot(self,H,edges,vmin=None,vmax=None,xlabel='',ylabel='',Source='',sname=None,clabel='',cmap=pl.cm.RdBu_r,show=True,Hstd=None,stip=None,pthresh=240,sf=1,bin_integral=False,stationary=None):

                fig,ax1 = pl.subplots(num=1)
		n0,n1 = H.shape
		clons = np.arange(-180 + 360./(2*n1), 180, 360./n1)
		ii0   = np.argmin(edges[::-1]**2)
	        s1    = round(H[0:ii0,:].sum(),0)
	        s2    = round(H[ii0:,:].sum(),0)	
		sign1,sign2 = '',''
		if s1 > 0: sign1 = '+'
		if s2 > 0: sign2 = '+'
	        H = np.ma.masked_where(H==0,H)
		# Plot
	        cf = ax1.imshow(H,aspect='auto',extent=[-180,180,edges[0],edges[-1]],interpolation='nearest',cmap=cmap,vmin=vmin,vmax=vmax)
		if Hstd != None:
			cseq      = 4#np.arange(10,40+10,10)*sf
			bins,lons = np.linspace(edges[0],edges[-1],n0),np.linspace(-180,180,n1)
			cl        = ax1.contour(lons,bins,Hstd[::-1,:],cseq,colors='k',alpha=0.5)
			ax1.clabel(cl,fmt='%2.0f',colors='k',fontsize=10)
		if stip != None:
			ax1.plot(stip[0],stip[1],'k.',alpha=0.25)
		if pthresh != None: ax1.plot([-180,180],[pthresh,pthresh],'k--',linewidth=1.5,alpha=0.6)
	        cbar  = fig.colorbar(cf,extend='both')
		cbar.set_label(clabel)
		
		if bin_integral == True:
			ax2     = ax1.twinx()
			ax2.plot(clons,1e12*H.sum(axis=0)/(10*self.xres*60*60*25*1e03),'k.-',linewidth=0.75,alpha=0.35,label='$\Delta{\overline{vq}}$')
			ax2.set_ylabel('Moisture flux bias [kg m$^{-1}$ s$^{-1}$]')
		if stationary != None:	
			ax1.plot(np.arange(-175,175+10,10),stationary*100,'r.-',linewidth=0.75,alpha=0.35,label='$\Delta{\overline{v}_{500hPa}}$')
			ax1.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)

		ax = pl.gca()
		ax.text(0.05, 0.95, '%s%s %s' % (sign1,s1,'Tg day$^{-1}$'), transform=ax.transAxes, fontsize=11, verticalalignment='top')
		ax.text(0.05, 0.05, '%s%s %s' % (sign2,s2,'Tg day$^{-1}$'), transform=ax.transAxes, fontsize=11, verticalalignment='top')

		ax1.set_xlabel(xlabel)
		ax1.set_ylabel(ylabel)
		pl.title('%s: %s %s' % (Source,s1+s2,'Tg day$^{-1}$'))
		ax1.set_yticks(edges[::4])	
		ax1.set_ylim(edges[0],edges[-1])
		ax1.set_xticks([-160,-100,-40,20,80,140],['-60','0','60','120','180','-120'])
		ax1.set_xlim(-180,180)
		if sname != None: pl.savefig(sname,format='pdf')
		if show == True: pl.show()
		pl.close()

	def stats(self):
		LevRange = (30,1000)
		case     = 'moist'
		Models   = [g[14:] for g in glob.glob('../historical/*')]
		for Model in Models:
			fname    = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2005.%s-%shPa.%sN.DJF.p' % (self.blat,Model,Model,case,LevRange[0],LevRange[1],self.blat)
			vq,dates = unpick(fname)
			vq       = self.interpolate(vq)
			l        = stats.scoreatpercentile(vq,90,axis=0)

        def fluxToYears(self,vq,dates):
                vqt = {}
                dt  = {}
                for i in range(len(vq)):
                        year,month,day,hour = dates[i]
                        if month > 8: year = year + 1
                        try:
                                vqt[year].append(vq[i])
                                dt[year].append(dates[i])
                        except:
                                vqt[year] = []
                                dt[year]  = []
                                vqt[year].append(vq[i])
                                dt[year].append(dates[i])
                return vqt.values(),dt.values(),vqt.keys()

	def CMIP5_Reanalyses_Biases(self):
		# Parameters
		show        = False
		integrate   =  True
		case        = 'moist'
                ylabel      = 'Moisture flux intensity [Tg day$^{-1}$ deg$^{-1}$]'
                clabel      = 'Dailymean transport [Tg day$^{-1}$]'
                vminf,vmaxf =  0,200
                vmind,vmaxd = -100,100
                pthresh     =  240
                cmap        =  pl.cm.RdBu_r
                colors      =  cmap(np.linspace(0.5, 1, cmap.N/2))
                cmap        =  LinearSegmentedColormap.from_list('Upper Half', colors)
                cint        =  75
                edges       =  np.arange(0,1500+cint,cint)

		# Compute biases with respect to each reanalysis, as well as inter-reanalysis biases
		YearR  = [[1981,2005],[1981,2005],[1981,2005],[1981,2005],[1975,2001]]
		Reanal = ['ERAInt','NCEP1','NCEP2','MERRA2','ERA40']
		Models = [x[27:] for x in glob.glob('/qnap/cian/cmip/historical/*')]
		Hm,Hr  = [],[]
	        # Model mean PDF
		for ii in range(len(Models)):
                	print '%s of %s ...' % (ii+1,len(Models))
                	Source     = Models[ii]
                	fmname     = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2005.300-1000hPa.%sN.DJF.p' % (Source,Source,case)
                	print fmname
                	vqm,datesm = unpick(fmname)
                	vqm        = self.interpolate(vqm)
                	hm         = self.pdf(vqm,edges,integrate=integrate)/len(datesm)
                	Hm.append(hm)
	        Hm = np.array(Hm)
                # Reanalysis mean PDF
                for ii in range(len(Reanal)):
                        print '%s of %s ...' % (ii+1,len(Reanal))
                        Source     = Reanal[ii]
                        frname     = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.%s-%s.300-1000hPa.%sN.DJF.p' % (Source,Source,case,YearR[ii][0],YearR[ii][1])
                        print frname
                        vqr,datesr = unpick(frname)
                        vqr        = self.interpolate(vqr)
                        hr         = self.pdf(vqr,edges,integrate=integrate)/len(datesr)
                        Hr.append(hr)
                Hr = np.array(Hr)
		# Arry containing mean biases wrt each reanalysis
        	DH = np.array([Hm-Hr[ii][np.newaxis,:,:] for ii in range(len(Reanal))])
                Hmm = Hm.mean(axis=0)
                Hrm = Hr.mean(axis=0)
		Hrs = Hr.std(axis=0)
		Hms = Hm.std(axis=0)
		for ii in range(len(DH)):
			dH = DH[ii]
        		# Stipling
        		n0,n1,n2    = dH.shape
        		bins,lons   = np.linspace(edges[0]+cint/2.,edges[-1]-cint/2.,n1),np.linspace(-180+5,180-5,n2)
        		ax_x,ax_y   = np.meshgrid(lons,bins)
        		ax_x,ax_y   = ax_x[::-1,:],ax_y[::-1,:]
        		stipx,stipy = stipling(dH,xx=ax_x,yy=ax_y,thresh=0.8)
        		# Means and standard deviation
        		Hs  = dH.std(axis=0)
        		dH  = dH.mean(axis=0)
        		# Plot
        		ff.plot(dH,edges,vmin=vmind/2.,vmax=vmaxd/2.,xlabel='Longitude',ylabel=ylabel,show=show,\
        		        Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.%s.diff.pdf' % (case,'all',Reanal[ii]),clabel=clabel,\
				cmap=pl.cm.RdBu_r,Hstd=Hs,stip=(stipx,stipy),pthresh=pthresh)
			ff.plot(Hr[ii]-Hrm,edges,vmin=vmind/2.,vmax=vmaxd/2.,xlabel='Longitude',ylabel=ylabel,show=show,\
                                Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.reanal.%s.diff.pdf' % (case,'all',Reanal[ii]),clabel=clabel,\
                                cmap=pl.cm.RdBu_r,pthresh=pthresh)
                	ff.plot(Hr[ii],edges,vmin=vminf/2.,vmax=vmaxf/2.,xlabel='Longitude',ylabel=ylabel,show=show,\
                 	        Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,Reanal[ii]),clabel=clabel,cmap=cmap,pthresh=pthresh)
		# Mean Model and Reanalysis plots
		ff.plot(Hmm,edges,vmin=vminf,vmax=vmaxf,xlabel='Longitude',ylabel=ylabel,show=show,\
                        Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.model.full.pdf' % (case,'all'),clabel=clabel,cmap=cmap,Hstd=Hms,pthresh=pthresh,sf=1.)
                ff.plot(Hrm,edges,vmin=vminf,vmax=vmaxf,xlabel='Longitude',ylabel=ylabel,show=show,\
                        Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.reanal.full.pdf' % (case,'all'),clabel=clabel,cmap=cmap,Hstd=Hrs,pthresh=pthresh,sf=0.5)
      
if __name__ == "__main__":

	from stipling import stipling
	import glob


	# Single model
	Source    = str(sys.argv[1])
	YearRange = (int(sys.argv[2]),int(sys.argv[3]))
	LevRange  = (int(sys.argv[4]),int(sys.argv[5]))
	case      = str(sys.argv[6])
	blat      = int(sys.argv[7])
	Season    = str(sys.argv[8])
	ff        = FluxFile(Source=Source,YearRange=YearRange,LevRange=LevRange,blat=blat,makefile=True,Season=Season)
	vq,dates  = ff.flux(case=case)
	print vq.shape

	sys.exit()


#	ff = FluxFile()
#	ff.CMIP5_Reanalyses_Biases()

	# All models
	case      = 'moist'
	blat      = 70
	YearRange = (1981,2016)
	show      = False
	integrate = True
	makefile  = False
        # Conversion factor; Tg deg**-1 day**-1 -> kg m**-1 s**-1
        con_fac = 1e9/((2*np.pi*6370e03*np.cos(np.pi*blat/180.)/360)*60*60*24)

	if case == 'moist':
		LevRange     = (0,1000)
		ylabel       = 'Meridional moisture flux intensity [kg m$^{-1}$ s$^{-1}$]'	#'Moisture flux intensity [Tg day$^{-1}$ deg$^{-1}$]'
		clabel       = 'Dailymean transport [Tg day$^{-1}$]'
		vminf,vmaxf  = -400,400
		vmind,vmaxd  = -200,200
		pthresh      =  0
		cmap         =  pl.cm.RdBu_r
		cint         =  25						# 75
		edges        =  np.arange(-250,650+cint,cint)			# np.arange(-900,2100+cint,cint)
		scale_factor =  1
	if case == 'vapor':
		LevRange     = (600,1000)
		clabel       = 'Precipitable water [kg m$^{-2}$]'
		ylabel       = 'Precipitable water [kg m$^{-2}$]'
		vminf,vmaxf  =  0,8
		vmind,vmaxd  = -4,4
		pthresh      =  None
		cmap         =  pl.cm.OrRd
		cint         =  0.75
		edges        =  np.arange(0,20+cint,cint)
		scale_factor = 1
	if case == 'mass':
		scale_factor = 1e03
		LevRange     = (30,1000)
		ylabel       = 'Vertically integrated mass flux [10$^{3}$ Tg day$^{-1}$ deg$^{-1}$]'
		clabel       = 'Dailymean transport [10$^{3}$ Tg day$^{-1}$]'
		vminf,vmaxf  =  -200,200
		vmind,vmaxd  =  -200,200
		pthresh      =  None
		cmap         =  pl.cm.RdBu_r
		#colors      =  cmap(np.linspace(0.5, 1, cmap.N/2))
		#cmap        =  LinearSegmentedColormap.from_list('Upper Half', colors)
		cint         =  100
		edges        =  np.arange(-1500,1500+cint,cint)

	"""
	# Single PDF
	Source1    = 'CAM4xCO2'
	Source2    = 'CAM8xCO2'
	Source3    = 'CAM16xCO2'
	ff         = FluxFile(Source=Source1,YearRange=YearRange,LevRange=LevRange,blat=blat,makefile=False,Season='NDJFM')
	# Open data
	vq1,dates1 = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1901-1919.0-1000hPa.%sN.NDJFM.p' % (blat,Source1,Source1,case))
	vq2,dates2 = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1901-1919.0-1000hPa.%sN.NDJFM.p' % (blat,Source2,Source2,case))
	vq3,dates3 = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1901-1919.0-1000hPa.%sN.NDJFM.p' % (blat,Source3,Source3,case))
	# Compute PDFs
	vq1         = ff.interpolate(vq1)/scale_factor
	vq2         = ff.interpolate(vq2)/scale_factor
	vq3         = ff.interpolate(vq3)/scale_factor
	hh1,hhmean1 = ff.pdf_new(vq1,edges,integrate=integrate,Field0=[])
	hh2,hhmean2 = ff.pdf_new(vq2,edges,integrate=integrate,Field0=[])
	hh3,hhmean3 = ff.pdf_new(vq3,edges,integrate=integrate,Field0=[])
	hh1,hh2,hh3 = hh1/len(dates1),hh2/len(dates2),hh3/len(dates3)
	# Plot full
        ff.plot(hh1,edges,vmin=vmind,vmax=vmaxd,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source=Source1,sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,Source1),clabel=clabel,cmap=cmap,pthresh=pthresh)
        ff.plot(hh2,edges,vmin=vmind,vmax=vmaxd,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source=Source2,sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,Source2),clabel=clabel,cmap=cmap,pthresh=pthresh)
        ff.plot(hh3,edges,vmin=vmind,vmax=vmaxd,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source=Source3,sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,Source3),clabel=clabel,cmap=cmap,pthresh=pthresh)
	# Plot differences
        ff.plot(hh2-hh1,edges,vmin=vminf/2.,vmax=vmaxf/2.,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='%s-%s' % (Source2,Source1),sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.diff.pdf' % (case,Source2),clabel=clabel,cmap=cmap,pthresh=pthresh)
        ff.plot(hh3-hh2,edges,vmin=vminf/2.,vmax=vmaxf/2.,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='%s-%s' % (Source3,Source2),sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.diff.pdf' % (case,Source3),clabel=clabel,cmap=cmap,pthresh=pthresh)
	"""

	"""
	# Reanalyses biases
	Modele,Modeln = 'ERAInt','ERA40'
	ff = FluxFile(Source='ERAInt',YearRange=(1981,2005),LevRange=LevRange,blat=blat,makefile=makefile)
	ename = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1980-2016.%s-%shPa.%sN.NDJFM.p' % (blat,Modele,Modele,case,LevRange[0],LevRange[1])
#	ename = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1950-2016.%s-%shPa.%sN.NDJFM.p' % (Modele,Modele,case,LevRange[0],LevRange[1])
#	ename = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2005.%s-%shPa.%sN.DJF.p' % (Modele,Modele,case,LevRange[0],LevRange[1])
#	nname = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1975-2001.%s-%shPa.%sN.DJF.p' % (Modeln,Modeln,case,LevRange[0],LevRange[1])
#	nname = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2005.%s-%shPa.%sN.DJF.p' % (Modeln,Modeln,case,LevRange[0],LevRange[1])
	vqe,de = unpick(ename)
	vqe    = vqe/scale_factor
	i0     = len(de)/2
	vqe    = ff.interpolate(vqe)
	he     = ff.pdf(vqe[:i0,:],edges,integrate=integrate)/len(de[:i0])
#	vqn,dn = unpick(nname)
#	vqn    = ff.interpolate(vqn)
#	hn     = ff.pdf(vqn,edges,integrate=integrate)/len(dn)
	hn     = ff.pdf(vqe[i0:,:],edges,integrate=integrate)/len(de[i0:])
        ff.plot(hn,edges,vmin=vminf,vmax=vmaxf,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source=Modele,sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,Modele),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh)
	ff.plot(hn-he,edges,vmin=vmind,vmax=vmaxd,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source=Modele,sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.diff.pdf' % (case,Modele),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh)	
	"""

	"""
	# Trends
        ff               = FluxFile(Source='ERAInt',YearRange=YearRange,LevRange=LevRange,blat=blat,makefile=makefile,Season='DJF')
        vq1,dates1       = ff.flux(case=case)
        vq1              = ff.interpolate(vq1)
        vq1,dates1,years = ff.fluxToYears(vq1,dates1)
        F                = np.array([ff.pdf(np.array(vq1[i]),edges,integrate=integrate)[0]/len(dates1[i]) for i in range(len(vq1))])
        Fm               = F.mean(axis=0)
        n0,n1,n2         = F.shape
        Ft,P             = np.zeros((n1,n2)),np.zeros((n1,n2))
        for i in range(n1):
                for j in range(n2):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(years,F[:,i,j])
                        Ft[i,j]=slope
                        P[i,j]=p_value/2.
        F = Ft*10
        #F = np.ma.masked_where(P>0.05,F)
        ff.plot(Fm,edges,vmin=vminf,vmax=vmaxf,xlabel='Longitude',ylabel=ylabel,show=show,\
               Source='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,'ERAInt'),clabel=clabel,cmap=cmap,pthresh=pthresh)
        ff.plot(F,edges,vmin=-20,vmax=20,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.trend.pdf' % (case,'ERAInt'),clabel='Dailymean transport trend [Tg day$^{-1}$ decade$^{-1}$]',\
                cmap=pl.cm.RdBu_r,Hstd=None,pthresh=pthresh)
	"""

	ff = FluxFile(Source='CMCC-CESM',YearRange=(2060,2095),LevRange=LevRange,blat=blat,makefile=makefile)
	# ERAInt flux
	fename     = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2016.%s-%shPa.%sN.DJF.p' % (blat,'ERAInt','ERAInt',case,LevRange[0],LevRange[1],blat)
	vqe,datese = unpick(fename)
	vqe        = ff.interpolate(con_fac*vqe)
	he,hemean  = ff.pdf_new(vqe,edges,integrate=integrate,Field0=[])
	he         = he/len(datese)

	"""
	# ERAInt stationary V
	dsV = MMDataServer(Field='V')
	v   = np.array([dsV.getSeason(Year=year,Season='DJF') for year in range(1981,2016+1,1)])
	if (type(v)==np.ma.core.MaskedArray) and (v.mask.shape == v.data.shape):
	        v  = ff.fill(v.data,invalid=v.mask)
	v  = (v - v.mean(axis=-1)[:,:,:,np.newaxis]).mean(axis=0)
	v  = ff.interpolateND(v,dsV.lev,500,axis=0)
	v  = ff.interpolateND(v,dsV.lat,blat,axis=0)
	ve = ff.interpolateND(v,dsV.lon,range(5,355+10,10),axis=0)
	"""

	# Model mean PDF
	#Models   = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*') if x[28:] not in Modelsno]
	Models   = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CMCC-CESM', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GFDL-CM3', 'GFDL-ESM2G',\
		    'GFDL-ESM2M', 'inmcm4', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'CMCC-CM',\
		    'CMCC-CMS', 'FGOALS-g2', 'MPI-ESM-MR', 'IPSL-CM5A-LR']
	Hh,Hr,Vh = [],[],[]
	for ii in range(len(Models)):

		print '\n'
		print '%s of %s (%s) ...' % (ii+1,len(Models),Models[ii])

		"""
		dsV = cmipDataServer(Field='va',LevType='plev',Source=Models[ii],ExpType='rcp85',DataFreq='day',LevRange=(0,1000))
		v   = np.array([dsV.getSeason(Year=year,Season='DJF') for year in range(1981,2016+1,1)])
		print v.shape
		if (type(v)==np.ma.core.MaskedArray) and (v.mask.shape == v.data.shape):
			v  = ff.fill(v.data,invalid=v.mask)
		v = (v - v.mean(axis=-1)[:,:,:,np.newaxis]).mean(axis=0)
		v = ff.interpolateND(v,dsV.lev,500,axis=0)
		v = ff.interpolateND(v,dsV.lat,blat,axis=0)
		v = ff.interpolateND(v,dsV.lon,range(5,355+10,10),axis=0)
		Vh.append(v)
		"""

		Source     = Models[ii]
		fhname     = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.1981-2016.%s-%shPa.%sN.DJF.p' % (blat,Source,Source,case,LevRange[0],LevRange[1],blat)
		frname     = '/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.%s.2060-2095.%s-%shPa.%sN.DJF.p' % (blat,Source,Source,case,LevRange[0],LevRange[1],blat)
		print fhname
		print frname
		vqh,datesh = unpick(fhname)
		vqh        = ff.interpolate(con_fac*vqh)
		hh,hhmean  = ff.pdf_new(vqh,edges,integrate=integrate,Field0=[])
		hh         = hh/len(datesh)
                vqr,datesr = unpick(frname)
                vqr        = ff.interpolate(con_fac*vqr)
                hr,hrmean  = ff.pdf_new(vqr,edges,integrate=integrate,Field0=[])
		hr         = hr/len(datesr)
		Hh.append(hh)
		Hr.append(hr)
	Hh,Hr,Vh = np.array(Hh),np.array(Hr),np.array(Vh)
	dH       = np.ma.masked_where(Hr==0,Hr) - np.ma.masked_where(Hh==0,Hh)
	# Stipling
	n0,n1,n2      = dH.shape
	bins,lons     = np.linspace(edges[0]+cint/2.,edges[-1]-cint/2.,n1),np.linspace(-180+5,180-5,n2)
	ax_x,ax_y     = np.meshgrid(lons,bins)
	ax_x,ax_y     = ax_x[::-1,:],ax_y[::-1,:]
	stipx,stipy   = stipling(dH,xx=ax_x,yy=ax_y,thresh=0.8)
	stipx1,stipy1 = stipling(Hh-he[np.newaxis,:,:],xx=ax_x,yy=ax_y,thresh=0.8)
	# Means and standard deviation
	Hm = Hh.mean(axis=0)
	Hs = dH.std(axis=0)
	#Vm = Vh.mean(axis=0)
	dH = dH.mean(axis=0)

	# Stationary wave bias
	#dV = list(Vm - ve)
	#dV = dV[18:] + dV[0:18]	# -175 -> 175
	#dV = dV[10:] + dV[0:10]	# -75  -> -85
	#dV = np.array(dV)

	# Plot mean and std of biases
	ff.plot(Hm,edges,vmin=vminf/2,vmax=vmaxf/2,xlabel='Longitude',ylabel=ylabel,show=show,\
		Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,'all'),clabel=clabel,cmap=cmap,pthresh=pthresh,bin_integral=True)
        ff.plot(dH/Hm,edges,vmin=-1,vmax=1,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.frac.pdf' % (case,'all'),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh)
	ff.plot(dH,edges,vmin=vmind/2,vmax=vmaxd/2,xlabel='Longitude',ylabel=ylabel,show=show,\
		Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.diff.pdf' % (case,'all'),clabel=clabel,cmap=pl.cm.RdBu_r,Hstd=None,stip=(stipx,stipy),pthresh=pthresh,bin_integral=True)
	ff.plot(Hs,edges,vmin=0,vmax=vmaxd/2,xlabel='Longitude',ylabel=ylabel,show=show,\
		Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.diff.std.pdf' % (case,'all'),clabel=clabel,cmap=pl.cm.OrRd,pthresh=pthresh)
        ff.plot(he,edges,vmin=vmind,vmax=vmaxd,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='all',sname='/mnt/climstorage/cian/scripts/figs/full/%s/%s.full.pdf' % (case,'ERAInt'),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh,bin_integral=False)
        ff.plot(Hm-he,edges,vmin=vmind/4,vmax=vmaxd/4,xlabel='Longitude',ylabel=ylabel,show=show,\
                Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.hdiff.pdf' % (case,'all'),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh)#,stip=(stipx1,stipy1))
#       ff.plot(Hm-he,edges,vmin=vmind/4,vmax=vmaxd/4,xlabel='Longitude',ylabel='%s, %s' % (ylabel,r'$\Delta{\overline{v}_{500hPa}}$ [cm s$^{-1}$]'),show=show,\
#               Source='all',sname='/mnt/climstorage/cian/scripts/figs/diff/%s/%s.hdiff.pdf' % (case,'all'),clabel=clabel,cmap=pl.cm.RdBu_r,pthresh=pthresh,bin_integral=True,stationary=dV)
