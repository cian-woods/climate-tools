from scipy import ndimage as nd
from scipy import interpolate
from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from toPick import *
from UnPickle import *
from stipling import stipling

import sys,glob
import numpy as np
import matplotlib.pyplot as pl

class Reynolds:

	def __init__(	self,
			Season = 'DJF',
			case   = 'vq'):
		# Attributes
		self.Season = Season
		self.case   = case
		if self.case == 'vq': self.Fieldf, self.orr, self.units, self.sf, self.Ef = ['q','hus'],    'moisture',            'kg s$^{-1}$ m$^{-1}$',     1,   1e03
		if self.case == 'vT': self.Fieldf, self.orr, self.units, self.sf, self.Ef = ['T', 'ta'], 'temperature', '10$^{3}$ K kg s$^{-1}$ m$^{-1}$', 1e-03, 334e03
		self.proj   = LambertProjector(boundinglat=25,resolution=200.)
		# Target grid
		self.lats = np.arange(0,90+1,1)
		self.lons = np.arange(0,360,1)
		self.levs = np.arange(0,1000+50,50)
		self.dlev = self.levs[-1]-self.levs[0]

	def extendField(self,field):
	        fieldend = field[:,:,:,-1][:,:,:,np.newaxis]
	        field    = np.append(field,fieldend,axis=3)
	        field    = np.append(fieldend,field,axis=3)
	        return field

        def moving_sum(self,x,y,n=3):
                if n%2 == 0: return 0,0
                ret = np.cumsum(y, dtype=float)
                ret[n:] = ret[n:] - ret[:-n]
                return x[n/2:-n/2 + 1],ret[n - 1:]

        def moving_average(self,x,y,n=3):
                if n%2 == 0: return 0,0
                ret = np.cumsum(y, dtype=float)
                ret[n:] = ret[n:] - ret[:-n]
                return x[n/2:-n/2 + 1],ret[n - 1:]/n

	def rotate_series(self,lons,data,lon0):
		x    = np.argmin((lons-lon0)**2)
		data = np.append(data[x:],data[:x],axis=0)
		lons = np.append(lons[x:]-360,lons[:x],axis=0)
		return lons,data

	def get_dM(self,levs):
		# Mass of girdbox on each level
		dP_ = []
		dP  = np.diff(levs)/2.
		dP_.append(dP[0])
		for i in range(len(dP)-1):
		        dP_.append(dP[i]+dP[i+1])
		dP_.append(dP[-1])
		dP_ = np.array(dP_)
		dM = 100.*dP_/9.80665
		return dM

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

	def za(self,field):
        	return field - field.mean(axis=-1)[...,np.newaxis]	

        def zm(self,field):
                return field.mean(axis=-1)

        def seasonDecompPeixotoOort(self,YearRange,Source):
                fname = '/mnt/climstorage/cian/PeixotoOort/%s.%s.%s-%s.hPa.%s-%s.%s.p' % (Source,self.case,self.levs[0],self.levs[-1],YearRange[0],YearRange[1],self.Season)
                print fname
                if not os.path.isfile(fname):
                        years = range(YearRange[0],YearRange[1]+1,1)
			dM    = self.get_dM(self.levs)
                        if Source == 'ERAInt':
                                vds = reDataServer(Field='V',LevType='plev',Source='ERAInt',LevRange=(0,1000),LatRange=(0,90))
                                fds = reDataServer(Field=self.Fieldf[0],LevType='plev',Source='ERAInt',LevRange=(0,1000),LatRange=(0,90))
                        else:
                                vds = cmipDataServer(Field='va' ,LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
                                fds = cmipDataServer(Field=self.Fieldf[1],LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
                	# Extend lon axis
                	lonres = np.abs(vds.lon[1]-vds.lon[0])
                	dlon   = np.append(vds.lon,vds.lon[-1]+lonres)
                	dlon   = np.append(vds.lon[0]-lonres,dlon)
                        VM,FM  = np.zeros((0,len(self.levs),len(self.lats),len(self.lons))),np.zeros((0,len(self.levs),len(self.lats),len(self.lons)))
                        VV,FF  = np.zeros((0,len(self.levs),len(self.lats),len(self.lons))),np.zeros((0,len(self.levs),len(self.lats),len(self.lons)))
                        for year in years:
                                print year
                                v = vds.getDataSnaps(Year=year,Season=self.Season,dailymean=True)
                                f = fds.getDataSnaps(Year=year,Season=self.Season,dailymean=True)
                                if (type(f)==np.ma.core.MaskedArray) and (f.mask.shape == f.data.shape):
                                        print 'Filling arrays ...'
                                        v  = self.fill(v.data,invalid=v.mask)
                                        f  = self.fill(f.data,invalid=f.mask)
				v = self.extendField(v)
                                v = self.interpolateND(v,dlon,self.lons,axis=3)
                                v = self.interpolateND(v,vds.lat,self.lats,axis=2)
                                v = self.interpolateND(v,vds.lev,self.levs,axis=1)
				f = self.extendField(f)
				f = self.interpolateND(f,dlon,self.lons,axis=3)
                                f = self.interpolateND(f,fds.lat,self.lats,axis=2)
                                f = self.interpolateND(f,fds.lev,self.levs,axis=1)
                                # Season mean and anomalies
                                vm,fm = v.mean(axis=0)[np.newaxis,:,:,:],f.mean(axis=0)[np.newaxis,:,:,:]
                                VM = np.append(VM, vm, axis=0)
                                FM = np.append(FM, fm, axis=0)
                                VV = np.append(VV,  v, axis=0)
                                FF = np.append(FF,  f, axis=0)
                        vm,fm = VM.mean(axis=0),FM.mean(axis=0)
                        va,fa = VV - vm[np.newaxis,:,:,:],FF - fm[np.newaxis,:,:,:]
                        # Products
			vfm = ((      (VV*FF).mean(axis=0)).mean(axis=-1) ).mean(axis=0)*(self.dlev*100/9.81) # Total
			t1  = ((      (va*fa).mean(axis=0)).mean(axis=-1) ).mean(axis=0)*(self.dlev*100/9.81) # Transients
			t2  = (( self.za(vm)*self.za(fm)).mean(axis=-1) ).mean(axis=0)*(self.dlev*100/9.81)   # Stationary
			t3  = (( self.zm(vm) - self.zm(vm).mean(axis=0)[np.newaxis,:])*(self.zm(fm) - self.zm(fm).mean(axis=0)[np.newaxis,:])).mean(axis=0)*(self.dlev*100/9.81) # Overturning
			t4  = (( self.zm(vm).mean(axis=0))*(self.zm(fm).mean(axis=0)) )*(self.dlev*100/9.81)	# residual
                        # Save file
                        lats = self.lats
                        toPick([lats,vfm,t1,t2,t3,t4],fname)
                else:
                        # Open file
                        lats,vfm,t1,t2,t3,t4 = unpick(fname)

		#pl.figure(1)
		#pl.plot(vfm-t4,lats,'k'  ,linewidth=2.0,alpha=0.75,label='$\{[\overline{vq}]\} - \{[\overline{v}]\}\{[\overline{q}]\}$')
		#pl.plot(t1 ,lats,'b--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v^{\prime}q^{\prime}}]\}$')
		#pl.plot(t2 ,lats,'r--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}^{*}\overline{q}^{*}]\}$')
		#pl.plot(t3 ,lats,'g--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}]^{\circ}[\overline{q}]^{\circ}\}$')
		#pl.plot(t4-t4,lats,'y--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}]\}\{[\overline{q}]\} = 0$')
		#pl.xlabel('Zonal-mean meridional moisture flux [kg s$^{-1}$ m$^{-1}$]')
		#pl.ylabel('Latitude')
		#pl.grid()
		#pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		#pl.show()

		return lats,vfm,t1,t2,t3,t4

	def peixotoOortBiases(self,YearRange):
		fname = '/mnt/climstorage/cian/PeixotoOort/%s.vq.%s-%s.hPa.%s-%s.%s.p' % ('ERAInt',self.levs[0],self.levs[-1],YearRange[0],YearRange[1],self.Season)
		print fname

		d_vq    = Dataset('/mnt/climstorage/cian/vq/vq.mon.mean.new.nc','r')
		vq_wne  = d_vq.variables['vq'][:]
		months  = np.tile(np.arange(1,12+1,1),vq_wne.shape[0]/12.)
		xs      = np.where((months==1)|(months==2)|(months==12)==True)[0]
		vq_wne  = vq_wne[xs,::2,:]
		lats_wn = d_vq.variables['lat'][::2]
		wvns    = d_vq.variables['wavenumber'][:]
		d_vq.close()
		vq_wne   = vq_wne[2:-1,:,:].reshape((-1,3,len(lats_wn),len(wvns))).mean(axis=1)
		years_wn = range(1980,2016+1,1)
		xs0,xs1  = years_wn.index(YearRange[0]),years_wn.index(YearRange[1])
		vq_wne   = vq_wne[xs0:xs1+1,:,1].mean(axis=0)#.sum(axis=-1)

		latse,vqme,t1e,t2e,t3e,t4e = unpick(fname)
		VQM,T1,T2,T3,T4,VQ_WN = [],[],[],[],[],[]
		Models = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*')]	
		for Model in Models:
			try:
				fname = '/mnt/climstorage/cian/PeixotoOort/%s.vq.%s-%s.hPa.%s-%s.%s.p' % (Model,self.levs[0],self.levs[-1],YearRange[0],YearRange[1],self.Season)	
				lats,vqm,t1,t2,t3,t4 = unpick(fname)

				"""
	                	pl.figure(1)
	                	ix = np.argmin((lats-0)**2)
	                	pl.plot(vqm[ix:],lats[ix:],'k',linewidth=2.0,alpha=0.75,label='$\{[\overline{vq}]\}$')
	                	pl.plot( t1[ix:],lats[ix:],'b--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v^{\prime}q^{\prime}}]\}$')
	                	pl.plot( t2[ix:],lats[ix:],'r--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}^{*}\overline{q}^{*}]\}$')
	                	pl.plot( t3[ix:],lats[ix:],'g--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}]^{\circ}[\overline{q}]^{\circ}\}$')
	                	pl.plot( t4[ix:],lats[ix:],'y--',linewidth=0.7,alpha=0.45,label='$\{[\overline{v}]\}\{[\overline{q}]\}$')
	                	pl.xlabel('Zonal-mean meridional moisture flux [kg s$^{-1}$ m$^{-1}$]')
	                	pl.ylabel('Latitude')
	                	pl.grid()
	                	pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
				pl.xlim(-80,40)
				pl.title(Model)
				pl.show()
				"""

		        	d_vq    = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/vq/vq.mon.mean.new.nc' % (Model),'r')
		        	vq_wn   = d_vq.variables['vq'][:]
				months = np.tile(np.arange(1,12+1,1),vq_wn.shape[0]/12.)
				xs     = np.where((months==1)|(months==2)|(months==12)==True)[0]
				vq_wn   = vq_wn[xs,::2,:]
				lats_wn = d_vq.variables['lat'][::2]
				wvns    = d_vq.variables['wavenumber'][:]
				d_vq.close()
				vq_wn    = vq_wn[2:-1,:,:].reshape((-1,3,len(lats_wn),len(wvns))).mean(axis=1)
				years_wn = range(1981,2100+1,1)
				xs0,xs1  = years_wn.index(YearRange[0]),years_wn.index(YearRange[1])
				vq_wn    = vq_wn[xs0:xs1+1,:,1].mean(axis=0)#.sum(axis=-1)

				VQM.append(vqm)
				T1.append(t1)
				T2.append(t2)
				T3.append(t3)
				T4.append(t4)
				VQ_WN.append(vq_wn)
				print Model
			except:
				print 'Model %s did not complete...' % (Model)
		VQM,T1,T2,T3,T4,VQ_WN = np.array(VQM).mean(axis=0),np.array(T1).mean(axis=0),np.array(T2).mean(axis=0),np.array(T3).mean(axis=0),np.array(T4).mean(axis=0),np.array(VQ_WN).mean(axis=0)

		lats_m,vqme_m = self.moving_average(lats,vqme,9)
		lats_m,VQM_m  = self.moving_average(lats, VQM,9)
		lats_m,t1e_m  = self.moving_average(lats,t1e,9)
		lats_m,t2e_m  = self.moving_average(lats,t2e,9)
		lats_m,t3e_m  = self.moving_average(lats,t3e,9)
		lats_m,t4e_m  = self.moving_average(lats,t4e,9)

		lats_m,dVQM_m  = self.moving_average(lats,VQM-vqme,9)
                lats_m,dT1_m  = self.moving_average(lats,T1-t1e,9)
                lats_m,dT2_m  = self.moving_average(lats,T2-t2e,9)
                lats_m,dT3_m  = self.moving_average(lats,T3-t3e,9)
                lats_m,dT4_m  = self.moving_average(lats,T4-t4e,9)

                pl.figure(1)
                pl.plot(vqme_m[:],lats_m[:],'k.-',linewidth=2.0,alpha=0.75,label='$\{[\overline{vq}]\}$')
                pl.plot( t1e_m[:],lats_m[:],'b--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v^{\prime}q^{\prime}}]\}$')
                pl.plot( t2e_m[:],lats_m[:],'r--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}^{*}\overline{q}^{*}]\}$')
                pl.plot( t3e_m[:],lats_m[:],'g--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}]^{\circ}[\overline{q}]^{\circ}\}$')
                pl.plot( t4e_m[:],lats_m[:],'y--',linewidth=1.0,alpha=0.55,label='$\{[\overline{v}]\}\{[\overline{q}]\}$')
                pl.xlabel('Zonal-mean meridional moisture flux [kg s$^{-1}$ m$^{-1}$]')
                pl.ylabel('Latitude')
                pl.grid()
                pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.title(vqme[:].sum())

		pl.figure(2)
                ix = np.argmin((lats-0)**2)
                pl.plot(dVQM_m[:],lats_m[:],'k'  ,linewidth=2.00,alpha=0.85,label='$\Delta{\{[\overline{vq}]\}}$')
                pl.plot( dT1_m[:],lats_m[:],'b'  ,linewidth=2.00,alpha=0.65,label='$\Delta{\{[\overline{v^{\prime}q^{\prime}}]\}}$')
                pl.plot( dT2_m[:],lats_m[:],'r'  ,linewidth=2.00,alpha=0.65,label='$\Delta{\{[\overline{v}^{*}\overline{q}^{*}]\}}$')
                pl.plot( dT3_m[:],lats_m[:],'g--',linewidth=1.10,alpha=0.80,label='$\Delta{\{[\overline{v}]^{\circ}[\overline{q}]^{\circ}\}}$')
                pl.plot( dT4_m[:],lats_m[:],'y--',linewidth=0.80,alpha=0.80,label='$\Delta{\{[\overline{v}]\}\{[\overline{q}]\}}$')
                pl.xlabel('Zonal-mean meridional moisture flux bias [kg s$^{-1}$ m$^{-1}$]')
                pl.ylabel('Latitude')
                pl.grid()
                pl.ylim(4,85)
                #pl.xlim(-1,2.5)
                pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/peixoto_oort.pdf',format='pdf')
                pl.show()

	def seasonDecompReynolds(self,YearRange,Source):
		fname = '/mnt/climstorage/cian/Reynolds/%s.%s.%s-%s.hPa.%s-%s.%s.p' % (Source,self.case,self.levs[0],self.levs[-1],YearRange[0],YearRange[1],self.Season)
		print fname
		if not os.path.isfile(fname):
			years = range(YearRange[0],YearRange[1]+1,1)
			if Source == 'ERAInt':
				vds = reDataServer(Field='V',LevType='plev',Source='ERAInt',LevRange=(0,1000),LatRange=(0,90))
				fds = reDataServer(Field=self.Fieldf[0],LevType='plev',Source='ERAInt',LevRange=(0,1000),LatRange=(0,90))
			else:
				vds = cmipDataServer(Field='va' ,LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
				fds = cmipDataServer(Field=self.Fieldf[1],LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
                        # Extend lon axis
                        lonres = np.abs(vds.lon[1]-vds.lon[0])
                        dlon   = np.append(vds.lon,vds.lon[-1]+lonres)
                        dlon   = np.append(vds.lon[0]-lonres,dlon)
			VM,FM = np.zeros((0,len(self.levs),len(self.lats),len(self.lons))),np.zeros((0,len(self.levs),len(self.lats),len(self.lons)))
			VV,FF = np.zeros((0,len(self.levs),len(self.lats),len(self.lons))),np.zeros((0,len(self.levs),len(self.lats),len(self.lons)))
			for year in years:
				print year
				v = vds.getDataSnaps(Year=year,Season=self.Season,dailymean=True)
				f = fds.getDataSnaps(Year=year,Season=self.Season,dailymean=True)
				if (type(f)==np.ma.core.MaskedArray) and (f.mask.shape == f.data.shape):
					print 'Filling arrays ...'
				        v  = self.fill(v.data,invalid=v.mask)
				        f  = self.fill(f.data,invalid=f.mask)	
					#vxs,fxs       = np.where(v.mask==True),np.where(f.mask==True)
					#v,f           = f.data,f.data
					#v[vxs],f[fxs] = 0,0
				v = self.extendField(v)
				v = self.interpolateND(v,dlon,self.lons,axis=3)
				v = self.interpolateND(v,vds.lat,self.lats,axis=2)
				v = self.interpolateND(v,vds.lev,self.levs,axis=1)
				f = self.extendField(f)
				f = self.interpolateND(f,dlon,self.lons,axis=3)
				f = self.interpolateND(f,fds.lat,self.lats,axis=2)
				f = self.interpolateND(f,fds.lev,self.levs,axis=1)
				# Season mean and anomalies
				vm,fm = v.mean(axis=0)[np.newaxis,:,:,:],f.mean(axis=0)[np.newaxis,:,:,:]
	                        VM = np.append(VM, vm, axis=0)
	                        FM = np.append(FM, fm, axis=0)
				VV = np.append(VV,  v, axis=0)
				FF = np.append(FF,  f, axis=0)
			vm,fm = VM.mean(axis=0),FM.mean(axis=0)
			va,fa = VV - vm[np.newaxis,:,:,:],FF - fm[np.newaxis,:,:,:]
			# Products
			vmfm = vm*fm
			vafa = (va*fa).mean(axis=0)
			vfm  = (VV*FF).mean(axis=0)
			# Save file
			levs,lats,lons = self.levs,self.lats,self.lons	
			toPick([levs,lats,lons,vmfm,vafa,vfm,vm,fm],fname)
                        # Close DataServer Files
			if Source != 'ERAInt':
                        	vds.closeFiles()
                        	fds.closeFiles()
		else:
			# Open file
			levs,lats,lons,vmfm,vafa,vfm,vm,fm = unpick(fname)
		return levs,lats,lons,vmfm,vafa,vfm,vm,fm

	def reynoldsBiases(self,YearRange1,YearRange2):
		# Open file(s) with ERAInt data
		POfname = '/mnt/climstorage/cian/PeixotoOort/%s.%s.%s-%s.hPa.%s-%s.%s.p' % ('ERAInt',self.case,self.levs[0],self.levs[-1],YearRange[0],YearRange[1],self.Season)
		latsPOe,vfmPOe,POt1e,POt2e,POt3e,POt4e = unpick(POfname)
		levs,lats,lons,vqse,vqte,vqme,vme,qme = self.seasonDecompReynolds(YearRange,'ERAInt')
		# Vertical integral
		dM             = self.get_dM(levs)	
		vqse,vqte,vqme = self.sf*vqse.mean(axis=0)*self.dlev*100/9.81,self.sf*vqte.mean(axis=0)*self.dlev*100/9.81,self.sf*vqme.mean(axis=0)*self.dlev*100/9.81
		vqme,vqse      = vqme - self.sf*POt3e[:,np.newaxis],vqse - self.sf*POt3e[:,np.newaxis]
		vmE,qmE        = vme,qme

		models = []
		# Open CMIP5 files
		VQS1,VQT1,VQM1,VM1,QM1 = [],[],[],[],[]
		VQS2,VQT2,VQM2,VM2,QM2 = [],[],[],[],[]
		Models = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*')]	
		for Model in Models:
			try:
				# Open file(s)
				POfname = '/mnt/climstorage/cian/PeixotoOort/%s.%s.%s-%s.hPa.%s-%s.%s.p' % (Model,self.case,self.levs[0],self.levs[-1],YearRange1[0],YearRange1[1],self.Season)
				latsPO,vfmPO,POt1,POt2,POt3,POt4 = unpick(POfname)
				levs1,lats1,lons1,vqs1,vqt1,vqm1,vm1,qm1 = self.seasonDecompReynolds(YearRange1,Model)
				levs2,lats2,lons2,vqs2,vqt2,vqm2,vm2,qm2 = self.seasonDecompReynolds(YearRange2,Model)
				# Vertical integral
				dM = self.get_dM(levs)
				vqs1,vqt1,vqm1 = self.sf*vqs1.mean(axis=0)*self.dlev*100/9.81,self.sf*vqt1.mean(axis=0)*self.dlev*100/9.81,self.sf*vqm1.mean(axis=0)*self.dlev*100/9.81
				vqs2,vqt2,vqm2 = self.sf*vqs2.mean(axis=0)*self.dlev*100/9.81,self.sf*vqt2.mean(axis=0)*self.dlev*100/9.81,self.sf*vqm2.mean(axis=0)*self.dlev*100/9.81
				vqm1,vqs1      = vqm1 - self.sf*POt3[:,np.newaxis],vqs1 - self.sf*POt3[:,np.newaxis]
				# Append to lists
				VQS1.append(vqs1)
				VQT1.append(vqt1)
				VQM1.append(vqm1)
				VM1.append(vm1)
				QM1.append(qm1)
                                VQS2.append(vqs2)
                                VQT2.append(vqt2)
                                VQM2.append(vqm2)
                                VM2.append(vm2)
                                QM2.append(qm2)
				models.append(Model)
				print Model
			except:
				print 'Model %s did not complete...' % (Model)
		print models
		n = len(VQS1)
		VQS1,VQT1,VQM1,VM1,QM1      = np.array(VQS1),np.array(VQT1),np.array(VQM1),np.array(VM1),np.array(QM1)
		VQS2,VQT2,VQM2,VM2,QM2      = np.array(VQS2),np.array(VQT2),np.array(VQM2),np.array(VM2),np.array(QM2)
		VQS1m,VQT1m,VQM1m,VM1m,QM1m = VQS1.mean(axis=0),VQT1.mean(axis=0),VQM1.mean(axis=0),VM1.mean(axis=0),QM1.mean(axis=0)
		VQS2m,VQT2m,VQM2m,VM2m,QM2m = VQS2.mean(axis=0),VQT2.mean(axis=0),VQM2.mean(axis=0),VM2.mean(axis=0),QM2.mean(axis=0)

		latx1,latx2 = np.argmin((lats-60)**2),np.argmin((lats-70)**2)
		px1         = np.argmin((levs-500)**2)
		# Terms (spherical) 1
		dvqm1      = VQM1m-vqme
		dvqs1      = VQS1m-vqse
		dvqt1      = VQT1m-vqte
		dvm1       = VM1m-vmE
		dqm1       = QM1m-qmE
                vbar_qstr1 = dM.sum()*(vmE*(QM1m - qmE)).mean(axis=0)
                qbar_vstr1 = dM.sum()*(qmE*(VM1m - vmE)).mean(axis=0)
                vstr_qstr1 = dM.sum()*((VM1m - vmE)*(QM1m - qmE)).mean(axis=0)
		epsilon1   = dvqm1 - vbar_qstr1 - qbar_vstr1 - dvqt1
                # Terms (spherical) 2
                dvqm2      = VQM2m-VQM1m
                dvqs2      = VQS2m-VQS1m
                dvqt2      = VQT2m-VQT1m
                vbar_qstr2 = dM.sum()*(VM1m*(QM2m - QM1m)).mean(axis=0)
                qbar_vstr2 = dM.sum()*(QM1m*(VM2m - VM1m)).mean(axis=0)
                vstr_qstr2 = dM.sum()*((VM2m - VM1m)*(QM2m - QM1m)).mean(axis=0)
                epsilon2   = dvqm2 - vbar_qstr2 - qbar_vstr2 - dvqt2
		# Meridonal section mean (Rotate 80W -> 280E) 1
		lons_m,dvqm_m1      =      self.rotate_series(lons,dvqm1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,dvqs_m1      =      self.rotate_series(lons,dvqs1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,dvqt_m1      =      self.rotate_series(lons,dvqt1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,vbar_qstr_m1 = self.rotate_series(lons,vbar_qstr1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,qbar_vstr_m1 = self.rotate_series(lons,qbar_vstr1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,vstr_qstr_m1 = self.rotate_series(lons,vstr_qstr1[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,epsilon_m1   =   self.rotate_series(lons,epsilon1[latx1:latx2+1,:].mean(axis=0),280)
		lons_m,dv_m1        =   self.rotate_series(lons,dvm1[px1,latx1:latx2+1,:].mean(axis=0),280)
                # Meridonal section mean (Rotate 80W -> 280E) 2
                lons_m,dvqm_m2      =      self.rotate_series(lons,dvqm2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,dvqs_m2      =      self.rotate_series(lons,dvqs2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,dvqt_m2      =      self.rotate_series(lons,dvqt2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,vbar_qstr_m2 = self.rotate_series(lons,vbar_qstr2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,qbar_vstr_m2 = self.rotate_series(lons,qbar_vstr2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,vstr_qstr_m2 = self.rotate_series(lons,vstr_qstr2[latx1:latx2+1,:].mean(axis=0),280)
                lons_m,epsilon_m2   =   self.rotate_series(lons,epsilon2[latx1:latx2+1,:].mean(axis=0),280)
		# Moving average 1
                lons_ms,dvqm_m1      = self.moving_average(lons_m,      dvqm_m1,9)
                lons_ms,dvqs_m1      = self.moving_average(lons_m,      dvqs_m1,9)
                lons_ms,dvqt_m1      = self.moving_average(lons_m,      dvqt_m1,9)
                lons_ms,vbar_qstr_m1 = self.moving_average(lons_m, vbar_qstr_m1,9)
                lons_ms,qbar_vstr_m1 = self.moving_average(lons_m, qbar_vstr_m1,9)
                lons_ms,vstr_qstr_m1 = self.moving_average(lons_m, vstr_qstr_m1,9)
                lons_ms,epsilon_m1   = self.moving_average(lons_m,   epsilon_m1,9)
		lons_ms,dv_m1        = self.moving_average(lons_m,        dv_m1,9)
                # Moving average 2
                lons_ms,dvqm_m2      = self.moving_average(lons_m,      dvqm_m2,9)
                lons_ms,dvqs_m2      = self.moving_average(lons_m,      dvqs_m2,9)
                lons_ms,dvqt_m2      = self.moving_average(lons_m,      dvqt_m2,9)
                lons_ms,vbar_qstr_m2 = self.moving_average(lons_m, vbar_qstr_m2,9)
                lons_ms,qbar_vstr_m2 = self.moving_average(lons_m, qbar_vstr_m2,9)
                lons_ms,vstr_qstr_m2 = self.moving_average(lons_m, vstr_qstr_m2,9)
                lons_ms,epsilon_m2   = self.moving_average(lons_m,   epsilon_m2,9)
		# Terms (cartesian)
                dvqmC1      = self.proj(     dvqm1,lons,lats)
                dvqsC1      = self.proj(     dvqs1,lons,lats)
                dvqtC1      = self.proj(     dvqt1,lons,lats)
                vbar_qstrC1 = self.proj(vbar_qstr1,lons,lats)
                qbar_vstrC1 = self.proj(qbar_vstr1,lons,lats)
                vstr_qstrC1 = self.proj(vstr_qstr1,lons,lats)
                epsilonC1   = self.proj(  epsilon1,lons,lats)
		# Project to Lambert
                vqse,vqte,vqmE,vmE,qmE = self.proj(vqse,lons,lats),self.proj(vqte,lons,lats),self.proj(vqme,lons,lats),self.proj(vmE,lons,lats),self.proj(qmE,lons,lats)
                VQS1,VQT1,VQM1,VM1,QM1 = self.proj(VQS1,lons,lats),self.proj(VQT1,lons,lats),self.proj(VQM1,lons,lats),self.proj(VM1,lons,lats),self.proj(QM1,lons,lats)
		VQS2,VQT2,VQM2,VM2,QM2 = self.proj(VQS2,lons,lats),self.proj(VQT2,lons,lats),self.proj(VQM2,lons,lats),self.proj(VM2,lons,lats),self.proj(QM2,lons,lats)
                # Stipling
                sxM,syM   = stipling(VQM1-vqmE,xx=self.proj.x,yy=self.proj.y,thresh=0.8)
                sxS,syS   = stipling(VQS1-vqse,xx=self.proj.x,yy=self.proj.y,thresh=0.8)
                sxT,syT   = stipling(VQT1-vqte,xx=self.proj.x,yy=self.proj.y,thresh=0.8)
                sxvM,syvM = stipling((VM1-vmE[np.newaxis,:,:,:]).mean(axis=1),xx=self.proj.x,yy=self.proj.y,thresh=0.8)
                sxqM,syqM = stipling((QM1-qmE[np.newaxis,:,:,:]).mean(axis=1),xx=self.proj.x,yy=self.proj.y,thresh=0.8)

		# Plot     
                cseq = np.arange(-28,28+4,4)
                pl.figure(1)
                cf   = pl.contourf(self.proj.x,self.proj.y,dvqmC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$\Delta{\overline{v.q}}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.plot(sxM,syM,'k.',alpha=0.5)
                pl.title('Total bias: %s' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.tot.bias.pdf',format='pdf')
		pl.close()

		pl.figure(2)
                cf   = pl.contourf(self.proj.x,self.proj.y,dvqsC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$\Delta{\overline{v}.\overline{q}}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.plot(sxS,syS,'k.',alpha=0.5)
                pl.title('Stationary bias: %s' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.stat.bias.pdf',format='pdf')
		pl.close()

                pl.figure(3)
		cseq0 = np.arange(-6,6+1,1)
                cf    = pl.contourf(self.proj.x,self.proj.y,dvqtC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label(r'$\Delta{\overline{v^{\prime}.q^{\prime}}}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('Transient bias: %s' % (n))
		pl.plot(sxT,syT,'k.',alpha=0.5)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.trans.bias.pdf',format='pdf')
		pl.close()

                pl.figure(4)
                cf   = pl.contourf(self.proj.x,self.proj.y,vbar_qstrC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$\overline{v_{R}}.q^{*}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.plot(sxqM,syqM,'k.',alpha=0.5)
                pl.title('Stationary bias: %s' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.qstat.bias.pdf',format='pdf')
		pl.close()

                pl.figure(5)
                cf   = pl.contourf(self.proj.x,self.proj.y,qbar_vstrC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$\overline{q_{R}}.v^{*}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
		pl.plot(sxvM,syvM,'k.',alpha=0.5)
                pl.title('Stationary bias: %s' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.vstat.bias.pdf',format='pdf')
		pl.close()

                pl.figure(6)
                cf   = pl.contourf(self.proj.x,self.proj.y,vstr_qstrC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$v^{*}.q^{*}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('Stationary bias: %s' % (n))
                pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.vqstr.bias.pdf',format='pdf')
		pl.close()

                pl.figure(7)
		cf   = pl.contourf(self.proj.x,self.proj.y,epsilonC1,cseq,cmap=pl.cm.coolwarm,extend='both')
                cbar = pl.colorbar(cf)
                cbar.set_label(r'$\epsilon$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('Stationary bias: %s' % (n))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/vq.epsilon.bias.pdf',format='pdf')
		pl.close()
		
		pl.figure(8)
		ix = 0#np.argmin((lats-0)**2)
		pl.plot(                 dvqm1[ix:,:].mean(axis=-1), lats[ix:], 'k'   , linewidth=2.00, alpha=0.85, label=r'$\Delta{\overline{v.q}}$')
		pl.plot(                 dvqt1[ix:,:].mean(axis=-1), lats[ix:], 'b'   , linewidth=2.00, alpha=0.65, label=r'$\Delta{\overline{v^{\prime}.q^{\prime}}}$')
		pl.plot(            vbar_qstr1[ix:,:].mean(axis=-1), lats[ix:], 'b--' , linewidth=1.40, alpha=0.80, label=r'$\overline{v_{R}}.q^{*}$')
		pl.plot(            qbar_vstr1[ix:,:].mean(axis=-1), lats[ix:], 'r--' , linewidth=1.40, alpha=0.80, label=r'$\overline{q_{R}}.v^{*}$')
		pl.plot(              epsilon1[ix:,:].mean(axis=-1), lats[ix:], 'k--', linewidth=0.45, alpha=0.80, label=r'$\epsilon$')
		pl.xlabel('Zonal-mean moisture flux bias [%s]' % (self.units))
		pl.ylabel('Latitude')
		pl.grid()
		pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/zm.pdf',format='pdf')
		pl.close()

                pl.figure(9)
		cseqf = np.arange(-100,100+10,10)
                cf    = pl.contourf(self.proj.x,self.proj.y,vqmE,cseqf,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label(r'$\overline{v.q}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2') 
                pl.title('Total climatology: ERAInt')
		pl.close()

                pl.figure(10)
                cseqf = np.arange(-50,50+10,10)
                cf    = pl.contourf(self.proj.x,self.proj.y,vqse,cseqf,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label(r'$\overline{v}.\overline{q}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('Stationary climatology: ERAInt')
		pl.close()

                pl.figure(11)
                cseqf = np.arange(-50,50+10,10)
                cf    = pl.contourf(self.proj.x,self.proj.y,vqte,cseqf,cmap=pl.cm.coolwarm,extend='both')
                cbar  = pl.colorbar(cf)
                cbar.set_label(r'$\overline{v^{\prime}.q^{\prime}}$ [%s]' % (self.units))
                self.proj.m.drawparallels([70,80],latmax=90)
                drawCoastlinesNoRivers(self.proj.m,linewidth=0.6,color='0.2')
                pl.title('Transient climatology: ERAInt')
		pl.close()

		fig,ax1 = pl.subplots(num=12)
                ax1.plot(lons_ms,      dvqm_m1, 'k'  , linewidth=2.00, alpha=0.85, label=r'$\Delta{\overline{v.q}}$')
                ax1.plot(lons_ms,      dvqt_m1, 'b'  , linewidth=2.00, alpha=0.65, label=r'$\Delta{\overline{v^{\prime}.q^{\prime}}}$')
                ax1.plot(lons_ms, vbar_qstr_m1, 'b--', linewidth=1.40, alpha=0.80, label=r'$\overline{v_{R}}.q^{*}$')
                ax1.plot(lons_ms, qbar_vstr_m1, 'r--', linewidth=1.40, alpha=0.80, label=r'$\overline{q_{R}}.v^{*}$')
		ax1.plot(lons_ms,   epsilon_m1, 'k--', linewidth=0.45, alpha=0.8, label=r'$\epsilon$')
		ax1.set_xlabel('Longitude')
		ax1.set_ylabel('Meridional moisture flux bias [%s]' % (self.units))
		ax1.grid()
		ax1.set_xlim(lons_ms[0],lons_ms[-1])
		ax1.set_ylim(-20,20)
		ax1.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
		ax1.set_xticks([-60,0,60,120,180,240])
		ax1.set_xticklabels(['-60','0','60','120','180','-120'])
		pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/mm1.pdf',format='pdf')

                fig,ax1 = pl.subplots(num=13)
                ax1.plot(lons_ms,dvqm_m2, 'k', linewidth=2.00, alpha=0.85, label=r'$\Delta{\overline{v.q}}$ (%s)' % (round(dvqm_m2.mean(),2)))
                ax1.plot(lons_ms,dvqt_m2, 'b', linewidth=2.00, alpha=0.65, label=r'$\Delta{\overline{v^{\prime}.q^{\prime}}}$ (%s)' % (round(dvqt_m2.mean(),2)))
                ax1.plot(lons_ms, vbar_qstr_m2, 'b--', linewidth=1.40, alpha=0.80, label=r'$\overline{v_{R}}.q^{*}$')
                ax1.plot(lons_ms, qbar_vstr_m2, 'r--', linewidth=1.40, alpha=0.80, label=r'$\overline{q_{R}}.v^{*}$')
		ax1.plot(lons_ms, vbar_qstr_m2+qbar_vstr_m2, 'y', linewidth=1.10, alpha=0.80, label=r'$\Delta{\overline{v}.\overline{q}}$')
                ax1.plot(lons_ms,   epsilon_m2, 'k--', linewidth=0.45, alpha=0.8, label=r'$\epsilon$')
                ax1.set_xlabel('Longitude')
                ax1.set_ylabel('Meridional moisture flux change [%s]' % (self.units))
                ax1.grid()
                ax1.set_xlim(lons_ms[0],lons_ms[-1])
                ax1.set_ylim(-20,20)
                ax1.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
                ax1.set_xticks([-60,0,60,120,180,240])
		ax1.set_xticklabels(['-60','0','60','120','180','-120'])
                pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/Reynolds/mm2.pdf',format='pdf')

		pl.show()


if __name__ == "__main__":

	from netCDF4 import Dataset

	YearRange = (int(sys.argv[1]),int(sys.argv[2]))
	Season    = str(sys.argv[3])
	Source    = str(sys.argv[4])
	case      = str(sys.argv[5])
	Re        = Reynolds(Season=Season,case=case)

	#Re.reynoldsBiases(YearRange,(2060,2095))
	Re.peixotoOortBiases(YearRange)
	#Re.seasonDecompReynolds(YearRange=YearRange,Source=Source)
	#Re.seasonDecompPeixotoOort(YearRange=YearRange,Source=Source)


	# Make yearly files
	#Models   = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*')]
	#for Model in Models[:15]:
	#	for year in range(YearRange[0],YearRange[1]+1,1):
	#		try:
	#			Re.seasonDecompReynolds(YearRange=(year,year),Source=Model)
	#		except:
	#			pass

	# Make files
	#Models   = [x[28:] for x in glob.glob('/mnt/climstorage/cian/rcp85/*')]
	#print Models
	#for Model in Models:
	#	try:
	#		levs,lats,lons,vmqm,vaqa,vqm,vm,qm = Re.seasonDecompReynolds(YearRange=YearRange,Source=Model)
	#		#lats,vqm,t1,t2,t3,t4 = Re.seasonDecompPeixotoOort(YearRange=YearRange,Source=Model)
	#	except:
	#		pass

