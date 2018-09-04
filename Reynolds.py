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
		return lats,vfm,t1,t2,t3,t4

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

