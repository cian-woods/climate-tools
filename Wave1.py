import sys
sys.path.insert(0, '/mnt/climstorage/cian/scripts')
from netCDF4 import Dataset as Dataset
from LambertProjector import *
from cmipDataServer import DataServer as cmipDataServer

import matplotlib.pyplot as pl
import numpy as np

class Wave1:

	def __init__(	self	):
		# Projector and Dataset
		self.proj = LambertProjector(boundinglat=50,resolution=100.)
		f         = Dataset('/mnt/climstorage/gabriele/Data/Warm_Arctic/slp_deseason_filt_y_NM.nc','r')
		# Data variables
		self.lats = f.variables['lat_slp'][::-1]
		self.lons = f.variables['lon_slp'][:]
		self.s    = f.variables['slp_a'][:,::-1,:]
		self.s    = np.rollaxis(self.s,2,0)
		self.s    = np.rollaxis(self.s,2,1)
		# fft attributes
		self.ftn = self.s.shape[-1]
		self.i0  = int(self.ftn/2) + 1
		# blat
		self.x   = np.argmin((self.lats-80)**2)
		self.y   = np.argmin((self.lats-90)**2)
		# Dates
		self.cmds = cmipDataServer(Field='psl',LevType='surface',Source='BNU-ESM',DataFreq='day')
		datelist  = []
		for year in range(1950,2017+1,1):
			datelist = datelist + self.cmds.getDateList(year,Season='NDJFM')[1:]
		self.datelist = datelist[60:-91+1]

	def wave1(self,v):
	        # Takes in (time)x(lat)x(lon) field and returns the wave1 pattern
	        # Divide by n becasue discrete transform and extract wave1 (i0 index after shift)
	        Fk  = np.fft.fft(v,axis=2)/self.ftn
	        Fk  = np.fft.fftshift(Fk,axes=(2,))[:,:,self.i0]
	        # Calcutae amplitide and phase of wave1
	        amp = 2*np.abs(Fk)
	        ang = np.arctan2(Fk.imag,Fk.real)
		# Make wave1 pattern
		field = np.array([-1*amp*np.cos(i + ang) for i in np.linspace(-np.pi,np.pi,self.ftn)])
		field = np.rollaxis(field,1,0)
		field = np.rollaxis(field,2,1)
	        return amp,ang,field

	def rankPhase(self,s):
		# Wave 1 properties and rank
		amp,ang,field = self.wave1(s)
		rank_index   = [i for (r,i) in sorted(zip((ang[:,self.x:self.y+1].mean(axis=1)+1.5)**2,range(len(amp))))]
#		rank_ang     = [r for (r,i) in sorted(zip((ang[:,self.x:self.y+1].mean(axis=1)+1.5)**2,range(len(amp))))]
#		rank_amp     = [amp[i,self.x:self.y+1].mean() for i in rank_index]
#		rank_index   = [i for (r,i) in sorted(zip(rank_amp,rank_index))][::-1]
		#rank_index  = [i for (r,i) in sorted(zip(ang[:,self.x:self.y+1].mean(axis=1),range(len(amp))))][::-1]
		# Plot
		cseq = np.arange(-3000,3000+300,300)
		for ii in range(100):
		#n    = np.argmin((ang[:,self.x:self.y+1].mean(axis=1)**2))
			n    = rank_index[ii]
			s0   = self.proj(field[n,:,:],self.lons,self.lats)
			s1   = self.proj(s[n,:,:],self.lons,self.lats)
			cf   = pl.contourf(self.proj.x,self.proj.y,s1,cseq,cmap=pl.cm.RdBu_r,extend='both')
			cl   = pl.contour(self.proj.x,self.proj.y,s0,levels=cseq,colors='k',linewidths=0.8,alpha=0.8)
			pl.clabel(cl,fmt='%2.0f',colors='k',fontsize=10)			
			cbar = pl.colorbar(cf)
			self.proj.m.drawcoastlines()
			self.proj.m.drawparallels([70,80],latmax=90)
			pl.title('%s radians, %s, %s' % (ang[n,self.x:self.y+1].mean(),self.datelist[n],n))
			pl.show()

	def makeFile(self,field,times,lats,lons,FileName,case):
		amp,ang,field = self.wave1(field)
		print amp.shape,ang.shape,field.shape

	print 'Creating %s ...' % (FileName)

	if case == 'vq': orr = 'meridional'
	if case == 'uq': orr = 'zonal'
	ntime,nlon,nlat = len(times),len(lons),len(lats)

	# Create file
	File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
	# Define some global attribs
File.Conventions='COARDS'
# Time is record dimension
File.createDimension('time',ntime)
var           = File.createVariable('time','d',('time',))
var.long_name = 'Time'
var.units     = 'hours since 1800-01-01 00:00:0.0'
var[:]        = times
# Horizontal axes
File.createDimension('lat',nlat)
var           = File.createVariable('lat','f',('lat',))
var.long_name = 'Latitude'
var.units     = 'degrees_north'
var[:]        = lats.astype('f')
File.createDimension('lon',nlon)
var           = File.createVariable('lon','f',('lon',))
var.long_name = 'Longitude'
var.units     = 'degrees_east'
var[:]        = lons.astype('f')
# Create Variables
var           = File.createVariable(case,'f',('time','lat','lon'))
var.long_name = '4xDaily vertically integrated %s moisture flux' % (orr)
var.units     = 'kg s**-1 m**-1'
var[:]        = field
# Close file
File.close()

if __name__ == "__main__":

	wv = Wave1()
	#wv.rankPhase(wv.s)
	wv.makeFile(wv.s)

