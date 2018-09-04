from netCDF4 import *
from scipy import stats
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer

import glob
import numpy as np
import matplotlib.pyplot as pl

SMs = {}
SMs['D']         = [12]
SMs['J']         = [1]
SMs['F']         = [2]
SMs['M']         = [3]
SMs['JF']        = [1,2]
SMs['DJ']        = [12,1]
SMs['DJF']       = [12,1,2]
SMs['JFD']       = [1,2,12]
SMs['JFM']       = [1,2,3]
SMs['JJA']       = [6,7,8]
SMs['SON']       = [9,10,11]
SMs['MAM']       = [3,4,5]
SMs['Annual']    = range(1,12+1,1)
SMs['Annual_wc'] = range(7,12+1,1) + range(1,6+1,1)
        
class ArgoServer:

	def __init__(	self,
			Field    =     'TEMP',
			LatRange =   (-90,90),
			LonRange = (-180,180),
			LevRange =   (0,2000)	):

		# Attributes
		self.Field    = Field
		self.LatRange = LatRange
		self.LonRange = LonRange
		self.LevRange = LevRange
		self.Files    = glob.glob('/mnt/climstorage/cian/Argo/ISAS15/*%s.nc' % (self.Field))
		self.dates    = [f.split('/')[-1].split('_')[2] for f in self.Files]
		self.dates    = [(int(date[0:4]),int(date[4:6])) for date in self.dates]
		print 'Data from %s to %s' % (self.dates[0],self.dates[-1])
		self.bmap     = Basemap(projection='cyl',llcrnrlat=LatRange[0],urcrnrlat=LatRange[-1],llcrnrlon=LonRange[0],urcrnrlon=LonRange[-1],resolution='i')
		# DataServer
		self.stds = MMDataServer(Field='sst',LatRange=LatRange,LonRange=(0,360))
		# Axes
		d_tmp    = Dataset(self.Files[0],'r')
		self.lat = d_tmp.variables['latitude'][:]
		self.lon = d_tmp.variables['longitude'][:]
		self.z   = d_tmp.variables['depth'][:]
		# Slice into ranges
		self.i0,self.i1,self.j0,self.j1,self.k0,self.k1,self.lat,self.lon,self.z = self.sliceAxes()
		d_tmp.close()
		self.nz,self.nlat,self.nlon = len(self.z),len(self.lat),len(self.lon)
		# Basemap coordinates
		self.j,self.k = self.bmap(self.lon,self.lat)
		# Mass of girdbox on each level
		dD_ = []
		dD  = np.diff(self.z)/2.
		dD_.append(self.z[0]+dD[0])
		for i in range(len(dD)-1):
		        dD_.append(dD[i]+dD[i+1])
		dD_.append(dD[-1])
		dD_ = np.array(dD_)
		self.dM = 1028.*dD_

	def sliceAxes(self):
		i0,i1 = np.argmin((self.lat-self.LatRange[0])**2),np.argmin((self.lat-self.LatRange[1])**2)
		j0,j1 = np.argmin((self.lon-self.LonRange[0])**2),np.argmin((self.lon-self.LonRange[1])**2)
		k0,k1 = np.argmin((self.z-self.LevRange[0])**2),np.argmin((self.z-self.LevRange[1])**2)
		return i0,i1,j0,j1,k0,k1,self.lat[i0:i1+1],self.lon[j0:j1+1],self.z[k0:k1+1]

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

	def getMonth(self,Year,Month):
		ix = self.dates.index((Year,Month))
		d  = Dataset(self.Files[ix],'r')
		f  = d.variables[self.Field][:,self.k0:self.k1+1,self.i0:self.i1+1,self.j0:self.j1+1].squeeze()
		d.close()
		return f

	def getSeason(self,Year,Season,mean=True):
		months     = SMs[Season]
		yearmonths = []
		for i in months:
			if i > months[-1]:
				yearmonths.append((Year-1,i))
			else:
				yearmonths.append((Year,i))
		print 'Getting data timeslice from %s --> %s ...' % (yearmonths[0],yearmonths[-1])
		f = np.ma.masked_array([self.getMonth(year,month) for year,month in yearmonths])
		if mean == True:
			return f.mean(axis=0)
		else:
			return f

	def OHC(self,YearRange,Season):
		years  = range(YearRange[0],YearRange[1]+1,1)
		years_ = range(YearRange[0],2012+1,1)
		# Get Argo data
		T     = np.ma.masked_array([self.getSeason(Year=year,Season=Season,mean=False) for year in years]) + 273.15
		T_sm  = T.mean(axis=0)
		T_ss  = T.std(axis=0)
		Ta    = T - T_sm[np.newaxis,:,:,:,:]
		Ta    = Ta.reshape((-1,self.nz,self.nlat,self.nlon))
		# Get sst data
		sst         = np.ma.masked_array([self.stds.getSeason(Year=year,Season=Season,mean=False) for year in years_])
		sst,lonsout = basemap.shiftgrid(180,sst,self.stds.lon,start=True,cyclic=360)
		lonsout     = lonsout - 360
		jds,kds     = self.bmap(lonsout,self.stds.lat)
		sst_sm      = sst.mean(axis=0)
		ssta        = sst - sst_sm[np.newaxis,:,:,:]
		ssta        = ssta.reshape((-1,self.stds.nlat,self.stds.nlon)) 

		# Slice
		latx0,lonx0     = np.argmin((self.lat-70)**2),np.argmin((self.lon-5)**2)
		latx1,lonx1     = np.argmin((self.lat-76)**2),np.argmin((self.lon-20)**2)
		T_tza           = Ta[:,:,latx0:latx1+1,lonx0:lonx1+1].mean(axis=-1).mean(axis=-1)
		latx0ds,lonx0ds = np.argmin((self.stds.lat-70)**2),np.argmin((lonsout-5)**2)
                latx1ds,lonx1ds = np.argmin((self.stds.lat-76)**2),np.argmin((lonsout-20)**2)
                sst_ta          = ssta[:,latx0ds:latx1ds+1,lonx0ds:lonx1ds+1].mean(axis=-1).mean(axis=-1)

		print self.lat[latx0],self.lat[latx1],self.lon[lonx0],self.lon[lonx1]
		print self.stds.lat[latx0ds],self.stds.lat[latx1ds],lonsout[lonx0ds],lonsout[lonx1ds]

		pl.figure(1)
                T_tz0a  = T_tza[:,0]
		T_smz0   = np.tile(T_sm[:,0,latx0:latx1+1,lonx0:lonx1+1].mean(axis=-1).mean(axis=-1),( len(years))  )
		T_ssz0   = np.tile(T_ss[:,0,latx0:latx1+1,lonx0:lonx1+1].mean(axis=-1).mean(axis=-1),( len(years))  )
		pl.plot(range(len(ssta)),sst_ta,'r',linewidth=1.75,alpha=0.65)
		pl.plot(range(len(Ta)),T_tz0a,'b',linewidth=1.75,alpha=0.65)
		pl.plot(range(len(Ta)),0.65*T_ssz0,'b--',linewidth=1.00,alpha=0.45)
		pl.plot(range(len(Ta)),-0.65*T_ssz0,'b--',linewidth=1.00,alpha=0.45)
		#pl.plot([0,len(Ta)],[T_tz0am+0.65*T_tz0as,T_tz0am+0.65*T_tz0as],'b--',linewidth=1.00,alpha=0.45)
		#pl.plot([0,len(Ta)],[T_tz0am-0.65*T_tz0as,T_tz0am-0.65*T_tz0as],'b--',linewidth=1.00,alpha=0.45)
		pl.grid()
		pl.ylabel('Sea surface temperature [K]')
		pl.xlabel('Year')
		#pl.title('Corr. coeff. = %s' % (np.corrcoef(sst_ta,T_tza[:,0])[0][1]))

                fig       = pl.figure(2)
                ax        = fig.add_subplot(111)
                cseq      = 13#np.arange(272,280+0.5,0.5)
                cf        = ax.contourf(jds,kds,sst.mean(axis=0).mean(axis=0),cseq,cmap=pl.cm.OrRd,extend='both')
                cbar      = pl.colorbar(cf,orientation='horizontal')
                pl.plot([self.j[lonx0],self.j[lonx1]],[self.k[latx0],self.k[latx0]],'k',alpha=0.75)
		pl.plot([self.j[lonx0],self.j[lonx1]],[self.k[latx1],self.k[latx1]],'k',alpha=0.75)
		pl.plot([self.j[lonx0],self.j[lonx0]],[self.k[latx0],self.k[latx1]],'k',alpha=0.75)
		pl.plot([self.j[lonx1],self.j[lonx1]],[self.k[latx0],self.k[latx1]],'k',alpha=0.75)
                cbar.set_label('Sea surface temperature [K]')
                self.bmap.drawcoastlines(linewidth=0.5,color='0.5')
                self.bmap.fillcontinents(color='0.85')
                #self.bmap.drawmeridians([0,20,40,60,80,100])
                #self.bmap.drawparallels([70,75,80])
                ax.set_aspect('auto')
                pl.title('ERAInt sea surface temperature')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/Argo/Ts.clim.ERAInt.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season), format='pdf')

                fig       = pl.figure(3)
                ax        = fig.add_subplot(111)
                cseq      = 13#np.arange(272,280+0.5,0.5)
                cf        = ax.contourf(self.j,self.k,T[:,:,0,:,:].mean(axis=0).mean(axis=0),cseq,cmap=pl.cm.OrRd,extend='both')
                cbar      = pl.colorbar(cf,orientation='horizontal')
                pl.plot([self.j[lonx0],self.j[lonx1]],[self.k[latx0],self.k[latx0]],'k',alpha=0.75)
                pl.plot([self.j[lonx0],self.j[lonx1]],[self.k[latx1],self.k[latx1]],'k',alpha=0.75)
                pl.plot([self.j[lonx0],self.j[lonx0]],[self.k[latx0],self.k[latx1]],'k',alpha=0.75)
                pl.plot([self.j[lonx1],self.j[lonx1]],[self.k[latx0],self.k[latx1]],'k',alpha=0.75)
                cbar.set_label('Surface level temperature [K]')
                self.bmap.drawcoastlines(linewidth=0.5,color='0.5')
                self.bmap.fillcontinents(color='0.85')
                #self.bmap.drawmeridians([0,20,40,60,80,100])
                #self.bmap.drawparallels([70,75,80])
                ax.set_aspect('auto')
		pl.title('Argo surface temperature')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/Argo/Ts.clim.argo.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season), format='pdf')

		pl.figure(4)
		T_tza = np.rollaxis(T_tza,1,0)
		cseq  = np.arange(-0.35,0.35+0.05,0.05)
		cf    = pl.contourf(range(len(Ta)),self.z,T_tza,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cbar  = pl.colorbar(orientation='horizontal')
		cbar.set_label('Temperature anomaly [K]')
		pl.ylabel('Depth [m]')
		pl.ylim(self.z[-1],self.z[0])
		pl.xticks(range(len(Ta))[::len(SMs[Season])],[r'%s/%02d' % (i,j) for i in years for j in SMs[Season]][::len(SMs[Season])],rotation=65,fontsize=9)
		#monthlist = np.tile(SMs[Season],(len(years)))
		#for i in range(len(Ta)):
		#	if monthlist[i] == 2:
		#		pl.plot([i,i],[self.z[0],self.z[-1]],'k--',linewidth=1.0,alpha=0.45)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/Argo/T.anom.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season), format='pdf')

		pl.figure(5)
		xs      = np.where((np.abs(T_tz0a)>=0.65*T_ssz0)==True)[0]
		T_pm    = (T_tza/T_tza[0,:])[:,xs].mean(axis=1)
		T_ps    = (T_tza/T_tza[0,:])[:,xs].std(axis=1)
		pl.plot([T_pm.mean(),T_pm.mean()],[self.z[0],self.z[-1]],'k--',linewidth=1,alpha=0.5)
		pl.plot(T_pm,self.z,'k',linewidth=2,alpha=0.7)
		pl.ylabel('Depth [m]')
		pl.xlabel('${T^{\prime}}_{Z}/{T^{\prime}}_{surf}$ [ ]')
		pl.ylim(self.z[-1],self.z[0])
		pl.xlim(0,1.5)
		pl.title('%s profiles' % (len(xs)))
		pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/Argo/T.anom.norm.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season), format='pdf')

		pl.show()

if __name__ == "__main__":

	AS = ArgoServer(  Field    =        'TEMP',
			  LatRange =       (68,83),
			  LonRange =       (0,100),
			  LevRange =       (0,250)	)

	f = AS.OHC((2002,2015),'JFD')

	
