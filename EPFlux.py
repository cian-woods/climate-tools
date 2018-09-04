import sys
sys.path.insert(0, '/mnt/climstorage/cian/scripts')
from ReanalysisDataServer import *
from UnPickle import *

import numpy as np
import matplotlib.pyplot as pl

class EPFlux:

	def __init__(	self	):
		# DataServers
		self.Fdiv_ds = DataServer(Field='Fdiv',LatRange=(60,89),LevRange=(10,500))
		self.Fphi_ds = DataServer(Field='Fphi',LatRange=(60,89),LevRange=(10,500))
		self.Fp_ds   = DataServer(Field='Fp',LatRange=(60,89),LevRange=(10,500))
		# Attributes/conversion factors
	        self.sid    = 24*3600		# seconds in day
	        self.lat    = self.Fdiv_ds.lat
	        self.lev    = self.Fdiv_ds.lev
	        self.loglev = np.exp(np.linspace(np.log(10),np.log(500),20))
	        self.xs     = [np.argmin((self.lev-self.loglev[ii])**2) for ii in range(len(self.loglev))]
	        self.N      = np.ones(len(self.lev))
	        self.N[np.where(self.lev<100)] = 3
	        self.N      = self.N[:,np.newaxis]
		# Intrusions
		self.G,self.Q,self.D       = unpick('../intrusions/ERAInt_intrusions.1980-2015.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p')
		self.maxdates,self.maxlons = self.datesAtMaxInten(self.G,self.Q,self.D)

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

	def getClims(self,years,Season='DJF'):
		# Climatologies
		Fdiv_clim = np.array([self.Fdiv_ds.getDataSnaps(Year=year,Season=Season).mean(axis=0) for year in years]).mean(axis=0)
		Fphi_clim = np.array([self.Fphi_ds.getDataSnaps(Year=year,Season=Season).mean(axis=0) for year in years]).mean(axis=0)
		Fp_clim   = np.array([  self.Fp_ds.getDataSnaps(Year=year,Season=Season).mean(axis=0) for year in years]).mean(axis=0)		
		return self.sid*Fdiv_clim, self.N*Fphi_clim, -1*self.N*Fp_clim

	def laggedComp(self):
		inds  = [i for i in range(len(self.maxlons)) if 150 <=self.maxlons[i] <= 220]
		dates = [self.maxdates[ind] for ind in inds]
		# Data
		case  = 'PAC'
        	days  = np.arange(-5,15+1,1)
        	#dates = unpick('../WarmArctic/dates.%s.50.ERAInt.p' % (case))
        	Fdiv,Fphi,Fp = [],[],[]
        	for date0 in dates:
        	        hour0    = self.Fdiv_ds.getHours(*date0)
        	        datelist = [self.Fdiv_ds.getDate(hour0 + ii*24) for ii in days]
        	        try:
        	                Fdiv.append([self.Fdiv_ds.snapshot(*date) for date in datelist])
        	                Fphi.append([self.Fphi_ds.snapshot(*date) for date in datelist])
        	                Fp.append([self.Fp_ds.snapshot(*date) for date in datelist])
        	        except:
				pass
		print '%s of %s dates used in lagged composite' % (len(Fdiv),len(dates))
        	Fdiv,Fphi,Fp = self.sid*np.array(Fdiv).mean(axis=0),self.N*np.array(Fphi).mean(axis=0),-1*self.N*np.array(Fp).mean(axis=0)
		# Climatology
		years = range(1991,2016+1,1)
		Fdiv_clim,Fphi_clim,Fp_clim = self.getClims(years,Season='NDJFM')
		self.plot(Fdiv_clim,Fphi_clim,Fp_clim,cseq=np.arange(-20,20+4,4),cmap=pl.cm.RdBu_r,extend='both',\
			  scale=400,ref_arr=200,title='E-P flux NDJFM clim',savename='figs/EP/clim.pdf')
		# Anomalies
		Fdiv,Fphi,Fp = Fdiv-Fdiv_clim[np.newaxis,:,:],Fphi-Fphi_clim[np.newaxis,:,:],Fp-Fp_clim[np.newaxis,:,:]
		# Plot
		cseq = np.arange(-5,5+1,1)
		for ii in range(len(days)):
			self.plot(Fdiv[ii],Fphi[ii],Fp[ii],cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',\
				  scale=100,ref_arr=20,title='day %s' % (days[ii]),savename='figs/EP/%s/%s.pdf' % (case,ii+1))

	def plot(self,Fdiv,Fphi,Fp,cseq=14,cmap=pl.cm.RdBu_r,extend='both',scale=100,ref_arr=20,title='E-P flux',savename=None):
	        cf   = pl.contourf(self.lat,self.lev,Fdiv,cseq,cmap=pl.cm.RdBu_r,extend='both')
	        cbar = pl.colorbar()
	        cbar.set_label('Zonal wind acceleration [m s$^{-1}$ day$^{-1}$]')
	        Q  = pl.quiver(self.lat[::5],self.lev[self.xs],Fphi[self.xs,::5],Fp[self.xs,::5],units='inches',scale=100,\
	                        scale_units='inches',headwidth=2.5,headlength=2.5,headaxislength=2.5,pivot='tail',alpha=0.6)
	        qk = pl.quiverkey(Q, 0.2, 1.02, ref_arr, '%s' % (ref_arr), labelpos='W',fontproperties={'weight': 'bold'})
	        pl.ylabel('Pressure [hPa]')
	        pl.xlabel(r'Latitude [$^{\circ}$ N]')
	        pl.ylim(1000,10)
	        pl.yscale('log')
	        pl.xlim(-90,90)
	        pl.title(title)
		if savename != None: pl.savefig(savename,format='pdf')
	        pl.show()

if __name__ == "__main__":

	EP = EPFlux()
	EP.laggedComp()
