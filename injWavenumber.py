from netCDF4 import Dataset
from ReanalysisDataServer import DataServer as reDataServer
from UnPickle import *
from scipy import stats

import numpy as np
import matplotlib.pyplot as pl

# Initial file for axes, plus other attributes
masktype    = 'injection'
years       = range(1980,2005+1,1)
dayinds     = range(60) + range(335,364+1,1)
dv          = Dataset('../v_decomp/v_decomp_1980.nc','r')
time        = dv.variables['time'][:]
lat         = dv.variables['lat'][:]
lev         = dv.variables['lev'][:]
wvn         = dv.variables['wavenumber'][:15]
dv.close()
# ERAInt DataServer
dsERA = reDataServer(Field='T2',LevType='surface_analysis')
# Mass of girdbox on each level
dP_ = []
dP  = np.diff(lev)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665
# ERAInt data
print 'ERAInt'
F        = np.zeros((0,len(lat),15))
datelist = []
for year in years:
        # Data files
        dq = Dataset('../q_decomp/q_decomp_%s.nc' % (year),'r')
        dv = Dataset('../v_decomp/v_decomp_%s.nc' % (year),'r')
        # Extract variables
        qAmp = dq.variables['Amp'][dayinds,:,:,:15]
        qPhi = dq.variables['Phi'][dayinds,:,:,:15]
        vAmp = dv.variables['Amp'][dayinds,:,:,:15]
        vPhi = dv.variables['Phi'][dayinds,:,:,:15]
        # Add dates to datelist
        datelist = datelist + [dsERA.getDate(24*day)[0:3] for day in dq.variables['time'][dayinds]]
        # Close files
        dq.close()
        dv.close()
        # Flux
	F = np.append(F,((vAmp*qAmp*np.cos(vPhi-qPhi)/2)*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1),axis=0)
if masktype == 'injection':
	#G,Q,D = unpick('/qnap/cian/cmip/intrusions/ERAInt_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p')
	G,Q,D = unpick('../intrusions/ERAInt_intrusions.1980-2015.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p')
	E     = np.array([sum([sum(q1) for q1 in q]) for q in Q])
	D     = [D[i] for i in np.where(E>=stats.scoreatpercentile(E,0))[0]]
	D     = [D[i][j][0:3] for i in range(len(D)) for j in range(len(D[i]))]
	xs_   = [datelist.index(injdate) for injdate in D if datelist.count(injdate)>0]
	xs    = []
	for i in range(len(xs_)):
	        if xs.count(xs_[i]) == 0:
	                xs.append(xs_[i])
elif masktype == 'sigma':
	F70 = F[:,70,:].sum(axis=-1)
	xs  = np.where(F70>F70.mean()+0.0*F70.std())[0]
xs_ = [i for i in range(len(F)) if i not in xs]
print len(F),len(xs),len(xs)/len(F)

sf = 2260e03/1e06
allF = sf*F[xs_,:,:].sum(axis=0)
injF = sf*F[xs,:,:].sum(axis=0)
dF   = injF - allF

pl.figure(1)
cf   = pl.contourf(wvn,lat,allF/len(F),np.arange(0,4+0.5,0.5),cmap=pl.cm.OrRd,extend='max')
cbar = pl.colorbar(cf)
cbar.set_label('Zonal-mean meridional mositure flux [10$^{6}$ W m$^{-1}$]')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Latitude')
pl.title('All days')

pl.figure(2)
cf   = pl.contourf(wvn,lat,injF/len(F),np.arange(0,4+0.5,0.5),cmap=pl.cm.OrRd,extend='max')
cbar = pl.colorbar(cf)
cbar.set_label('Zonal-mean meridional mositure flux [10$^{6}$ W m$^{-1}$]')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Latitude')
pl.title('Injection days')

pl.figure(3)
cf   = pl.contourf(wvn,lat,dF/len(F),np.arange(-3,3+0.5,0.5),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('Zonal-mean meridional mositure flux [10$^{6}$ W m$^{-1}$]')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Latitude')
pl.title('Difference')

pl.show()
