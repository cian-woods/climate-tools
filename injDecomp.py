from netCDF4 import Dataset
from ReanalysisDataServer import DataServer as reDataServer
from UnPickle import *
from scipy import stats

import numpy as np
import matplotlib.pyplot as pl

def datesAtMaxInten(G,Q,Dates):
        dates,lons = [],[]
        for i in range(len(Q)):
                ymax    = [max(Q[i][j]) for j in range(len(Q[i]))]
                yargmax = [np.argmax(Q[i][j]) for j in range(len(Q[i]))]
                k       = np.argmax(ymax)
                l       = yargmax[k]
                dates.append(Dates[i][k])
                lons.append(G[i][k][l])
        return dates,lons

# ReanalysisDataServer
dsERA = reDataServer(Field='V',LevType='plev')
# Injections
G,Q,D      = unpick('../intrusions/ERAInt_intrusions.1980-2015.NDJFM.6x6hr.9deg.200.6dt.20.5.filtered.p')
#G,Q,D      = unpick('/qnap/cian/cmip/intrusions/ERAInt_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p')
E          = np.array([sum([sum(q1) for q1 in q]) for q in Q])
injxs      = np.where(E>=stats.scoreatpercentile(E,80))[0]
G,Q,D      = [G[ix] for ix in injxs],[Q[ix] for ix in injxs],[D[ix] for ix in injxs]
dates,lons = datesAtMaxInten(G,Q,D)
dates      = [dates[i] for i in range(len(dates)) if (lons[i] < 90) or (330 < lons[i] < 360)]
#dates     = unpick('../WarmArctic/dates.warm.1000.ERAInt.p')[:100]
dates      = [(date[0],date[1],date[2],0) for date in dates]
days       = [dsERA.getDays(*date) for date in dates]
print days
# Initial file for axes, plus other attributes
years   = range(1979,2016+1,1)
dv      = Dataset('../v_decomp/v_decomp_1980.nc','r')
lat     = dv.variables['lat'][:]
lev     = dv.variables['lev'][:]
wvn     = dv.variables['wavenumber'][:15]
levx    = 9
latx    = 70
wvnx    = range(0,14+1,1)
dayinds = range(60) + range(335,364+1,1)

# ERAInt data
print 'ERAInt'
v_Amp0 = np.zeros((0,15))
v_Amp1 = np.zeros((0,15))
v_Phi0 = np.zeros((0,15))
v_Phi1 = np.zeros((0,15))
q_Amp0 = np.zeros((0,15))
q_Amp1 = np.zeros((0,15))
q_Phi0 = np.zeros((0,15))
q_Phi1 = np.zeros((0,15))
for year in years:
        # Data files
        dv = Dataset('../v_decomp/v_decomp_%s.nc' % (year),'r')
	dq = Dataset('../q_decomp/q_decomp_%s.nc' % (year),'r')
        # Extract variables 
        vAmp = dv.variables['Amp'][dayinds,levx,latx,wvnx]
        vPhi = dv.variables['Phi'][dayinds,levx,latx,wvnx]
        qAmp = dq.variables['Amp'][dayinds,levx,latx,wvnx]
        qPhi = dq.variables['Phi'][dayinds,levx,latx,wvnx]
        # Open time axis
	daylist = list(dv.variables['time'][dayinds])
        # Close files
        dv.close()
	dq.close()
	# Extract timesteps corresponding to injections at max intensity
	xdays  = list(set(daylist)&set(days))
	xs     = [daylist.index(xday) for xday in xdays]
	v_Amp0 = np.append(v_Amp0,vAmp[:,:],axis=0)
	v_Amp1 = np.append(v_Amp1,vAmp[xs,:],axis=0)
        v_Phi0 = np.append(v_Phi0,vPhi[:,:],axis=0)
        v_Phi1 = np.append(v_Phi1,vPhi[xs,:],axis=0)
        q_Amp0 = np.append(q_Amp0,qAmp[:,:],axis=0)
        q_Amp1 = np.append(q_Amp1,qAmp[xs,:],axis=0)
        q_Phi0 = np.append(q_Phi0,qPhi[:,:],axis=0)
        q_Phi1 = np.append(q_Phi1,qPhi[xs,:],axis=0)
v_Amp0m,v_Amp0s = v_Amp0.mean(axis=0),v_Amp0.std(axis=0)
v_Amp1m,v_Amp1s = v_Amp1.mean(axis=0),v_Amp1.std(axis=0)
v_Phi0m,v_Phi0s = v_Phi0.mean(axis=0),v_Phi0.std(axis=0)
v_Phi1m,v_Phi1s = v_Phi1.mean(axis=0),v_Phi1.std(axis=0)
q_Amp0m,q_Amp0s = q_Amp0.mean(axis=0),q_Amp0.std(axis=0)
q_Amp1m,q_Amp1s = q_Amp1.mean(axis=0),q_Amp1.std(axis=0)
q_Phi0m,q_Phi0s = q_Phi0.mean(axis=0),q_Phi0.std(axis=0)
q_Phi1m,q_Phi1s = q_Phi1.mean(axis=0),q_Phi1.std(axis=0)

pl.figure(1)
#pl.plot(range(1,15+1,1),v_Amp0m,'r--')
#pl.fill_between(range(1,15+1,1),v_Amp0m-v_Amp0s,v_Amp0m+v_Amp0s,color='r',alpha=0.1)
#pl.errorbar(range(1,15+1,1),v_Amp1m,v_Amp1s,fmt='o',color='k')
pl.plot(range(1,15+1,1),v_Amp0m-v_Amp0m,'r--')
pl.fill_between(range(1,15+1,1),-v_Amp0s,v_Amp0s,color='r',alpha=0.1)
pl.errorbar(range(1,15+1,1),v_Amp1m-v_Amp0m,v_Amp1s,fmt='o',color='k')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Meridional wind amplitiude [m s$^{-1}$]')
pl.xlim(0.5,15.5)
#pl.ylim(0,13)
pl.title('Climatology; %s days, Injections; %s days' % (len(v_Amp0),len(v_Amp1)))
pl.grid()
pl.savefig('figs/inj_Vamp_anom_P80.pdf', format='pdf')

pl.figure(2)
pl.plot(range(1,15+1,1),v_Phi0m,'r--')
pl.fill_between(range(1,15+1,1),v_Phi0m-v_Phi0s,v_Phi0m+v_Phi0s,color='r',alpha=0.1)
pl.errorbar(range(1,15+1,1),v_Phi1m,v_Phi1s,fmt='o',color='k')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Meridional wind phase [rad]')
pl.xlim(0.5,15.5)
#pl.ylim(0,13)
pl.title('Climatology; %s days, Injections; %s days' % (len(v_Phi0),len(v_Phi1)))
pl.grid()
pl.savefig('figs/inj_Vphi.pdf', format='pdf')

pl.figure(3)
pl.plot(range(1,15+1,1),q_Amp0m,'r--')
pl.fill_between(range(1,15+1,1),q_Amp0m-q_Amp0s,q_Amp0m+q_Amp0s,color='r',alpha=0.1)
pl.errorbar(range(1,15+1,1),q_Amp1m,q_Amp1s,fmt='o',color='k')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Specific humidity amplitiude [kg/kg]')
pl.xlim(0.5,15.5)
#pl.ylim(0,13)
pl.title('Climatology; %s days, Injections; %s days' % (len(q_Amp0),len(q_Amp1)))
pl.grid()
pl.savefig('figs/inj_Qamp.pdf', format='pdf')

pl.figure(4)
pl.plot(range(1,15+1,1),q_Phi0m,'r--')
pl.fill_between(range(1,15+1,1),q_Phi0m-q_Phi0s,q_Phi0m+q_Phi0s,color='r',alpha=0.1)
pl.errorbar(range(1,15+1,1),q_Phi1m,q_Phi1s,fmt='o',color='k')
pl.xlabel('Zonal wavenumber')
pl.ylabel('Specific humidity phase [rad]')
pl.xlim(0.5,15.5)
#pl.ylim(0,13)
pl.title('Climatology; %s days, Injections; %s days' % (len(q_Phi0),len(q_Phi1)))
pl.grid()
pl.savefig('figs/inj_Qphi.pdf', format='pdf')

pl.show()
