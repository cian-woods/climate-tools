import sys
sys.path.insert(0, '/mnt/climstorage/cian/scripts')
from UnPickle import *
from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *

import numpy as np
import matplotlib.pyplot as pl

# Attributes
case  = str(sys.argv[1])
proj  = LambertProjector(boundinglat=60,resolution=80.)
dates = unpick('dates.%s.50.ERAInt.p' % (case))[:20]
dsu   = reDataServer(Field='U',LevType='plev')
dsz   = reDataServer(Field='Z',LevType='plev')
dsT   = reDataServer(Field='T',LevType='plev')
mmdsu = MMDataServer(Field='U',LevType='plev')
mmdsz = MMDataServer(Field='Z',LevType='plev')
mmdsT = MMDataServer(Field='T',LevType='plev')
latx  = np.argmin((dsz.lat-80)**2)
px    = np.argmin((dsu.lev-500)**2)
days  = np.arange(-10,0+0.25,0.25)
# Climatology
Uclim  = proj(np.array([mmdsu.getSeason(Year=year,Season='NDJFM')[[9,14],:,:] for year in range(1980,2012+1,1)]).mean(axis=0),mmdsu.lon,mmdsu.lat)
dUclim = -1*np.diff(Uclim,axis=0).squeeze()
Zclim  = proj(np.array([mmdsz.getSeason(Year=year,Season='NDJFM')[[9,14],:,:] for year in range(1980,2012+1,1)]).mean(axis=0),mmdsz.lon,mmdsz.lat)
dZclim = -1*np.diff(Zclim,axis=0).squeeze()/(1000*9.8)
Tclim  = proj(np.array([mmdsT.getSeason(Year=year,Season='NDJFM')[14,:,:] for year in range(1980,2012+1,1)]).mean(axis=0),mmdsT.lon,mmdsT.lat)
dUdZclim = dUclim/dZclim
# Lagged composite
U,U0,UP,Z,T = [],[],[],[],[]
for date0 in dates:
	print date0
	hour0 = dsu.getHours(*date0)
	dlist = [dsu.getDate(hour0 + k*24) for k in days]
	U.append([dsu.snapshot(*date)[[9,14],:,:] for date in dlist])
	U0.append([dsz.snapshot(*date)[:,latx:,:] for date in dlist])
	UP.append([dsu.snapshot(*date)[:,:,:] for date in dlist])
	Z.append([dsz.snapshot(*date)[[9,14],:,:] for date in dlist])
	T.append([dsT.snapshot(*date)[14,:,:] for date in dlist])
U,U0,UP,Z,T = np.array(U),np.array(U0),np.array(UP),np.array(Z)/(1000*9.8),np.array(T)
UP          = proj(UP[:,:,[9,14],:,:].mean(axis=0),dsu.lon,dsu.lat)
U0          = U0.mean(axis=0).mean(axis=-2)
T           = proj(T.mean(axis=0),dsT.lon,dsT.lat)
dU,dZ       = -1*np.diff(U,axis=2).squeeze(),-1*np.diff(Z,axis=2).squeeze()
dUdZ        = dU/dZ
dUdZ        = proj(dUdZ.mean(axis=0),dsu.lon,dsu.lat)
#dUdZ       = dUdZ - dUdZclim[np.newaxis,:,:]
# Plot
cseq  = np.arange(-2.5,2.5+0.25,0.25)
cseq0 = np.arange(244,274+2,2)
cseq1 = np.arange(-2800,2800+400,400)
for i in range(dUdZ.shape[0]-8):
	pl.figure(1)
	#mask    = np.invert(((np.sign(UP[i,0,:,:]) + 1).astype(bool)) ^ ((np.sign(UP[i,1,:,:]) + 1).astype(bool)))
	#xx,yy   = np.ma.masked_array(proj.x,mask=mask),np.ma.masked_array(proj.y,mask=mask)
	xx,yy   = np.ma.masked_where(UP[i:i+9,0,:,:].mean(axis=0)>0,proj.x),np.ma.masked_where(UP[i:i+9,0,:,:].mean(axis=0)>0,proj.y)
	cf      = pl.contourf(proj.x,proj.y,dUdZ[i:i+9].mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
	cbar    = pl.colorbar(cf)
	pl.plot(xx[::2,::2],yy[::2,::2],'k.',alpha=0.5)
	cbar.set_label(r'$\frac{U{_{500}}-U{_{925}}}{Z{_{500}}-Z{_{925}}}$ [m s$^{-1}$ km$^{-1}$]')
	proj.m.drawcoastlines()
	proj.m.drawparallels([70,80],latmax=90)
	pl.title('day %s to %s' % (0.25*(i-40),0.25*(i-40)+2))
	pl.savefig('figs/lagged/baroclin/%s/%s.pdf' % (case,i+1), format='pdf')
	pl.close()
	pl.figure(2)
        cf      = pl.contourf(proj.x,proj.y,T[i:i+9].mean(axis=0),cseq0,cmap=pl.cm.OrRd,extend='max')
	cl      = pl.contour(proj.x,proj.y,T[i:i+9].mean(axis=0),[273],colors='k',linewidths=1.5)
        cbar    = pl.colorbar(cf) 
        cbar.set_label('Temperature 925hPa [K]')
        proj.m.drawcoastlines()
        proj.m.drawparallels([70,80,85],latmax=90)
        pl.title('day %s to %s' % (0.25*(i-40),0.25*(i-40)+2))
        pl.savefig('figs/lagged/T/%s/%s.pdf' % (case,i+1), format='pdf')
        pl.close()
	pl.figure(3)
	cf0  = pl.contourf(dsu.lon,dsu.lev,U0[i:i+9].mean(axis=0)-U0[0],cseq1,cmap=pl.cm.RdBu_r,extend='both')
	cl0  = pl.contour(dsu.lon,dsu.lev,U0[i:i+9].mean(axis=0),[0],colors='k',linewidths=1.5)
	cbar = pl.colorbar(cf0)
	cbar.set_label(r'Meridional-mean zonal wind speed [m s$^{-1}$]')
	pl.ylim(1000,30)
	pl.title('day %s to %s' % (0.25*(i-40),0.25*(i-40)+2))
	pl.savefig('figs/lagged/U/%s/%s.pdf' % (case,i+1), format='pdf')
	pl.close()
# Plot clim
pl.figure(1)
cseq  = np.arange(-2.5,2.5+0.25,0.25)
xx,yy = np.ma.masked_where(Uclim[0,:,:]>0,proj.x),np.ma.masked_where(Uclim[0,:,:]>0,proj.y)
cf    = pl.contourf(proj.x,proj.y,dUdZclim,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar  = pl.colorbar(cf)
pl.plot(xx[::2,::2],yy[::2,::2],'k.',alpha=0.5)
cbar.set_label(r'$\frac{U{_{500}}-U{_{925}}}{Z{_{500}}-Z{_{925}}}$ [m s$^{-1}$ km$^{-1}$]')
proj.m.drawcoastlines()
proj.m.drawparallels([70,80],latmax=90)
pl.title('Climatology')
pl.savefig('figs/lagged/baroclin/%s/clim.pdf' % (case), format='pdf')
pl.close()
pl.figure(2)
cseq  = np.arange(244,274+2,2)
cf    = pl.contourf(proj.x,proj.y,Tclim,cseq,cmap=pl.cm.OrRd,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Temperature 925hPa [K]')
proj.m.drawcoastlines()
proj.m.drawparallels([70,80,85],latmax=90)
pl.title('Climatology')
pl.savefig('figs/lagged/T/%s/clim.pdf' % (case), format='pdf')
pl.close()
