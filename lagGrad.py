from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from cmipDataServer import DataServer as cmipDataServer
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from LambertProjector import *
from scipy import interpolate
from scipy import ndimage as nd
from stipling import *

import glob
import numpy as np
import matplotlib.pyplot as pl

def extendField(field):
        fieldend = field[:,:,-1][:,:,np.newaxis]
        field    = np.append(field,fieldend,axis=2)
        field    = np.append(fieldend,field,axis=2)
        return field

def interpolateND(field,xold,xnew,axis):
        # field in (time)x(lev)x(lon)
        f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=False,fill_value='extrapolate')
        field = f(xnew)
        return field

def fill(data,invalid=None):
         if invalid is None: invalid = np.isnan(data)
         ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
         return data[tuple(ind)]

def reGrid(F,levold,latold,lonold,levnew,latnew,lonnew):
        # Extend lon axis
        lonres = np.abs(lonold[1]-lonold[0])
        dlon   = np.append(lonold,lonold[-1]+lonres)
        dlon   = np.append(lonold[0]-lonres,dlon)
	if (type(F)==np.ma.core.MaskedArray) and (F.mask.shape == F.data.shape):
	        F = fill(F.data,invalid=F.mask)
	F   = extendField(F)
	F   = interpolateND(F,  dlon,lonnew,axis=2)
	F   = interpolateND(F,latold,latnew,axis=1)
	F   = interpolateND(F,levold,levnew,axis=0)
	return F

years = range(1981,2005+1,1)
# Lambert projector, interpolation axes and other attributes
proj  = LambertProjector(boundinglat=40,resolution=200.)
lats1 = np.arange(10,30+1,1)
lats2 = np.arange(70,90+1,1)
lons  = np.arange(0,359+1,1)
levs  = np.arange(150,400+50,50)
years = range(1981,2005+1,1)
# Cosine weights
cosw1  = np.cos(np.pi*lats1/180)
cosw1  = cosw1/cosw1.sum()
cosw2  = np.cos(np.pi*lats2/180)
cosw2  = cosw2/cosw2.sum()
# Mass of girdbox on each level
dP_ = []
dP  = np.diff(levs)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665
dM  = dM/dM.sum()

# ERAInt DataServers
dTERA  = MMDataServer(Field='T',LevType='plev')
dvERA  = MMDataServer(Field='V',LevType='plev')
# Temperature data
sera1   = np.array([reGrid(dTERA.getSeason(year,Season='DJF'),dTERA.lev,dTERA.lat,dTERA.lon,levs,lats1,lons) for year in years]) - 273.15
sera2   = np.array([reGrid(dTERA.getSeason(year,Season='DJF'),dTERA.lev,dTERA.lat,dTERA.lon,levs,lats2,lons) for year in years]) - 273.15
sera1   = (sera1*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw1[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
sera2   = (sera2*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw2[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
sera    = sera1-sera2
m,c     = np.polyfit(years,sera,1)
lineE   = np.array([m*x + c for x in years])
smE,ssE = sera.mean(),sera.std()
years1  = np.array(years)[np.where(sera>(lineE+ssE))]
years2  = np.array(years)[np.where(sera<(lineE-ssE))]

#pl.plot(years,sera,'r')
#pl.plot(years,lineE,'k--')
#pl.plot(years,lineE+ssE,'k--')
#pl.plot(years,lineE-ssE,'k--')
#pl.title('ERAInt')
#pl.show()

# V field
s1 = np.array([reGrid(dvERA.getSeason(year,Season='DJF'),dvERA.lev,dvERA.lat,dvERA.lon,[500],np.arange(1,90+1,1),lons) for year in years1]).mean(axis=0)
s2 = np.array([reGrid(dvERA.getSeason(year,Season='DJF'),dvERA.lev,dvERA.lat,dvERA.lon,[500],np.arange(1,90+1,1),lons) for year in years2]).mean(axis=0)
# Stationary eddy
s1za  = s1 - s1.mean(axis=-1)[:,:,np.newaxis]
s2za  = s2 - s2.mean(axis=-1)[:,:,np.newaxis]
dsERA = s1za - s2za
dsERA = proj(dsERA,lons,np.arange(1,90+1,1)).squeeze()

# All models
#Models = [g[9:] for g in glob.glob('../rcp85/*')]
#Models = ['MRI-ESM1','bcc-csm1-1','bcc-csm1-1-m','BNU-ESM','CMCC-CESM','IPSL-CM5A-LR','IPSL-CM5A-MR']
Models  = ['MIROC5','CNRM-CM5','incmcm4','ACCESS1-0','CanESM2','GFDL-ESM2M','GFDL-ESM2G','MPI-ESM-LR','MPI-ESM-MR']
dS      = []
for Model in Models:
	try:
		print Model
		# CMIP5 DataServers
		dTCMIP = cmipDataServer(Field='ta',LevType='plev',Source=Model,DataFreq='mon',ExpType='rcp85')
		dvCMIP = cmipDataServer(Field='va',LevType='plev',Source=Model,DataFreq='mon',ExpType='rcp85')
		# Temperature data
		scmip1 = np.array([reGrid(dTCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dTCMIP.lev,dTCMIP.lat,dTCMIP.lon,levs,lats1,lons) for year in years]) - 273.15
		scmip2 = np.array([reGrid(dTCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dTCMIP.lev,dTCMIP.lat,dTCMIP.lon,levs,lats2,lons) for year in years]) - 273.15
		scmip1 = (scmip1*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw1[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
		scmip2 = (scmip2*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw2[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
		scmip  = scmip1-scmip2
		m,c    = np.polyfit(years,scmip,1)
		line   = np.array([m*x + c for x in years])
		sm,ss  = scmip.mean(),scmip.std()
		years1 = np.array(years)[np.where(scmip>(line+ss))]
		years2 = np.array(years)[np.where(scmip<(line-ss))]

		#pl.plot(years,scmip,'r')
		#pl.plot(years,line,'k--')
		#pl.plot(years,line+ss,'k--')
		#pl.plot(years,line-ss,'k--')
		#pl.title(Model)
		#pl.show()

		# V field
		s1 = np.array([reGrid(dvCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dvCMIP.lev,dvCMIP.lat,dvCMIP.lon,[500],np.arange(1,90+1,1),lons) for year in years1]).mean(axis=0).squeeze()
		s2 = np.array([reGrid(dvCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dvCMIP.lev,dvCMIP.lat,dvCMIP.lon,[500],np.arange(1,90+1,1),lons) for year in years2]).mean(axis=0).squeeze()
		# Stationary eddy
		s1za = s1 - s1.mean(axis=-1)[:,np.newaxis]
		s2za = s2 - s2.mean(axis=-1)[:,np.newaxis]
		ds   = s1za - s2za
		ds   = proj(ds,lons,np.arange(1,90+1,1))
		dS.append(ds)
	except:
		pass
N           = len(dS)
dS          = np.array(dS)
dSm         = dS.mean(axis=0)
stipx,stipy = stipling(dS,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)

pl.figure(1)
cf   = pl.contourf(proj.x,proj.y,dsERA,np.arange(-6,6+1,1),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar()
cbar.set_label('%s [%s]' % (dvERA.long_name,dvERA.units))
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('ERAInt')
pl.savefig('figs/eddy/stationary/ERAInt_trop.pdf', format='pdf')

pl.figure(2)
cf   = pl.contourf(proj.x,proj.y,dSm,np.arange(-3,3+0.5,0.5),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar()
pl.plot(stipx[::1,::1],stipy[::1,::1],'k.',alpha=0.5)
cbar.set_label('%s [%s]' % (dvCMIP.long_name,dvCMIP.units))
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5; %s models' % (N))
pl.savefig('figs/eddy/stationary/all_trop.pdf', format='pdf')

pl.figure(3)
cf   = pl.contourf(proj.x,proj.y,dSm-dsERA,12,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar()
cbar.set_label('%s [%s]' % (dvCMIP.long_name,dvCMIP.units))
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5-ERAInt')
pl.savefig('figs/eddy/stationary/bias_trop.pdf', format='pdf')

pl.show()
