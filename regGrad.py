from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from cmipDataServer import DataServer as cmipDataServer
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from LambertProjector import *
from scipy import interpolate,stats
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

# Lambert projector, interpolation axes and other attributes
proj  = LambertProjector(boundinglat=40,resolution=200.)
lats1 = np.arange(-15,15+1,1)
lats2 = np.arange(70,90+1,1)
lons  = np.arange(0,359+1,1)
levs  = np.arange(150,400+50,50)
years = np.arange(1852,2014+1,1)
# Cosine weights
cosw1 = np.cos(np.pi*lats1/180)
cosw1 = cosw1/cosw1.sum()
cosw2 = np.cos(np.pi*lats2/180)
cosw2 = cosw2/cosw2.sum()
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
sera    = sera2#-sera2
m,c,r_value,p_value,std_err = stats.linregress(years,sera)
#m,c     = np.polyfit(years,sera,1)
lineE   = np.array([m*x + c for x in years])
sera    = sera - lineE
smE,ssE = sera.mean(),sera.std()

pl.plot(years,sera+lineE,'k.-',linewidth=2,alpha=0.75)
pl.plot(years,lineE,'k--',linewidth=1,alpha=0.5)
pl.grid()
pl.xlabel('Years')
pl.ylabel('Temperature [C]')
pl.title('p = %s' % (p_value/2.))
pl.show()

# V field
s   = proj(np.array([reGrid(dvERA.getSeason(year,Season='DJF'),dvERA.lev,dvERA.lat,dvERA.lon,[500],np.arange(1,90+1,1),lons) for year in years]).squeeze(),lons,np.arange(1,90+1,1))
# Detrend
mERA,pERA = np.zeros((proj.nx,proj.ny)),np.zeros((proj.nx,proj.ny))
for i in range(proj.nx):
	for j in range(proj.ny):
		slope, intercept, r_value, p_value, std_err = stats.linregress(sera,s[:,i,j])
		mERA[i,j] = slope
		pERA[i,j] = p_value/2
stxERA,styERA = np.ma.masked_where(pERA>0.1,proj.x),np.ma.masked_where(pERA>0.1,proj.y)



# Plot
pl.figure(1)
cseq = np.arange(-0.6,0.6+0.1,0.1)
#cseq = np.arange(-2,2+0.25,0.25)
cf   = pl.contourf(proj.x,proj.y,mERA,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
pl.plot(stxERA,styERA,'k.',alpha=0.5)
cbar.set_label('dV$^{*}$/dT [m s$^{-1}$ K$^{-1}$]')
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/ERAInt.full.reg.pdf',format='pdf')
pl.show()


Models = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')][:]
mCMIP,pCMIP = [],[]
for Model in Models:
        try:
                # CMIP5 DataServers
                dTCMIP = cmipDataServer(Field='ta',LevType='plev',Source=Model,DataFreq='mon',ExpType='rcp85')
                dvCMIP = cmipDataServer(Field='va',LevType='plev',Source=Model,DataFreq='mon',ExpType='rcp85')
                # Temperature data
                scmip1 = np.array([reGrid(dTCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dTCMIP.lev,dTCMIP.lat,dTCMIP.lon,levs,lats1,lons) for year in years]) - 273.15
                scmip2 = np.array([reGrid(dTCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dTCMIP.lev,dTCMIP.lat,dTCMIP.lon,levs,lats2,lons) for year in years]) - 273.15
                scmip1 = (scmip1*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw1[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
                scmip2 = (scmip2*dM[np.newaxis,:,np.newaxis,np.newaxis]*cosw2[np.newaxis,np.newaxis,:,np.newaxis]).mean(axis=-1).sum(axis=-1).sum(axis=-1)
                scmip  = scmip2#-scmip2
                m,c    = np.polyfit(years,scmip,1)
                line   = np.array([m*x + c for x in years])
		scmip  = scmip - line
                sm,ss  = scmip.mean(),scmip.std()
		# V field
		s = proj(np.array([reGrid(dvCMIP.getDataSnaps(year,Season='DJF').mean(axis=0),dvCMIP.lev,dvCMIP.lat,dvCMIP.lon,[500],np.arange(1,90+1,1),lons) for year in years]).squeeze(),lons,np.arange(1,90+1,1))
		# Detrend
		m,p = np.zeros((proj.nx,proj.ny)),np.zeros((proj.nx,proj.ny))
		for i in range(proj.nx):
		        for j in range(proj.ny):
        		        slope, intercept, r_value, p_value, std_err = stats.linregress(scmip,s[:,i,j])
        		        m[i,j] = slope
        		        p[i,j] = p_value/2
		mCMIP.append(m)
		pCMIP.append(p)
	except:
		pass
n               = len(mCMIP)
mCMIP,pCMIP     = np.array(mCMIP),np.array(pCMIP)
dm              = mCMIP - mERA[np.newaxis,:,:]
mCMIPm,pCMIPm   = mCMIP.mean(axis=0),pCMIP.mean(axis=0)
stxCMIP,styCMIP = np.ma.masked_where(pCMIPm>0.2,proj.x),np.ma.masked_where(pCMIPm>0.2,proj.y)
stipx,stipy     = stipling(dm,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)

cseq = np.arange(-0.6,0.6+0.1,0.1)
#cseq = np.arange(-2,2+0.25,0.25)

# Plot
pl.figure(1)
cf   = pl.contourf(proj.x,proj.y,mERA,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
pl.plot(stxERA,styERA,'k.',alpha=0.5)
cbar.set_label('dV$^{*}$/dT [m s$^{-1}$ K$^{-1}$]')
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/ERAInt.full.reg.pdf',format='pdf')

pl.figure(2)
cf   = pl.contourf(proj.x,proj.y,mCMIPm,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
pl.plot(stxCMIP,styCMIP,'k.',alpha=0.5)
cbar.set_label(r'dV$^{*}$/dT [m s$^{-1}$ K$^{-1}$]')
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (%s models)' % (n))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.full.reg.pdf',format='pdf')

pl.figure(3)
cf   = pl.contourf(proj.x,proj.y,mCMIPm-mERA,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
pl.plot(stipx,stipy,'k.',alpha=0.5)
cbar.set_label(r'dV$^{*}$/dT [m s$^{-1}$ K$^{-1}$]')
drawCoastlinesNoRivers(proj.m,color='k',linewidth=1)
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5-ERAInt (%s models)' % (n))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.diff.reg.pdf',format='pdf')

pl.show()
