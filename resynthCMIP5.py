from netCDF4 import Dataset
from LambertProjector import *
from drawCoastlinesNoRivers import *
from stipling import stipling

import glob
import numpy as np
import matplotlib.pyplot as pl

def moving_average(a,n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:]/n

wn0,wn1    = int(sys.argv[1]),int(sys.argv[2])
LevRange   = (600,1000)
wvsec      = np.arange(wn0,wn1+1,1) - 1
# LambertProjector
proj = LambertProjector(boundinglat=60,resolution=100.)
# Get axes from a given file (same for all files)
f         = Dataset('/mnt/climstorage/cian/%s_resynth/%s_resynth_1979.nc' % ('v','v'),'r')
long_name = f.variables['v'].long_name
units     = f.variables['v'].units
levs      = f.variables['lev'][:]
lats      = f.variables['lat'][:]
lons      = f.variables['lon'][:]
wvns      = f.variables['wavenumber'][:]
f.close()
# Indices
levx1,levx2         = np.argmin((levs-LevRange[0])**2),np.argmin((levs-LevRange[1])**2)
latx                = np.argmin((lats-70)**2)
levs                = levs[levx1:levx2+1]
wvns                = wvns[wvsec]
nlev,nlon,nlat,nwvn = len(levs),len(lons),len(lats),len(wvns)
# Mass of girdbox on each level
dP_ = []
dP  = np.diff(levs)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665
# Years
years = range(1980,2005+1,1)

# Get ERAInt data
print 'ERAInt'
SERAv,SERAq,SERAvkqk  = [],[],[]
for year in years:
	fv    = Dataset('/mnt/climstorage/cian/%s_resynth/%s_resynth_%s.nc' % ('v','v',year),'r')
	fq    = Dataset('/mnt/climstorage/cian/%s_resynth/%s_resynth_%s.nc' % ('q','q',year),'r')
	fvkqk = Dataset('/mnt/climstorage/cian/%s_resynth/%s_resynth_%s.nc' % ('vkqk','vkqk',year),'r')
	v     =    fv.variables[ 'v'][[0,1,11],levx1:levx2+1,:,:,wvsec]
	q     =    fq.variables[ 'q'][[0,1,11],levx1:levx2+1,:,:,wvsec]
	vkqk  = fvkqk.variables['vq'][[0,1,11],levx1:levx2+1,:,:,wvsec]
	SERAv.append(v)
	SERAq.append(q)
	SERAvkqk.append(vkqk)
	fv.close()
	fq.close()
	fvkqk.close()
SERAv,SERAq = np.array(SERAv),np.array(SERAq)
SERAvkqk    = np.array(SERAvkqk)
SERAv       =    SERAv.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
SERAq       =    SERAq.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
SERAvkqk    = SERAvkqk.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
# Time mean and vertical integral, summed over all wavenumbers
SERAv    = ((SERAv[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
SERAv    = proj(SERAv,lons,lats)
SERAq    = ((SERAq[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
SERAq    = proj(SERAq,lons,lats)
SERAvkqk = ((SERAvkqk[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
SERAvkqk = proj(SERAvkqk,lons,lats)

# Get CMIP5 data
Models        = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')]
SCMIPv,SCMIPq = [],[]
SCMIPvkqk     = []
for Model in Models:
	print Model
	Scmipv,Scmipq,Scmipvkqk = [],[],[]
	for year in years:
	        fv    = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/plev/%s_resynth/%s_resynth_%s.nc' % (Model,'v','v',year),'r')
		fq    = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/plev/%s_resynth/%s_resynth_%s.nc' % (Model,'q','q',year),'r')
		fvkqk = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/plev/%s_resynth/%s_resynth_%s.nc' % (Model,'vkqk','vkqk',year),'r')
	        v     = fv.variables['v'][[0,1,11],levx1:levx2+1,:,:,wvsec]
		q     = fq.variables['q'][[0,1,11],levx1:levx2+1,:,:,wvsec]
		vkqk  = fvkqk.variables['vq'][[0,1,11],levx1:levx2+1,:,:,wvsec]
	        Scmipv.append(v)
		Scmipq.append(q)
		Scmipvkqk.append(vkqk)
	        fv.close()
		fq.close()
	Scmipv,Scmipq = np.array(Scmipv),np.array(Scmipq)
	Scmipvkqk     = np.array(Scmipvkqk)
	Scmipv        = Scmipv.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
	Scmipq        = Scmipq.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
	Scmipvkqk     = Scmipvkqk.reshape((-1,nlev,nlat,nlon,nwvn))[2:-1].reshape((-1,3,nlev,nlat,nlon,nwvn)).mean(axis=1)
	# Time mean and vertical integral, summed over all wavenumbers
	Scmipv    = ((Scmipv[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
	Scmipv    = proj(Scmipv,lons,lats)
        Scmipq    = ((Scmipq[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
        Scmipq    = proj(Scmipq,lons,lats)
        Scmipvkqk = ((Scmipvkqk[:,:,:,:,:].mean(axis=0).sum(axis=-1))*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
        Scmipvkqk = proj(Scmipvkqk,lons,lats)
	SCMIPv.append(Scmipv)
	SCMIPq.append(Scmipq)
	SCMIPvkqk.append(Scmipvkqk)
SCMIPv,SCMIPq   = np.array(SCMIPv),np.array(SCMIPq)
SCMIPvkqk       = np.array(SCMIPvkqk)
SCMIPvm,SCMIPqm = SCMIPv.mean(axis=0),SCMIPq.mean(axis=0)
SCMIPvkqkm      = SCMIPvkqk.mean(axis=0)
# Stiplings and mask where both v and q are same sign
stipx,stipy     = stipling(SCMIPv-SERAv,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
maskERA         = np.sign(SERAv)==np.sign(SERAq)
stxERA,styERA   = np.ma.masked_array(proj.x,mask=maskERA),np.ma.masked_array(proj.y,mask=maskERA)
maskCMIP        = np.sign(SCMIPvm)==np.sign(SCMIPqm)
stxCMIP,styCMIP = np.ma.masked_array(proj.x,mask=maskCMIP),np.ma.masked_array(proj.y,mask=maskCMIP)

# Plot
pl.figure(1)
cf   = pl.contourf(proj.x,proj.y,SERAv/1000,14,cmap=pl.cm.RdBu_r,extend='both')
cl   = pl.contour(proj.x,proj.y,SERAq,14,colors='0.65',linewidths=1.2)
cbar = pl.colorbar(cf)
pl.plot(stxERA[::1,::1],styERA[::1,::1],'k.',alpha=0.3)
pl.clabel(cl,fmt='%0.2f',colors='k',fontsize=9)
cbar.set_label('Meridional mass flux [10$^{3}$ kg s$^{-1}$ m$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.ERAInt.full.%d-%dhPa.%d-%dk.pdf' % ('v.q',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(2)
cf   = pl.contourf(proj.x,proj.y,SERAv*SERAq/1000,np.arange(-8,8+1,1),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('vq [10$^{3}$ kg s$^{-1}$ m$^{-1}$ kg m$^{-2}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.ERAInt.full.%d-%dhPa.%d-%dk.pdf' % ('vk.qk',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(3)
cf   = pl.contourf(proj.x,proj.y,SERAvkqk,np.arange(-3,3+0.5,0.5),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('Meridional moisture flux [kg s$^{-1}$ m$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.ERAInt.full.%d-%dhPa.%d-%dk.pdf' % ('vkqk',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(4)
cf   = pl.contourf(proj.x,proj.y,SCMIPvm/1000,14,cmap=pl.cm.RdBu_r,extend='both')
cl   =  pl.contour(proj.x,proj.y,SCMIPqm,14,colors='0.65',linewidths=1.2)
cbar = pl.colorbar(cf)
pl.plot(stxCMIP[::1,::1],styCMIP[::1,::1],'k.',alpha=0.3)
pl.clabel(cl,fmt='%0.2f',colors='k',fontsize=9)
cbar.set_label('Meridional mass flux [10$^{3}$ kg s$^{-1}$ m$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('CMIP5')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.all.full.%d-%dhPa.%d-%dk.pdf' % ('v.q',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(5)
cf   = pl.contourf(proj.x,proj.y,SCMIPvm*SCMIPqm/1000,np.arange(-8,8+1,1),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('vq [10$^{3}$ kg s$^{-1}$ m$^{-1}$ kg m$^{-2}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('CMIP5')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.all.full.%d-%dhPa.%d-%dk.pdf' % ('vk.qk',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(6)
cf   = pl.contourf(proj.x,proj.y,SCMIPvkqkm - SERAvkqk,np.arange(-1.5,1.5+0.25,0.25),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('Meridional moisture flux [kg s$^{-1}$ m$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('CMIP5-ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.all.bias.%d-%dhPa.%d-%dk.pdf' % ('vkqk',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(7)
cf   = pl.contourf(proj.x,proj.y,((SCMIPvm*SCMIPqm)-(SERAv*SERAq))/1000,np.arange(-3,3+0.5,0.5),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('vq [10$^{3}$ kg s$^{-1}$ m$^{-1}$ kg m$^{-2}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('CMIP5-ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.all.diff.%d-%dhPa.%d-%dk.pdf' % ('vk.qk',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.figure(8)
cf   = pl.contourf(proj.x,proj.y,(SCMIPvm-SERAv)/1000,14,cmap=pl.cm.RdBu_r,extend='both')
cl   =  pl.contour(proj.x,proj.y,SCMIPqm-SERAq,8,colors='0.65',linewidths=1.2)
cbar = pl.colorbar(cf)
pl.clabel(cl,fmt='%0.2f',colors='k',fontsize=9)
cbar.set_label('Meridional mass flux [10$^{3}$ kg s$^{-1}$ m$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmaz=90)
pl.title('CMIP5 - ERAInt')
pl.savefig('/mnt/climstorage/cian/scripts/figs/wavenumber/resynth/%s.all.diff.%d-%dhPa.%d-%dk.pdf' % ('v.q',levs[0],levs[-1],wvns[0],wvns[-1]),format='pdf')

pl.show()
