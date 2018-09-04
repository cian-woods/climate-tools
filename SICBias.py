from netCDF4 import Dataset as Dataset
from mpl_toolkits.basemap import Basemap as Basemap
from LambertProjector import *
from scipy.interpolate import griddata
from stipling import stipling

import numpy as np
import glob
import matplotlib.pyplot as pl
import sys

# ------------------------ Input ------------------------ #
Field      = 'sic'
Projection = 'nplaea'
blat       = 55
proj       = LambertProjector(boundinglat=blat,resolution=80.,Projection=Projection)
bm         = Basemap(boundinglat=blat,lon_0=0,projection=Projection)
notM       = ['BNU-ESM','IPSL-CM5A-LR','bcc-csm1','bcc-csm1-m','MIROC5']
Models     = [g[9:] for g in glob.glob('../rcp85/*') if g[9:] not in notM]
# ------------------------------------------------------- #

# CMIP5
C1,C2,cseq = [],[],np.arange(-80,80+10,10)
for Model in Models[0:2]:
	print Model
	Dir      = '/mnt/climstorage/cian/rcp85/%s/mon/surface/%s' % (Model,Field)
	Handles  = [Dataset(Name) for Name in sorted(glob.glob('%s/*.nc' % (Dir)))]
	Year0s   = [int(Handle.variables['time'].units.split()[2].split('-')[0]) for Handle in Handles]
	lon,lat  = Handles[0].variables['lon'][:],Handles[0].variables['lat'][:]
	fillval = Handles[0].variables[Field]._FillValue
	missval = Handles[0].variables[Field].missing_value
	s        = np.concatenate(np.array([Handle.variables[Field][:] for Handle in Handles]))
	n0,n1,n2 = s.shape
	s        = s.reshape(n0/12.,12,n1,n2)
	s        = s[-25:,[0,10,11],:,:]
	s[np.where(s==fillval)] = 0
	s[np.where(s==missval)] = 0
	s[np.where(s<0)] = 0
	s[np.where(s>100)] = 100
	s = s.mean(axis=0).mean(axis=0)
	if len(lon.shape)==2: j,k = bm(lon,lat)
	if len(lon.shape)==1: j,k = bm(*np.meshgrid(lon,lat))
	data  = griddata((j.reshape(-1),k.reshape(-1)), s.reshape(-1), (proj.x.reshape(-1), proj.y.reshape(-1)), method='linear')
	data1 = data.reshape((proj.nx,proj.ny))
	C1.append(data1)

        Dir      = '/mnt/oldstorage/cian/cmip/historical/%s/mon/surface/%s' % (Model,Field)
        Handles  = [Dataset(Name) for Name in sorted(glob.glob('%s/*.nc' % (Dir)))]
        Year0s   = [int(Handle.variables['time'].units.split()[2].split('-')[0]) for Handle in Handles]
        lon,lat  = Handles[0].variables['lon'][:],Handles[0].variables['lat'][:]
        fillval = Handles[0].variables[Field]._FillValue
        missval = Handles[0].variables[Field].missing_value
        s        = np.concatenate(np.array([Handle.variables[Field][:] for Handle in Handles]))
        n0,n1,n2 = s.shape
        s        = s.reshape(n0/12.,12,n1,n2)
        s        = s[-25:,[0,10,11],:,:]
        s[np.where(s==fillval)] = 0
        s[np.where(s==missval)] = 0
        s[np.where(s<0)] = 0
        s[np.where(s>100)] = 100
        s = s.mean(axis=0).mean(axis=0)
        if len(lon.shape)==2: j,k = bm(lon,lat)
        if len(lon.shape)==1: j,k = bm(*np.meshgrid(lon,lat))
        data  = griddata((j.reshape(-1),k.reshape(-1)), s.reshape(-1), (proj.x.reshape(-1), proj.y.reshape(-1)), method='linear')
        data2 = data.reshape((proj.nx,proj.ny))
        C2.append(data2)

	# Plot
	cf   = pl.contourf(proj.x,proj.y,data1-data2,cseq,cmap=pl.cm.RdBu_r,extend='both')
	cl1  = pl.contour(proj.x,proj.y,data1,levels=[15],linewidths=2,colors='k',linestyles='dashed',alpha=0.8)
	cl2  = pl.contour(proj.x,proj.y,data2,levels=[15],linewidths=2,colors='k',alpha=0.8)
	cbar = pl.colorbar(cf)
	cbar.set_label('Sea ice concentration change [%]')
	bm.drawcoastlines(color='0.4',linewidth=0.7)
	bm.drawparallels([70,80],latmax=90)
	bm.fillcontinents(color='0.8', lake_color='0.8', zorder=None)
	pl.title(Model)
	pl.savefig('figs/bias/%s/%s.%s.pdf' % (Field,Model,Field),format='pdf')
	pl.close()

C1,C2       = np.array(C1),np.array(C2)
C           = C1-C2
stipx,stipy = stipling(C,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
C           = C.mean(axis=0)
cseq        = np.arange(-50,50+10,10)
C           = np.ma.masked_where(C==0,C)
cf          = pl.contourf(proj.x,proj.y,C,cseq,cmap=pl.cm.RdBu_r,extend='both')
cl1         = pl.contour(proj.x,proj.y,C1.mean(axis=0),levels=[15],linewidths=2,colors='k',linestyles='dashed',alpha=0.8)
cl2         = pl.contour(proj.x,proj.y,C2.mean(axis=0),levels=[15],linewidths=2,colors='k',alpha=0.8)
pl.plot(stipx[::2,::2],stipy[::2,::2],'k.',alpha=0.5)
cbar = pl.colorbar(cf)
cbar.set_label('Sea ice concentration change [%]')
bm.drawcoastlines(color='0.4')
bm.drawparallels([70,80],latmax=90)
bm.fillcontinents(color='0.8', lake_color='0.8', zorder=None)
pl.title('all')
pl.savefig('figs/bias/%s/%s.%s.pdf' % (Field,'all',Field),format='pdf')
pl.close()

