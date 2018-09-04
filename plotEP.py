from netCDF4 import *

import matplotlib.pyplot as pl
import numpy as np
import glob

Files  = glob.glob('/mnt/climstorage/cian/EPFD/EPFD_????_01.nc')[:]
d      = [Dataset(File,'r') for File in Files]

sid    = 24*3600
Fphi   = np.array([d[ii].variables['Fphi'][:,:,0:-1].mean(axis=0) for ii in range(len(d))])
Fp     = np.array([d[ii].variables['Fp'][:,:,0:-1].mean(axis=0)   for ii in range(len(d))])
div    = np.array([d[ii].variables['Fdiv'][:,:,0:-1].mean(axis=0) for ii in range(len(d))])

lat    = d[0].variables['lat'][0:-1]
lev    = d[0].variables['lev'][:]
loglev = np.exp(np.linspace(np.log(10),np.log(500),20))
xs     = [np.argmin((lev-loglev[ii])**2) for ii in range(len(loglev))]
N      = np.ones(len(lev))
N[np.where(lev<100)] = 3
N      = N[:,np.newaxis]

Fphi0 = Fphi[:,:,:].mean(axis=0)*N
Fp0   = -1*Fp[:,:,:].mean(axis=0)*N
div0  = -1*sid*div[:,:,:].mean(axis=0)

cseq = np.arange(-25,25+5,5)
cf   = pl.contourf(lat,lev,div0,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar()
cbar.set_label('Zonal wind acceleration [m s$^{-1}$ day$^{-1}$]')

Q  = pl.quiver(lat[::5],lev[xs],Fphi0[xs,::5],Fp0[xs,::5],units='inches',scale=500,\
		scale_units='inches',headwidth=2.5,headlength=2.5,headaxislength=2.5,pivot='tail',alpha=0.6)
qk = pl.quiverkey(Q, 0.2, 1.02, 300, '%s%s' % (300,''), labelpos='W',fontproperties={'weight': 'bold'})

pl.ylabel('Pressure [hPa]')
pl.xlabel(r'Latitude [$^{\circ}$ N]')
pl.ylim(1000,10)
pl.yscale('log')
pl.xlim(-90,90)
pl.title('E-P divergence')
pl.show()
