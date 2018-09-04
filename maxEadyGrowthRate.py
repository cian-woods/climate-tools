import numpy as np

def brunt_vaisala(T,Z,P):
	# Compute brunt-viasala frequency N
	# N = sqrt((g/<theta>)*(dtheta/dP)
	P0P_k  = (1000/P)**0.286
	theta  = T*(P0P_k[:,np.newaxis,np.newaxis])
	# Column average potential temperature (weighted mean using pressure)
	theta0 = (theta*P[:,np.newaxis,np.newaxis]).sum(axis=0)/P.sum()	# K
	# Calculate dtheta/dP
	dtheta   = -1*np.diff(theta[[0,-1],:,:],axis=0).squeeze()
	dz       = -1*np.diff(Z[[0,-1],:,:],axis=0).squeeze()/9.81
	dthetadz = dtheta/dz	# K/m
	# Calculate N
	N = np.sqrt((9.81/theta0)*dthetadz)
	return N

def dudz(U,Z):
	dU   = -1*np.diff(U[[0,-1],:,:],axis=0).squeeze()
	dz   = -1*np.diff(Z[[0,-1],:,:],axis=0).squeeze()/9.81
	dUdz = dU/dz
	return dUdz

def megr(U,Z,T,P,lat):
	# Compute theoretical maximum eady growtn rate
	# gr = 0.31*g*f*(dUdz)/N
	N       = brunt_vaisala(T,Z,P)
	dUdz    = dudz(U,Z)
	f       = (2*(7.29e-5)*np.sin(lat*np.pi/180))[:,np.newaxis]
	gr      = (0.31*9.81*f*dUdz/N)
	return gr,N,dUdz,f

if __name__ == "__main__":

	import sys
	sys.path.insert(0, '/home/cian/scripts')
	from ReanalysisDataServer import *
	from LambertProjector import *

	import matplotlib.pyplot as pl

	proj = LambertProjector(boundinglat=60,resolution=100.)
	dsu  = DataServer(Field='U',LevType='plev',LevRange=(500,925),LatRange=(30,90))
	dsz  = DataServer(Field='Z',LevType='plev',LevRange=(500,925),LatRange=(30,90))
	dst  = DataServer(Field='T',LevType='plev',LevRange=(500,925),LatRange=(30,90))
	latx = np.argmin((dst.lat-70)**2)

	for date in dsu.getDateList(Year=2016,Month=12):
		u,z,t       = dsu.snapshot(*date),dsz.snapshot(*date),dst.snapshot(*date)
		gr,N,dUdz,f = megr(u,z,t,dsu.lev,dsu.lat)
		cseq  = np.arange(-4,4+0.5,0.5)
		cseq0 = np.arange(255,275+2,2)
		gr    = proj(1000*dUdz,dsu.lon,dsu.lat)
		t0    = proj(t[-1,latx:,:],dst.lon,dst.lat[latx:])
		t0    = np.ma.masked_where(proj.lat<70,t0)
		cf    = pl.contourf(proj.x,proj.y,gr,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cl    = pl.contour(proj.x,proj.y,t0,cseq0,colors='k',linewidths=1,alpha=0.5)
		pl.clabel(cl,fmt='%3.0f',colors='k',fontsize=8)
		cbar  = pl.colorbar(cf)
		cbar.set_label('dU/dz [m s$^{-1}$ km$^{-1}$]')
		proj.m.drawcoastlines(color='0.5',linewidth=0.7)
		proj.m.drawparallels([70,80],latmax=90)
		pl.title(str(date))
		pl.show()
