from netCDF4 import *
from LambertProjector import *
from drawCoastlinesNoRivers import *
from scipy import stats,interpolate
from ReanalysisMonthlyMeanDataServer import *
from toPick import *
from UnPickle import *

import sys
import numpy as np
import matplotlib.pyplot as pl

SMs = {}
SMs['D']	 = [12]
SMs['J']         = [1]
SMs['DJF']       = [12,1,2]
SMs['NDJF']      = [11,12,1,2]
SMs['JJA']       = [6,7,8]
SMs['SON']       = [9,10,11]
SMs['MAM']       = [3,4,5]
SMs['Annual']    = range(1,12+1,1)
SMs['Annual_wc'] = range(7,12+1,1) + range(1,6+1,1)
SMs['Annual_MF'] = range(3,12+1,1) + range(1,2+1,1)
SMs['MAMJJASO']  = [3,4,5,6,7,8,9,10]
SMs['MAMJJASO_l1']  = [3,4,5,6,7,8,9,10]

def interpolateND(field,xold,xnew,axis):
        # field in (time)x(lev)x(lon)
        #f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
        f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=False,fill_value='extrapolate')
        field = f(xnew)
        return field

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

def detrend2d(years,field):
        n0,n1,n2 = field.shape
        m,c,p    = np.zeros((n1,n2)),np.zeros((n1,n2)),np.zeros((n1,n2))
        for i in range(n1):
                for j in range(n2):
                        slope, intercept, r_value, p_value, std_err = stats.linregress(years,field[:,i,j])
                        m[i,j] = slope
                        c[i,j] = intercept
                        p[i,j] = p_value
        line  = np.array([m*year + c for year in years])
        field = field - line
        return field,m,p

def xres(lat):
	xr = 6370000*2*np.pi*np.cos(np.pi*lat/180.)/360.
	print xr
	return xr

def bathy(lons,lats):
	bd         = Dataset('/mnt/climstorage/cian/Bathymetry_etopo2-cian.nc','r')
	bath_lon   = bd.variables['lon'][:]
	bath_lat   = bd.variables['lat'][:]
	bathy      = -1*bd.variables['Bathymetry'][:].squeeze()
	bathy      = interpolateND(bathy,bath_lon,lons,axis=1)
	bathy      = interpolateND(bathy,bath_lat,lats,axis=0)
	return bathy

# Inputs
Season     = str(sys.argv[1])
YearRange  = (int(sys.argv[2]),int(sys.argv[3]))
Field      = str(sys.argv[4])
# Attributes
#ymoor,htmoor = unpick('/mnt/climstorage/cian/BSO/ht_moor.bso.full.%s-%s.%s.p' % (YearRange[0],YearRange[1],Season))
years        = range(1979,2016+1,1)
yres         = 6370000*2*np.pi/360.
proj         = LambertProjector(boundinglat=50,resolution=200.)
SeasMonths   = SMs[Season]
if SeasMonths[0]>SeasMonths[-1]: years = years[1:]
# Open files
du = Dataset('/mnt/climstorage/cian/ORAS4/monthly/lev/uo.1979-2016.mon.mean.deeper.nc','r')
dv = Dataset('/mnt/climstorage/cian/ORAS4/monthly/lev/vo.1979-2016.mon.mean.deeper.nc','r')
dT = Dataset('/mnt/climstorage/cian/ORAS4/monthly/lev/thetao.1979-2016.mon.mean.deeper.nc','r')
# Axes
times  = du.variables['time']
months = np.tile(np.arange(1,12+1,1),(len(times)/12.))
xs     = [range(i,i+len(SeasMonths)) for i in range(len(months)) if (list(months[i:i+len(SeasMonths)]) == SeasMonths)]
xsD    = [range(i,i+len(SMs['D'])) for i in range(len(months)) if (list(months[i:i+len(SMs['D'])]) == SMs['D'])]
xsJ    = [range(i,i+len(SMs['J'])) for i in range(len(months)) if (list(months[i:i+len(SMs['J'])]) == SMs['J'])]
lons   = du.variables['longitude'][:]
lats   = du.variables['latitude'][:]
levs   = du.variables['LEV'][:]
# LandSea mask
cids = DataServer(Field='ci' )
lsm  = cids.lsm
lsm  = interpolateND(lsm,cids.lsm_lon,lons,axis=1)
lsm  = interpolateND(lsm,cids.lsm_lat,lats,axis=0)
lsm[np.where(lsm>0)]  = 1
lsm[np.where(lsm<=0)] = 0
# Mass of girdbox on each level
dD_ = []
dD  = np.diff(levs)/2.
dD_.append(levs[0]+dD[0])
for i in range(len(dD)-1):
        dD_.append(dD[i]+dD[i+1])
dD_.append(dD[-1])#+(1000-levs[-1]))
dD_ = np.array(dD_)
dM  = 1028.*dD_
# Data (and regrid)
#sln0,sln1 = np.argmin((lons-0.5)**2),np.argmin((lons-100.5)**2)
#slt0,slt1 = np.argmin((lats-68.5)**2),np.argmin((lats-83.5)**2)

#sln0,sln1 = np.argmin((lons-20.5)**2),np.argmin((lons-35.5)**2)
#slt0,slt1 = np.argmin((lats-70.5)**2),np.argmin((lats-75.5)**2)

sln0,sln1 = np.argmin((lons-0.5)**2),np.argmin((lons-360)**2)
slt0,slt1 = np.argmin((lats-68.5)**2),np.argmin((lats-90)**2)

lons    = lons[sln0:sln1+1]
lats    = lats[slt0:slt1+1]
uo      = du.variables[    'uo'][:,:,slt0:slt1+1,sln0:sln1+1]
vo      = dv.variables[    'vo'][:,:,slt0:slt1+1,sln0:sln1+1]
uv      = np.sqrt(uo**2 + vo**2)
T       = dT.variables['thetao'][:,:,slt0:slt1+1,sln0:sln1+1] + 1# + 273.15
lsm     = lsm[slt0:slt1+1,sln0:sln1+1]
# Slices
levx0  = np.argmin((levs-0 )**2)
levx1  = np.argmin((levs-5000 )**2)
# Right
lonxr  = np.argmin((lons-90.5)**2)
latx0r = np.argmin((lats-70.5)**2)
latx1r = np.argmin((lats-81.5)**2)
# Left
lonxl  = np.argmin((lons-20.5)**2)
latx0l = np.argmin((lats-70.5)**2)
latx1l = np.argmin((lats-81.5)**2)
# Bottom
lonx0b = np.argmin((lons-20.5)**2)
lonx1b = np.argmin((lons-90.5)**2)
latxb  = np.argmin((lats-70.5)**2)
# Top
lonx0t = np.argmin((lons-20.5)**2)
lonx1t = np.argmin((lons-90.5)**2)
latxt  = np.argmin((lats-81.5)**2)

df_tmp = DataServer(Field='ci',Source='ERAInt')
ci     = np.array([df_tmp.getSeason(Year=year,Season=Season) for year in range(1980,2016+1,1)])

yx0,yx1          = years.index(YearRange[0]),years.index(YearRange[1])
Hux,Huy          = (uo*T)[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=0).mean(axis=0),(vo*T)[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=0).mean(axis=0)

ToD               = T[xsD,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=1)
ToDa              = ToD - ToD.mean(axis=0)[np.newaxis,:,:,:]
OHCD              = (3985*ToD*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)

ToD0              = T[xsD,:,:,:][yx0:yx1+1-1,:,:,:,:].mean(axis=1)
ToD0a             = ToD0 - ToD0.mean(axis=0)[np.newaxis,:,:,:]
OHCD0             = (3985*ToD0*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)
ToD1              = T[xsD,:,:,:][yx0+1:yx1+1,:,:,:,:].mean(axis=1)
ToD1a             = ToD1 - ToD1.mean(axis=0)[np.newaxis,:,:,:]
OHCD1             = (3985*ToD1*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)

ToJ               = T[xsJ,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=1)
ToJa              = ToJ - ToJ.mean(axis=0)[np.newaxis,:,:,:]
OHCJ              = (3985*ToJ*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)

deltaOHC          = OHCD1 - OHCD0

To             = T[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=1)
Toa            = To - To.mean(axis=0)[np.newaxis,:,:,:]
OHC            = (3985*Toa*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)
dOHC,mOHC,pOHC = detrend2d(years[yx0:yx1+1],OHC)

Hux_y,Huy_y    = (uo*T)[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=1),(vo*T)[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=1)
Hux_y,Huy_y     = 3985*(Hux_y*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1),3985*(Huy_y*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)

dHux,mHux,pHux = detrend2d(years[yx0:yx1+1],Hux_y)
dHuy,mHuy,pHuy = detrend2d(years[yx0:yx1+1],Huy_y)

UV          = uv[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=0).mean(axis=0).mean(axis=1).mean(axis=1)
Uux,Uuy     = uo[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=0).mean(axis=0).mean(axis=0),vo[xs,:,:,:][yx0:yx1+1,:,:,:,:].mean(axis=0).mean(axis=0).mean(axis=0)
Hux,Huy     = 3985*(Hux*dM[:,np.newaxis,np.newaxis]).sum(axis=0),3985*(Huy*dM[:,np.newaxis,np.newaxis]).sum(axis=0)

# Save climatology fields
fname_h = '/mnt/climstorage/cian/ht.clim.reanal.%s-%s.%s.p' % (YearRange[0],YearRange[1],Season)
fname_u = '/mnt/climstorage/cian/vt.clim.reanal.%s-%s.%s.p' % (YearRange[0],YearRange[1],Season)
if not os.path.isfile(fname_h):
	toPick([Hux,Huy,lats,lons],fname_h)
if not os.path.isfile(fname_u):
        toPick([Uux,Uuy,lats,lons],fname_u)

#bmap    = Basemap(projection='cyl',llcrnrlat=lats[0]-2,urcrnrlat=lats[-1]+2,llcrnrlon=lons[0]-4,urcrnrlon=lons[-1]+4,resolution='i')
bmap    = Basemap(projection='cyl',llcrnrlat=lats[0],urcrnrlat=lats[-1],llcrnrlon=lons[0],urcrnrlon=lons[-1],resolution='i')
j,k     = bmap(lons,lats)
jj,kk   = bmap(df_tmp.lon,df_tmp.lat)

fig     = pl.figure(1)
ax      = fig.add_subplot(111)
cseq    = 13#np.arange(0,0.06+0.005,0.005)
cf      = ax.contourf(j,k,np.sqrt(Hux**2 + Huy**2)/1e08,cseq,cmap=pl.cm.OrRd,extend='both')
for i in ci:
	cl = pl.contour(jj,kk,i,[0.15],colors='k',linewidths=0.3,alpha=0.4)
Q       = pl.quiver(j,k,Hux,Huy,pivot='tail',alpha=0.5)
cbar    = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Ocean heat transport [10$^{8}$ W m$^{-1}$]')
bmap.drawcoastlines(linewidth=0.5,color='0.5')
bmap.fillcontinents(color='0.85')
bmap.drawmeridians([0,20,40,60,80,100])
bmap.drawparallels([70,75,80])
ax.set_aspect('auto')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.clim.reanal.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

fig     = pl.figure(2)
ax      = fig.add_subplot(111)
cseq    = 13#np.arange(0,0.06+0.005,0.005)
cf      = ax.contourf(j,k,np.sqrt(Uux**2 + Uuy**2),cseq,cmap=pl.cm.OrRd,extend='both')
for i in ci:
        cl      = pl.contour(jj,kk,i,[0.15],colors='k',linewidths=0.3,alpha=0.4)
Q       = pl.quiver(j,k,Uux,Uuy,pivot='tail',alpha=0.5)
cbar    = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Ocean velocity [m s$^{-1}$]')
bmap.drawcoastlines(linewidth=0.5,color='0.5')
bmap.fillcontinents(color='0.85')
bmap.drawmeridians([0,20,40,60,80,100])
bmap.drawparallels([70,75,80])
ax.set_aspect('auto')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/vt.clim.reanal.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

pl.figure(3)
pl.plot(UV,levs[levx0:levx1+1],'b',linewidth=2)
pl.xlabel('Ocean velocity [m s$^{-1}$]')
pl.ylabel('Depth [m]')
pl.ylim(400,0)

pl.show()
sys.exit()

fig       = pl.figure(3)
ax        = fig.add_subplot(111)
cseq      = 13#np.arange(0,0.06+0.005,0.005)
cf        = ax.contourf(j,k,10*np.sqrt(mHux**2 + mHuy**2)/1e08,cseq,cmap=pl.cm.OrRd,extend='both')
mHux,mHuy = np.ma.masked_where(pHux>0.1,mHux),np.ma.masked_where(pHuy>0.1,mHuy)
Q         = pl.quiver(j,k,10*mHux,10*mHuy,pivot='tail',alpha=0.5)
cbar      = pl.colorbar(cf,orientation='horizontal')
pl.plot(lons[lonxl:lonxr+1],[lats[latxb] for i in range(lonxl,lonxr+1)],'b',linewidth=1.5,alpha=0.6)
pl.plot(lons[lonxl:lonxr+1],[lats[latxt] for i in range(lonxl,lonxr+1)],'b',linewidth=1.5,alpha=0.6)
pl.plot([lons[lonxl] for i in range(latxb,latxt+1)],lats[latxb:latxt+1],'b',linewidth=1.5,alpha=0.6)
pl.plot([lons[lonxr] for i in range(latxb,latxt+1)],lats[latxb:latxt+1],'b',linewidth=1.5,alpha=0.6)
cbar.set_label('Ocean heat transport trend [10$^{8}$ W m$^{-1}$ decade$^{-1}$]')
bmap.drawcoastlines(linewidth=0.5,color='0.5')
bmap.fillcontinents(color='0.85')
bmap.drawmeridians([0,20,40,60,80,100])
bmap.drawparallels([70,75,80])
ax.set_aspect('auto')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.trend.reanal.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

fig       = pl.figure(4)
ax        = fig.add_subplot(111)
cseq      = np.arange(-10,10+1.25,1.25)
#cf        = ax.contourf(j,k,10*np.ma.masked_where(lsm==1,mOHC)/1e08,cseq,cmap=pl.cm.RdBu_r,extend='both')
cf        = ax.contourf(j,k,10*mOHC/1e08,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar      = pl.colorbar(cf,orientation='horizontal')
pl.plot(lons[lonxl:lonxr+1],[lats[latxb] for i in range(lonxl,lonxr+1)],'b',linewidth=1.5,alpha=0.6)
pl.plot(lons[lonxl:lonxr+1],[lats[latxt] for i in range(lonxl,lonxr+1)],'b',linewidth=1.5,alpha=0.6)
pl.plot([lons[lonxl] for i in range(latxb,latxt+1)],lats[latxb:latxt+1],'b',linewidth=1.5,alpha=0.6)
pl.plot([lons[lonxr] for i in range(latxb,latxt+1)],lats[latxb:latxt+1],'b',linewidth=1.5,alpha=0.6)
cbar.set_label('Ocean heat content trend [10$^{8}$ J m$^{-2}$ decade$^{-1}$]')
bmap.drawcoastlines(linewidth=0.5,color='0.5')
bmap.fillcontinents(color='0.85')
bmap.drawmeridians([0,20,40,60,80,100])
bmap.drawparallels([70,75,80])
ax.set_aspect('auto')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ohc.trend.reanal.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

pl.figure(5)
OHC_w = 1.85e12*((deltaOHC[:,latxb:latxt+1,lonxl:lonxr+1].mean(axis=-1))*(np.cos(np.pi*lats[latxb:latxt+1]/180)[np.newaxis,:])).sum(axis=1)/(1e18*np.cos(np.pi*lats[latxb:latxt+1]/180.).sum())
OHC_w = np.cumsum(OHC_w)
OHC_w = OHC_w - OHC_w.mean()
m,c   = np.polyfit(years[yx0+1:yx1+1],OHC_w,1)
line  = [m*x + c for x in years[yx0+1:yx1+1]]
pl.plot(years[yx0+1:yx1+1],OHC_w,'k.-',linewidth=1.75,alpha=0.7)
pl.plot(years[yx0+1:yx1+1],line ,'k--',linewidth=1.00,alpha=0.5)
pl.xlabel('Year')
pl.ylabel('Ocean heat content anomaly [10$^{18}$ J]')
pl.grid()
pl.xlim(years[yx0],years[yx1])
pl.title('trend = %s 10$^{18}$ J year$^{-1}$' % (round(m,3)))

pl.figure(6)
OHC_w = 1.85e12*((OHC[:,latxb:latxt+1,lonxl:lonxr+1].mean(axis=-1))*(np.cos(np.pi*lats[latxb:latxt+1]/180)[np.newaxis,:])).sum(axis=1)/(1e18*np.cos(np.pi*lats[latxb:latxt+1]/180.).sum())
m,c   = np.polyfit(years[yx0:yx1+1],OHC_w,1)
line  = [m*x + c for x in years[yx0:yx1+1]]
pl.plot(years[yx0:yx1+1],OHC_w,'k.-',linewidth=1.75,alpha=0.7)
pl.plot(years[yx0:yx1+1],line ,'k--',linewidth=1.00,alpha=0.5)
pl.xlabel('Year')
pl.ylabel('Ocean heat content anomaly [10$^{18}$ J]')
pl.grid()
pl.xlim(years[yx0],years[yx1])
pl.title('trend = %s 10$^{18}$ J year$^{-1}$' % (round(m,3)))

"""
ux,uy = uo[xs,:,:,:].mean(axis=0).mean(axis=0),vo[xs,:,:,:].mean(axis=0).mean(axis=0)
ux,uy = (ux*dM[:,np.newaxis,np.newaxis]).sum(axis=0),(uy*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
cseq  = np.arange(-20000,20000+2500,2500)
cf    = pl.contourf(lons,lats,np.sqrt(ux**2 + uy**2),cseq,cmap=pl.cm.coolwarm,extend='both')
Q     = pl.quiver(lons,lats,ux,uy)
cbar  = pl.colorbar(cf)
cbar.set_label('Ocean mass transport [kg s$^{-1}$]')
pl.show()
"""


# Extract SeasMonths (left)
uol = uo[xs,levx0:levx1+1,latx0l:latx1l+1,lonxl]
Tl  =  T[xs,levx0:levx1+1,latx0l:latx1l+1,lonxl]
HTl = ((3985*((uol*Tl).mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*yres).sum(axis=1)
VTl = ((      (uol.mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*yres).sum(axis=1)/1028.
# Extract SeasMonths (right)
uor = uo[xs,levx0:levx1+1,latx0r:latx1r+1,lonxr]
Tr  =  T[xs,levx0:levx1+1,latx0r:latx1r+1,lonxr]
HTr = ((3985*((uor*Tr).mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*yres).sum(axis=1)
VTr = ((      (uor.mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*yres).sum(axis=1)/1028.
# Extract SeasMonths (bottom)
vob = vo[xs,levx0:levx1+1,latxb,lonx0b:lonx1b+1]
Tb  =  T[xs,levx0:levx1+1,latxb,lonx0b:lonx1b+1]
HTb = ((3985*((vob*Tb).mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*xres(70.5)).sum(axis=1)
VTb = ((      (vob.mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*xres(70.5)).sum(axis=1)/1028.
# Extract SeasMonths (top)
vot = vo[xs,levx0:levx1+1,latxt,lonx0t:lonx1t+1]
Tt  =  T[xs,levx0:levx1+1,latxt,lonx0t:lonx1t+1]
HTt = ((3985*((vot*Tt).mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*xres(81.5)).sum(axis=1)
VTt = ((      (vot.mean(axis=1))*dM[np.newaxis,levx0:levx1+1,np.newaxis]).sum(axis=1)*xres(81.5)).sum(axis=1)/1028.
# Extract SeasMonths (all)
HTa = HTl-HTr+HTb-HTt
VTa = VTl-VTr+VTb-VTt

# Cut to YearRange
yx0,yx1 = years.index(YearRange[0]),years.index(YearRange[1])
years   = years[yx0:yx1+1]
HTl,HTr,HTb,HTt = HTl[yx0:yx1+1],HTr[yx0:yx1+1],HTb[yx0:yx1+1],HTt[yx0:yx1+1]
VTl,VTr,VTb,VTt = VTl[yx0:yx1+1],VTr[yx0:yx1+1],VTb[yx0:yx1+1],VTt[yx0:yx1+1]
HTa,VTa = HTa[yx0:yx1+1],VTa[yx0:yx1+1]
# Trend (left)
slopeHTl,interceptHTl,r_valueHTl,p_valueHTl,std_errHTl = stats.linregress(years,HTl)
line = np.array([slopeHTl*x + interceptHTl for x in years])
htl  = HTl - line
slopeVTl,interceptVTl,r_valueVTl,p_valueVTl,std_errVTl = stats.linregress(years,VTl)
line = np.array([slopeVTl*x + interceptVTl for x in years])
vtl  = VTl - line
# Trend (right)
slopeHTr,interceptHTr,r_valueHTr,p_valueHTr,std_errHTr = stats.linregress(years,HTr)
line = np.array([slopeHTl*x + interceptHTl for x in years])
htl  = HTl - line
slopeVTr,interceptVTr,r_valueVTr,p_valueVTr,std_errVTr = stats.linregress(years,VTr)
line = np.array([slopeVTr*x + interceptVTr for x in years])
vtr  = VTr - line
# Trend (bottom)
slopeHTb,interceptHTb,r_valueHTb,p_valueHTb,std_errHTb = stats.linregress(years,HTb)
line = np.array([slopeHTb*x + interceptHTb for x in years])
htb  = HTb - line
slopeVTb,interceptVTb,r_valueVTb,p_valueVTb,std_errVTb = stats.linregress(years,VTb)
line = np.array([slopeVTb*x + interceptVTb for x in years])
vtb  = VTb - line
# Trend (top)
slopeHTt,interceptHTt,r_valueHTt,p_valueHTt,std_errHTt = stats.linregress(years,HTt)
line = np.array([slopeHTt*x + interceptHTt for x in years])
htt  = HTt - line
slopeVTt,interceptVTt,r_valueVTt,p_valueVTt,std_errVTt = stats.linregress(years,VTt)
line = np.array([slopeVTt*x + interceptVTt for x in years])
vtt  = VTt - line
# Trend (all)
slopeHTa,interceptHTa,r_valueHTa,p_valueHTa,std_errHTa = stats.linregress(years,HTa)
linea = np.array([slopeHTa*x + interceptHTa for x in years])
hta   = HTa - linea
slopeVTa,interceptVTa,r_valueVTa,p_valueVTa,std_errVTa = stats.linregress(years,VTa)
line = np.array([slopeVTa*x + interceptVTa for x in years])
vta  = VTa - line

#pl.figure(999)
#pl.plot(years,365*24*3600*np.cumsum(HTa)/1e18)
#pl.show()

"""
# Save timeseries
fname1 = '/mnt/climstorage/cian/BSO/ht.bso.full.%s-%s.%s.p' % (years[0],years[-1],Season)
fname2 = '/mnt/climstorage/cian/BSO/ht.bso.dtrend.%s-%s.%s.p' % (years[0],years[-1],Season)
if not os.path.isfile(fname1):
	toPick([years,HT],fname1)
if not os.path.isfile(fname2):
        toPick([years,ht],fname2)
"""

"""
# Fields for regression
df          = DataServer(Field=Field,Source='ERAInt')
F           = proj(np.array([df.getSeason(Year=year,Season=Season) for year in years]),df.lon,df.lat)
F,mf,pf     = detrend2d(years,F)
m,p         = np.zeros((proj.nx,proj.ny)),np.zeros((proj.nx,proj.ny))
mmoor,pmoor = np.zeros((proj.nx,proj.ny)),np.zeros((proj.nx,proj.ny))
for i in range(proj.nx):
	for j in range(proj.ny):
		#slope,intercept,r_value,p_value,std_err = stats.linregress(htmoor,F[:,i,j])
		#mmoor[i,j] = slope
		#pmoor[i,j] = p_value/2.
                slope,intercept,r_value,p_value,std_err = stats.linregress(ht,F[:,i,j])
                m[i,j] = slope
                p[i,j] = p_value/2
"""

# Plot
fig,ax1   = pl.subplots(num=7)
ax2       = ax1.twinx()
ax1.plot(years,HTl/1e12,'k.-',linewidth=1.75,alpha=0.7,label='Heat transport')
ax2.plot(years,VTl/1e06,'k.-',linewidth=0.6,alpha=0.5,label='Volume transport')
ax1.set_xlabel('Year')
ax1.set_ylabel('$HT_{left}$ [TW]')
ax2.set_ylabel('$VT_{left}$ [Sv]')
ax1.grid()
ax1.legend(loc=1,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
#pl.title('%s $HT_{left}$ (trend = %s TW year$^{-1}$; p = %s)' % (Season,round(slopeHTl/1e12,2),round(p_valueHTl/2.,6)))
ax1.set_xlim(years[0],years[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.left.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

fig,ax1   = pl.subplots(num=8)
ax2       = ax1.twinx()
ax1.plot(years,HTr/1e12,'k.-',linewidth=1.75,alpha=0.7,label='Heat transport')
ax2.plot(years,VTr/1e06,'k.-',linewidth=0.6,alpha=0.5,label='Volume transport')
ax1.set_xlabel('Year')
ax1.set_ylabel('$HT_{right}$ [TW]')
ax2.set_ylabel('$VT_{right}$ [Sv]')
ax1.grid()
ax1.legend(loc=1,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
#pl.title('%s $HT_{left}$ (trend = %s TW year$^{-1}$; p = %s)' % (Season,round(slopeHTl/1e12,2),round(p_valueHTl/2.,6)))
ax1.set_xlim(years[0],years[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.right.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

fig,ax1   = pl.subplots(num=9)
ax2       = ax1.twinx()
ax1.plot(years,HTb/1e12,'k.-',linewidth=1.75,alpha=0.7,label='Heat transport')
ax2.plot(years,VTb/1e06,'k.-',linewidth=0.6,alpha=0.5,label='Volume transport')
ax1.set_xlabel('Year')
ax1.set_ylabel('$HT_{bottom}$ [TW]')
ax2.set_ylabel('$VT_{bottom}$ [Sv]')
ax1.grid()
ax1.legend(loc=1,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
#pl.title('%s $HT_{left}$ (trend = %s TW year$^{-1}$; p = %s)' % (Season,round(slopeHTl/1e12,2),round(p_valueHTl/2.,6)))
ax1.set_xlim(years[0],years[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.bottom.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

fig,ax1   = pl.subplots(num=10)
ax2       = ax1.twinx()
ax1.plot(years,HTt/1e12,'k.-',linewidth=1.75,alpha=0.7,label='Heat transport')
ax2.plot(years,VTt/1e06,'k.-',linewidth=0.6,alpha=0.5,label='Volume transport')
ax1.set_xlabel('Year')
ax1.set_ylabel('$HT_{top}$ [TW]')
ax2.set_ylabel('$VT_{top}$ [Sv]')
ax1.grid()
ax1.legend(loc=1,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
#pl.title('%s $HT_{left}$ (trend = %s TW year$^{-1}$; p = %s)' % (Season,round(slopeHTl/1e12,2),round(p_valueHTl/2.,6)))
ax1.set_xlim(years[0],years[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.top.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

# Plot
fig,ax1   = pl.subplots(num=11)
ax2       = ax1.twinx()
ax1.plot(years,HTa/1e12,'k.-',linewidth=1.75,alpha=0.7,label='Heat transport')
ax1.plot(years,linea/1e12,'k--',linewidth=1,alpha=0.5)
ax2.plot(years,VTa/1e06,'k.-',linewidth=0.6,alpha=0.5,label='Volume transport')
print VTa.mean()/1e06
print HTa.mean()/1e12
ax1.set_xlabel('Year')
ax1.set_ylabel('$HT_{all}$ [TW]')
ax2.set_ylabel('$VT_{all}$ [Sv]')
ax1.grid()
ax1.legend(loc=1,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.title('%s $HT_{all}$ (trend = %s TW year$^{-1}$; p = %s)' % (Season,round(slopeHTa/1e12,2),round(p_valueHTa/2.,6)))
ax1.set_xlim(years[0],years[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht.all.%s-%s.%s.pdf' % (YearRange[0],YearRange[1],Season),format='pdf')

"""
pl.figure(2)
cseq  = 13#np.arange(-30,30+5,5)
sx,sy = np.ma.masked_where(pmoor>0.05,proj.x),np.ma.masked_where(pmoor>0.05,proj.y)
bx,by = proj.m([20 for i in range(len(lats[latx0:latx1+1]))],lats[latx0:latx1+1])
cf    = pl.contourf(proj.x,proj.y,1e12*mmoor,cseq,cmap=pl.cm.coolwarm,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (df.units))
pl.plot(sx,sy,'k.',alpha=0.5)
pl.plot(bx,by,'g-',linewidth=1.75,alpha=0.6)
proj.m.drawparallels([70,80],latmax=90)
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
pl.title('BSO mooring data')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.moor.regress.%s-%s.%s.pdf' % (Field,YearRange[0],YearRange[1],Season),format='pdf')
"""

"""
pl.figure(3)
cseq  = 13#np.arange(-30,30+5,5)
sx,sy = np.ma.masked_where(p>0.05,proj.x),np.ma.masked_where(p>0.05,proj.y)
bx,by = proj.m([20 for i in range(len(lats[latx0:latx1+1]))],lats[latx0:latx1+1])
cf    = pl.contourf(proj.x,proj.y,1e12*m,cseq,cmap=pl.cm.coolwarm,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Reg. coeff. [%s TW$^{-1}$]' % (df.units))
pl.plot(sx,sy,'k.',alpha=0.5)
pl.plot(bx,by,'g-',linewidth=1.75,alpha=0.6)
proj.m.drawparallels([70,80],latmax=90)
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
pl.title('ORAS4')
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/ht_%s.reanal.regress.%s-%s.%s.pdf' % (Field,YearRange[0],YearRange[1],Season),format='pdf')
"""

pl.show()

