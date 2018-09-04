from netCDF4 import Dataset
from ReanalysisMonthlyMeanDataServer import DataServer
from scipy import stats
from mpl_toolkits.basemap import Basemap

import numpy as np
import matplotlib.pyplot as pl

SMs = {}
SMs['D']         = [12]
SMs['J']         = [1]
SMs['DJF']       = [12,1,2]
SMs['JJA']       = [6,7,8]
SMs['SON']       = [9,10,11]
SMs['MAM']       = [3,4,5]
SMs['Annual']    = range(1,12+1,1)
SMs['Annual_wc'] = range(7,12+1,1) + range(1,6+1,1)

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
        return xr

Season    = 'Annual'
YearRange = (1980,2015)

# Data and axes
d    = Dataset('/mnt/climstorage/cian/EN4_T_S_BS.nc','r')
time = d.variables['time'][:]
T    = d.variables['temperature'][:]
Z    = d.variables['depth'][:]
lat  = d.variables['lat'][:]
lon  = d.variables['lon'][:]
# Basemap
bmap = Basemap(projection='cyl',llcrnrlat=68.5,urcrnrlat=83.5,llcrnrlon=0.5,urcrnrlon=100.5,resolution='i')
j,k  = bmap(lon,lat)
# Time indexes
years      = range(1950,2015+1,1)
SeasMonths = SMs[Season]
months     = np.tile(np.arange(1,12+1,1),(len(time)/12.))
xs         = [range(i,i+len(SeasMonths)) for i in range(len(months)) if (list(months[i:i+len(SeasMonths)]) == SeasMonths)]
if SeasMonths[0]>SeasMonths[-1]: years = years[1:]
yx0,yx1    = years.index(YearRange[0]),years.index(YearRange[1])
# Mass of girdbox on each level
dD_ = []
dD  = np.diff(Z)/2.
dD_.append(Z[0]+dD[0])
for i in range(len(dD)-1):
        dD_.append(dD[i]+dD[i+1])
dD_.append(dD[-1])#+(1000-levs[-1]))
dD_ = np.array(dD_)
dM  = 1028.*dD_

# OHC anomalies
OHC            = (3985*T[xs,:,:,:].mean(axis=1)*dM[np.newaxis,:,np.newaxis,np.newaxis]).sum(axis=1)
dOHC,mOHC,pOHC = detrend2d(years[yx0:yx1+1],OHC[yx0:yx1+1])
mOHC           = np.ma.masked_array(mOHC,mask=OHC[0].mask)
# Time series
OHC_w = 1.15e12*(OHC.mean(axis=-1)*np.cos(np.pi*lat/180)[np.newaxis,:]).sum(axis=1)/np.cos(np.pi*lat/180).sum()
#OHC_w = (( (OHC)*(xres(lat)[np.newaxis,:,np.newaxis])).sum(axis=-1)*xres(0)).sum(axis=-1)
OHC_w = OHC_w - OHC_w.mean()

# Plot
fig       = pl.figure(1)
ax        = fig.add_subplot(111)
cseq      = np.arange(-10,10+1.25,1.25)
cf        = ax.contourf(j,k,10*mOHC/1e8,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar      = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Ocean heat content trend [10$^{8}$ J m$^{-2}$ decade$^{-1}$]')
bmap.drawcoastlines(linewidth=0.5,color='0.5')
bmap.fillcontinents(color='0.85')
bmap.drawmeridians([0,20,40,60,80,100])
bmap.drawparallels([70,75,80])
ax.set_aspect('auto')

pl.figure(2)
m,c  = np.polyfit(years[yx0:yx1+1],OHC_w[yx0:yx1+1],1)
line = np.array([m*x + c for x in years[yx0:yx1+1]])
pl.plot(years[yx0:yx1+1],OHC_w[yx0:yx1+1]/1e18,'k',linewidth=2,alpha=0.65)
pl.plot(years[yx0:yx1+1],line/1e18,'k--',linewidth=1,alpha=0.5)
pl.grid()
pl.ylabel('Ocean heat content anomaly [10$^{18}$ J]')
pl.xlabel('Year')
pl.title('trend = %s 10$^{18}$ year$^{-1}$' % (round(m/1e18,3)))
pl.xlim(years[yx0],years[yx1])

pl.show()
