from netCDF4 import Dataset as Dataset
from UnPickle import *

import sys
import numpy as np
import matplotlib.pyplot as pl

# FluxFile
vq     = unpick('newfluxfiles/70/ERAInt/ERAInt.moist.1980-2015.30-1000hPa.70N.DJF.p')[0]
inds   = np.arange(-4,4+1,1)
vq_    = []
for i in range(len(vq[0])):
	indx = i + inds
	for il in range(len(indx)):
		if indx[il]<0:    indx[il] = indx[il] + 360
		if indx[il]>=360: indx[il] = indx[il] - 360
	vq_.append(vq[:,indx].mean(axis=1))
vq = np.array(vq_)
vq = np.rollaxis(vq,1,0).reshape(-1)
vq     = [vq[i] for i in range(len(vq)) if vq[i] > 0]
# DataFile attributes
Dir    = '../Tracks/back1_crossings'
maskc  = 9.96921e36/2.
#Years1 = [1982,1989,1998,1999,2003,2004]
#Years2 = [1985,1990,1991,1995,2000,2002]
Years1 = range(1981,1998+1,1)
Years2 = range(1999,2016+1,1)
Months = [1,2,3,11,12]
flux1,lonc1,step1,pw1,latst1,lonst1 = [],[],[],[],[],[]
flux2,lonc2,step2,pw2,latst2,lonst2 = [],[],[],[],[],[]
# Get data
for Year in Years1:
	for Month in Months:
		fname = '%s/cross_%s_%02d.nc' % (Dir,Year,Month)
		d = Dataset(fname,'r')
	flux1  =  flux1 + [i for i in     d.variables['Flux'][:].reshape(-1) if i <= maskc]
	lonc1  =  lonc1 + [i for i in d.variables['LonCross'][:].reshape(-1) if i <= maskc]
	step1  =  step1 + [i for i in     d.variables['Step'][:].reshape(-1) if i <= maskc]
	pw1    =    pw1 + [i for i in       d.variables['PW'][:].reshape(-1) if i <= maskc]
	latst1 = latst1 + [i for i in d.variables['LatStart'][:].reshape(-1) if i <= maskc]
	#lonst1 = lonst1 + [i for i in      d.variables['lon'][:].reshape(-1) if i <= maskc]
for Year in Years2:
        for Month in Months:
                fname = '%s/cross_%s_%02d.nc' % (Dir,Year,Month)
                d = Dataset(fname,'r')
        flux2  =  flux2 + [i for i in     d.variables['Flux'][:].reshape(-1) if i <= maskc]
        lonc2  =  lonc2 + [i for i in d.variables['LonCross'][:].reshape(-1) if i <= maskc]
        step2  =  step2 + [i for i in     d.variables['Step'][:].reshape(-1) if i <= maskc]
        pw2    =    pw2 + [i for i in       d.variables['PW'][:].reshape(-1) if i <= maskc]
        latst2 = latst2 + [i for i in d.variables['LatStart'][:].reshape(-1) if i <= maskc]
	#lonst2 = lonst2 + [i for i in np.tile(d.variables['lon'][np.newaxis,:],(nt,1)).reshape(-1) if i <= maskc]
maxstep = 41
lonc1_  = [ lonc1[i] for i in range(len(flux1)) if (latst1>80) and (step1[i]<=maxstep)]
lonc2_  = [ lonc2[i] for i in range(len(flux2)) if (latst2>80) and (step2[i]<=maxstep)]
print len(lonc1_),len(lonc2_)
h1,e1   = np.histogram(lonc1_,50,normed=False)
h2,e2   = np.histogram(lonc2_,50,normed=False)
pl.plot(e1[0:-1],h1,'b')
pl.plot(e2[0:-1],h2,'r')
pl.ylabel('Frequency')
pl.xlabel('Longitude')
pl.grid()
pl.show()

# Masking criterion
LonRange  = [int(sys.argv[1]),int(sys.argv[2])]
maxstep   = 4*10 + 1
blat      = 80
pwbs,bins = zip(range(0,7+1,1),range(1,8+1,1)),np.arange(0.5,7.5+1,1)
# LonRange hacking
#lonc_n = np.array(lonc)
#lonc_n[np.where(lonc_n<0)] = lonc_n[np.where(lonc_n<0)] + 360 
# Apply masks
FLUX1m,STEP1m,LONC1m = [],[],[]
FLUX2m,STEP2m,LONC2m = [],[],[]
FLUX1s,STEP1s,LONC1s = [],[],[]
FLUX2s,STEP2s,LONC2s = [],[],[]
N1,N2                = [],[]
for pwb0,pwb1 in pwbs:
	indices1   = [ i for i in range(len(flux1)) if (latst1[i]>blat) and (flux1[i]>0) and (LonRange[0]<lonc1[i]<LonRange[1]) and (pwb0<=pw1[i]<pwb1) and (step1[i]<maxstep)]
	N1.append(len(indices1)/1000.)
	f1,st1,lc1 = np.array([ flux1[i] for i in indices1]),np.array([ step1[i] for i in indices1]),np.array([ lonc1[i] for i in indices1])
	FLUX1m.append(f1.mean())
	STEP1m.append(st1.mean())
	LONC1m.append(lc1.mean())
        FLUX1s.append(f1.std())
        STEP1s.append(st1.std())
        LONC1s.append(lc1.std())
	indices2   = [ i for i in range(len(flux2)) if (latst2[i]>blat) and (flux2[i]>0) and (LonRange[0]<lonc2[i]<LonRange[1]) and (pwb0<=pw2[i]<pwb1) and (step2[i]<maxstep)]
	N2.append(len(indices2)/1000.)
	f2,st2,lc2 = np.array([ flux2[i] for i in indices2]),np.array([ step2[i] for i in indices2]),np.array([ lonc2[i] for i in indices2])
        FLUX2m.append(f2.mean())
        STEP2m.append(st2.mean())
        LONC2m.append(lc2.mean())
        FLUX2s.append(f2.std())
        STEP2s.append(st2.std())
        LONC2s.append(lc2.std())
FLUX1m,STEP1m,LONC1m = np.array(FLUX1m),np.array(STEP1m),np.array(LONC1m)
FLUX2m,STEP2m,LONC2m = np.array(FLUX2m),np.array(STEP2m),np.array(LONC2m)
FLUX1s,STEP1s,LONC1s = np.array(FLUX1s),np.array(STEP1s),np.array(LONC1s)
FLUX2s,STEP2s,LONC2s = np.array(FLUX2s),np.array(STEP2s),np.array(LONC2s)
# Plotting
fig,ax1 = pl.subplots(num=1)
ax2     = ax1.twinx()
ax1.plot(bins,FLUX1m,'b',linewidth=1.75,alpha=0.7)
ax1.plot(bins,FLUX2m,'r',linewidth=1.75,alpha=0.7)
ax1.grid()
ax2.plot(bins,N1,'b--',linewidth=1,alpha=0.5)
ax2.plot(bins,N2,'r--',linewidth=1,alpha=0.5)
pl.title('%sE to %sE' % (LonRange[0],LonRange[1]))
ax1.set_xlabel('Preciptable water [kg m$^{-2}$]')
ax1.set_ylabel('Characteristic flux [Tg day$^{-1}$ deg$^{-1}$]')
ax2.set_ylabel('Frequency [10$^{3}$]')
#pl.ylim(0,300)
pl.xlim(0,9)
pl.show()
