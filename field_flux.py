from scipy import interpolate,stats
from scipy import ndimage as nd
from cmipDataServer import DataServer as cmipDataServer
from UnPickle import *

import glob
import numpy as np
import matplotlib.pyplot as pl

def interp(vq):
        n0,n1 = vq.shape
        # Add wrap value to vq
        ends = vq[:,-1]
        vq   = np.append(vq,ends[:,np.newaxis],1)
        # Interpolate vq to 1 degree resolution
        yold = range(n0)
        ynew = range(n0)
        xold = np.array(list(np.arange(0,360,360./n1))+[360])
        xnew = np.arange(0,360,1.)
        f    = interpolate.interp2d(xold,yold,vq,kind='linear',bounds_error=True)
        vq   = f(xnew,ynew)
        return vq

def fill(self,data,invalid=None):
        if invalid is None: invalid = np.isnan(data)
        ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
        return data[tuple(ind)]

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

def fluxToYears(vq,dates):
        vqt = {}
        dt  = {}
        for i in range(len(vq)):
                year,month,day,hour = dates[i]
                if month > 8: year = year + 1
                try:
                        vqt[year].append(vq[i])
                        dt[year].append(dates[i])
                except:
                        vqt[year] = []
                        dt[year]  = []
                        vqt[year].append(vq[i])
                        dt[year].append(dates[i])
        return vqt.values(),dt.values(),vqt.keys()

Re       = 6371
xres     = 2*np.pi*np.cos(70*np.pi/180)*Re/360.
sf       = 1e09/(24.*60*60*xres*1000)
lats1    = np.linspace(50,80,20)
#levs    = np.linspace(700,1000,20)
lats2    = np.linspace(-15,15,20)
levs     = np.linspace(150,450,20)
coslats1 = np.abs(np.cos(lats1*np.pi/180.))
coslats2 = np.abs(np.cos(lats2*np.pi/180.))
years1   = range(1981,2005+1,1)
years2   = range(2075,2095+1,1)
thresh   = 0
Modelsno = ['GFDL-CM3','inmcm4','FGOALS-g2']
Models   = [i[9:] for i in glob.glob('../rcp85/*') if i[9:] not in Modelsno]
print Models
dT,dF    = [],[]
for Source in Models:
	# DataServer
	dsRCP = cmipDataServer(Field='ta',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='mon')
	dsHST = cmipDataServer(Field='ta',LevType='plev',Source=Source,ExpType='historical',DataFreq='mon')

	T_HST = np.array([dsHST.getSeason(Year=year,Season='DJF') for year in years1])#.mean(axis=0)
	if (type(T_HST)==np.ma.core.MaskedArray) and (T_HST.mask.shape == T_HST.data.shape):
	        T_HST = self.fill(T_HST.data,invalid=T_HST.mask)
	T_HST = T_HST.mean(axis=0)
	T_HST = interpolateND(T_HST,dsHST.lev,levs,axis=0)

	T_HST1 = interpolateND(T_HST,dsHST.lat,lats1,axis=1).mean(axis=0).mean(axis=1)
	T_HST1 = (T_HST1*coslats1).sum()/coslats1.sum()
        T_HST2 = interpolateND(T_HST,dsHST.lat,lats2,axis=1).mean(axis=0).mean(axis=1)
        T_HST2 = (T_HST2*coslats2).sum()/coslats2.sum()

        T_RCP = np.array([dsRCP.getSeason(Year=year,Season='DJF') for year in years2])#.mean(axis=0)
        if (type(T_RCP)==np.ma.core.MaskedArray) and (T_RCP.mask.shape == T_RCP.data.shape):
                T_RCP = self.fill(T_RCP.data,invalid=T_RCP.mask)
	T_RCP = T_RCP.mean(axis=0)
        T_RCP = interpolateND(T_RCP,dsRCP.lev,levs,axis=0)

        T_RCP1 = interpolateND(T_RCP,dsRCP.lat,lats1,axis=1).mean(axis=0).mean(axis=1)
	T_RCP1 = (T_RCP1*coslats1).sum()/coslats1.sum()
        T_RCP2 = interpolateND(T_RCP,dsRCP.lat,lats2,axis=1).mean(axis=0).mean(axis=1)
        T_RCP2 = (T_RCP2*coslats2).sum()/coslats2.sum()

        dT.append((T_RCP2 - T_RCP1) - (T_HST2 - T_HST1))

        # Open historical flux files
        hdir               = '/qnap/cian/cmip/scripts'
        vqht0,datesh       = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))
        vht0,datesh        = unpick('%s/newfluxfiles/%s/%s/%s.mass.1981-2005.600-1000hPa.%sN.DJF.p'  % (hdir,70,Source,Source,70))
        qht0               = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))[0]/(40000/9.8)
        vqht0,vht0,qht0    = interp(vqht0),interp(vht0),interp(qht0)
        vqht,vht,qht       = sf*np.ma.masked_where(vqht0<thresh,vqht0),sf*np.ma.masked_where(vqht0<thresh,vht0),np.ma.masked_where(vqht0<thresh,qht0)
        # Open rcp85 flux files
        rdir               = '/mnt/climstorage/cian/scripts'
        vqrt0,datesr       = unpick('%s/newfluxfiles/%s/%s/%s.moist.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))
        vrt0,datesr        = unpick('%s/newfluxfiles/%s/%s/%s.mass.2075-2100.600-1000hPa.%sN.DJF.p'  % (rdir,70,Source,Source,70))
        qrt0               = unpick('%s/newfluxfiles/%s/%s/%s.vapor.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))[0]/(40000/9.8)
        vqrt0,vrt0,qrt0    = interp(vqrt0),interp(vrt0),interp(qrt0)
        vqrt,vrt,qrt       = sf*np.ma.masked_where(vqrt0<thresh,vqrt0),sf*np.ma.masked_where(vqrt0<thresh,vrt0),np.ma.masked_where(vqrt0<thresh,qrt0)
        # Changes in means
        delvq = vqrt.mean() - vqht.mean()
	dF.append(delvq)

slope, intercept, r_value, p_value, std_err = stats.linregress(dF,dT)
xs = np.linspace(6,20,10)
line = [slope*x + intercept for x in xs]
pl.plot(dF,dT,'k+',markersize=9,mew=2)
pl.plot(xs,line,'k--',linewidth=1.75,alpha=0.5,label=r'$y = %sx + %s$' % (round(slope,3),round(intercept,3)))
for label, x, y in zip(range(1,len(Models)+1), dF, dT):
	print label,Models[label-1]
	pl.annotate(   label,xy=(x, y),xytext=(5.5,2.5),color='k',
                       textcoords='offset points', ha='center', va='bottom',alpha=1,fontsize=8)
pl.ylabel('%d-%dN %d-%dhPa temp change [K]' % (lats1[0],lats1[-1],levs[0],levs[-1]))
pl.xlabel('70N vq+ change [kg m$^{-1}$ s$^{-1}$]')
pl.title('(RCP: 2075 - 2095) - (HIST: 1981 - 2005) [DJF] (p=%s)' % (round(p_value/2.,4)))
pl.grid()
pl.xlim(6,20)
pl.legend(loc=2)
#pl.savefig('figs/%d-%dN.%d-%dhPa.pdf' % (lats[0],lats[-1],levs[0],levs[-1]),format='pdf')
pl.show()
