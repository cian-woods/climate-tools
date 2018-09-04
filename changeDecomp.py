from UnPickle import unpick
from scipy import interpolate

import glob,sys
import numpy as np
import matplotlib.pyplot as pl

def moving_sum(a,n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:]

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

# Parameters
Re       = 6371
xres     = 2*np.pi*np.cos(70*np.pi/180)*Re/360.
sf       = 1e09/(24.*60*60*xres*1000)
lons     = np.arange(0,360,1.)
thresh   = int(sys.argv[1])

T1,T2,T3 = [],[],[]
T4,T5,T6 = [],[],[]
Models   = [i[9:] for i in glob.glob('../rcp85/*')]
for Source in Models:
        print Source
        # Open historical flux files
        hdir               = '/qnap/cian/cmip/scripts'
        vqht0,datesh       = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))
        vht0,datesh        = unpick('%s/newfluxfiles/%s/%s/%s.mass.1981-2005.600-1000hPa.%sN.DJF.p'  % (hdir,70,Source,Source,70))
        qht0               = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))[0]/(40000/9.8)
	vqht0,vht0,qht0    = interp(vqht0),interp(vht0),interp(qht0)
	vqht,vht,qht       = sf*np.ma.masked_where(vqht0<thresh,vqht0),sf*np.ma.masked_where(vqht0<thresh,vht0),np.ma.masked_where(vqht0<thresh,qht0)
	#vqhbar,vhbar,qhbar = vqht.mean(axis=0),vht.mean(axis=0),qht.mean(axis=0)
	#vqhprm,vhprm,qhprm = vqht-vqhbar[np.newaxis,:],vht-vhbar[np.newaxis,:],qht-qhbar[np.newaxis,:]
        # Open rcp85 flux files
        rdir               = '/mnt/climstorage/cian/scripts'
        vqrt0,datesr       = unpick('%s/newfluxfiles/%s/%s/%s.moist.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))
        vrt0,datesr        = unpick('%s/newfluxfiles/%s/%s/%s.mass.2075-2100.600-1000hPa.%sN.DJF.p'  % (rdir,70,Source,Source,70))
        qrt0               = unpick('%s/newfluxfiles/%s/%s/%s.vapor.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))[0]/(40000/9.8)
	vqrt0,vrt0,qrt0    = interp(vqrt0),interp(vrt0),interp(qrt0)
	vqrt,vrt,qrt       = sf*np.ma.masked_where(vqrt0<thresh,vqrt0),sf*np.ma.masked_where(vqrt0<thresh,vrt0),np.ma.masked_where(vqrt0<thresh,qrt0)
	#vqrbar,vrbar,qrbar = vqrt.mean(axis=0),vrt.mean(axis=0),qrt.mean(axis=0)
	#vqrprm,vrprm,qrprm = vqrt-vqrbar[np.newaxis,:],vrt-vrbar[np.newaxis,:],qrt-qrbar[np.newaxis,:]

        # Changes in means
        #delvq,delv,delq = vqrbar-vqhbar,vrbar-vhbar,qrbar-qhbar
	# Product of primes
	#vhqhbar = (vhprm*qhprm).mean(axis=0)
	#vrqrbar = (vrprm*qrprm).mean(axis=0)
	# Terms
	#term1 = delv*qhbar
	#term2 = delq*vhbar
	#term3 = delv*delq
	#term4 = (vrqrbar - vhqhbar)
	#term5 = term1+term2+term3+term4

	vqhbar = vqht.sum(axis=0)/len(vqht)#np.ma.masked_where(vqt1<thresh,vqt1).sum(axis=0)/len(vqt1)
	vhbar  = vht.sum(axis=0)/len(vqht)#np.ma.masked_where(vqt1<thresh, vt1).sum(axis=0)/len(vqt1)
	qhbar  = qht.mean(axis=0)#np.ma.masked_where(vqt1<thresh, qt1).mean(axis=0)

	vqrbar = vqrt.sum(axis=0)/len(vqrt)#np.ma.masked_where(vqt2<thresh,vqt2).sum(axis=0)/len(vqt2)
	vrbar  =  vrt.sum(axis=0)/len(vqrt)#np.ma.masked_where(vqt2<thresh, vt2).sum(axis=0)/len(vqt2)
	qrbar  =  qrt.mean(axis=0)#np.ma.masked_where(vqt2<thresh, qt2).mean(axis=0)

	delvq,delv,delq = vqrbar-vqhbar,vrbar-vhbar,qrbar-qhbar

	term1 = qhbar*delv
	term2 = vhbar*delq
	term3 = delv*delq
	term4 = delvq
	term5 = delvq-(term1+term2+term3)

	T1.append(term1)
	T2.append(term2)
	T3.append(term3)
	T4.append(term4)
	T5.append(term5)
	#T6.append(delvq)

T1m,T2m,T3m = np.array(T1).mean(axis=0),np.array(T2).mean(axis=0),np.array(T3).mean(axis=0)
T4m,T5m     = np.array(T4).mean(axis=0),np.array(T5).mean(axis=0)
T4s = np.array(T4).std(axis=0)

T1m = moving_sum(T1m,11)
T2m = moving_sum(T2m,11)
T3m = moving_sum(T3m,11)
T4m = moving_sum(T4m,11)
T5m = moving_sum(T5m,11)
T4s = moving_sum(T4s,11)
lons = lons[5:360-5]

pl.figure(1)
pl.plot(lons,T1m,'r',linewidth=1.5,alpha=0.8,label='dV.Q')
pl.plot(lons,T2m,'b',linewidth=1.5,alpha=0.8,label='dQ.V')
pl.plot(lons,T3m,'g',linewidth=1.25,alpha=0.6,label='dV.dQ')
pl.plot(lons,T4m,'DarkMagenta',linewidth=2,alpha=1,label='Total')
pl.fill_between(lons,T4m-T4s,T4m+T4s,color='DarkMagenta',alpha=0.1)
pl.plot(lons,T5m,'k--',linewidth=1.25,alpha=0.6,label='residual')
pl.xlabel('Longitude')
pl.ylabel('Change in moisture flux [kg m s$^{-1}$]')
pl.title('(RCP: 2075-2100) - (HIST: 1981-2005) [DJF]')
pl.xlim(0,359)
if thresh == 0:
        pl.ylim(-30,160)
else:
        pl.ylim(-150,150)
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.savefig('figs/decomp.CMIP5.%s.pdf' % (thresh),format='pdf')

"""
# ERAInt
# Open ERAInt flux file
Dir            = '/mnt/climstorage/cian/scripts'
vqt0,datesh    = unpick('%s/newfluxfiles/%s/%s/%s.moist.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))
vt0,datesh     = unpick('%s/newfluxfiles/%s/%s/%s.mass.1980-2016.600-1000hPa.%sN.DJF.p'  % (Dir,70,'ERAInt','ERAInt',70))
qt0            = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))[0]/(40000/9.8)
date1,date2    = (1991,12,31,18),(2005,1,1,0)
ii1,ii2        = datesh.index(date1),datesh.index(date2)
vqt0,vht0,qht0 = interp(vqt0),interp(vt0),interp(qt0)
vqt,vt,qt      = sf*np.ma.masked_where(vqt0<thresh,vqt0),sf*np.ma.masked_where(vqt0<thresh,vt0),np.ma.masked_where(vqt0<thresh,qt0)

vqt1,vt1,qt1 = vqt[0:ii1,:],vt[0:ii1,:],qt[0:ii1,:]
vqt2,vt2,qt2 = vqt[ii2:,:],vt[ii2:,:],qt[ii2:,:]

vqbar1,vbar1,qbar1 = vqt1.mean(axis=0),vt1.mean(axis=0),qt1.mean(axis=0)
vqprm1,vprm1,qprm1 = vqt1-vqbar1[np.newaxis,:],vt1-vbar1[np.newaxis,:],qt1-qbar1[np.newaxis,:]

vqbar2,vbar2,qbar2 = vqt2.mean(axis=0),vt2.mean(axis=0),qt2.mean(axis=0)
vqprm2,vprm2,qprm2 = vqt2-vqbar2[np.newaxis,:],vt2-vbar2[np.newaxis,:],qt2-qbar2[np.newaxis,:]

#pl.plot(lons,vqbar1,'b')
pl.plot(lons,vqbar2-vqbar1,'r')
pl.show()

# Changes in means
delvq,delv,delq = vqbar2-vqbar1,vbar2-vbar1,qbar2-qbar1
# Product of primes
v1q1bar = (vprm1*qprm1).mean(axis=0)
v2q2bar = (vprm2*qprm2).mean(axis=0)
# Terms
term1 = delv*qbar1
term2 = delq*vbar1
term3 = delv*delq
term4 = (v2q2bar - v1q1bar)
term5 = term1+term2+term3+term4
print delvq[np.where(delvq>0)].sum(),delvq[np.where(delvq<0)].sum()

pl.plot(lons,term1,'r',linewidth=1.5,alpha=0.8,label='dV.Q')
pl.plot(lons,term2,'b',linewidth=1.5,alpha=0.8,label='dQ.V')
pl.plot(lons,term3,'g',linewidth=1.25,alpha=0.6,label='dV.dQ')
pl.plot(lons,term4,'y',linewidth=1.25,alpha=0.6,label='d(VpQp)')
pl.plot(lons,term5,'DarkMagenta',linewidth=2,alpha=1,label='Total')
pl.xlabel('Longitude')
pl.ylabel('Change in moisture flux')
pl.title('(ERAInt: 2005-2016) - (ERAInt: 1980-1991) [DJF]')
pl.xlim(0,359)
pl.ylim(-15,15)
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.savefig('figs/decomp.ERAInt.%s.pdf' % (thresh),format='pdf')
pl.show()
"""

"""
# ERAInt
# Open ERAInt flux file
Dir            = '/mnt/climstorage/cian/scripts'
vqt0,datesh    = unpick('%s/newfluxfiles/%s/%s/%s.moist.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))
vt0,datesh     = unpick('%s/newfluxfiles/%s/%s/%s.mass.1980-2016.600-1000hPa.%sN.DJF.p'  % (Dir,70,'ERAInt','ERAInt',70))
qt0            = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))[0]/(40000/9.8)
date1,date2    = (1999,12,31,18),(1999,12,31,18)
#date1,date2   = (1991,12,31,18),(2005,1,1,0)
ii1,ii2        = datesh.index(date1),datesh.index(date2)
vqt0,vht0,qht0 = interp(vqt0),interp(vt0),interp(qt0)
vqt,vt,qt      = sf*np.ma.masked_where(vqt0<thresh,vqt0),sf*np.ma.masked_where(vqt0<thresh,vt0),np.ma.masked_where(vqt0<thresh,qt0)

vqt1,vt1,qt1 = vqt[0:ii1,:],vt[0:ii1,:],qt[0:ii1,:]
vqt2,vt2,qt2 = vqt[ii2:,:],vt[ii2:,:],qt[ii2:,:]

vqbar1 = vqt1.sum(axis=0)/len(vqt1)#np.ma.masked_where(vqt1<thresh,vqt1).sum(axis=0)/len(vqt1)
vbar1  = vt1.sum(axis=0)/len(vqt1)#np.ma.masked_where(vqt1<thresh, vt1).sum(axis=0)/len(vqt1)
qbar1  = qt1.mean(axis=0)#np.ma.masked_where(vqt1<thresh, qt1).mean(axis=0)

vqbar2 = vqt2.sum(axis=0)/len(vqt2)#np.ma.masked_where(vqt2<thresh,vqt2).sum(axis=0)/len(vqt2)
vbar2  =  vt2.sum(axis=0)/len(vqt2)#np.ma.masked_where(vqt2<thresh, vt2).sum(axis=0)/len(vqt2)
qbar2  =  qt2.mean(axis=0)#np.ma.masked_where(vqt2<thresh, qt2).mean(axis=0)

delvq,delv,delq = vqbar2-vqbar1,vbar2-vbar1,qbar2-qbar1

term1 = qbar1*delv
term2 = vbar1*delq
term3 = delv*delq
term4 = delvq
term5 = delvq-(term1+term2+term3)

term1 = moving_sum(term1,11)
term2 = moving_sum(term2,11)
term3 = moving_sum(term3,11)
term4 = moving_sum(term4,11)
term5 = moving_sum(term5,11)
#lons = lons[5:360-5]

pl.figure(2)
pl.plot(lons,term1,'r',linewidth=1.5,alpha=0.8,label='dV.Q')
pl.plot(lons,term2,'b',linewidth=1.5,alpha=0.8,label='dQ.V')
pl.plot(lons,term3,'g',linewidth=1.25,alpha=0.6,label='dV.dQ')
pl.plot(lons,term4,'DarkMagenta',linewidth=2,alpha=1,label='Total')
pl.plot(lons,term5,'k--',linewidth=1.25,alpha=0.6,label='residual')
pl.xlabel('Longitude')
pl.ylabel('Change in moisture flux [kg m s$^{-1}$]')
pl.title('(ERAInt: 1999-2016) - (ERAInt: 1980-1998) [DJF]')
pl.xlim(0,359)
if thresh == 0:
	pl.ylim(-30,160)
else:
	pl.ylim(-150,150)
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.savefig('figs/decomp.ERAInt.%s.pdf' % (thresh),format='pdf')
pl.show()
"""

"""
# ERAInt
# Open ERAInt flux file
Dir            = '/mnt/climstorage/cian/scripts'
vqt0,datesh    = unpick('%s/newfluxfiles/%s/%s/%s.moist.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))
vt0,datesh     = unpick('%s/newfluxfiles/%s/%s/%s.mass.1980-2016.600-1000hPa.%sN.DJF.p'  % (Dir,70,'ERAInt','ERAInt',70))
qt0            = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1980-2016.600-1000hPa.%sN.DJF.p' % (Dir,70,'ERAInt','ERAInt',70))[0]/(40000/9.8)

# Date to yearly format
vqt0,datesht,years = fluxToYears(vqt0,datesh)
vt0,datesht,years = fluxToYears(vt0,datesh)
qt0,datesht,years = fluxToYears(qt0,datesh)
# Means for each year
vqt = np.array([np.ma.masked_where(np.array(vqt0[i])<thresh,np.array(vqt0[i])).mean(axis=0) for i in range(len(years))])
vt  = np.array([np.ma.masked_where(np.array(vqt0[i])<thresh,np.array( vt0[i])).mean(axis=0) for i in range(len(years))])
qt  = np.array([np.ma.masked_where(np.array(vqt0[i])<thresh,np.array( qt0[i])).mean(axis=0) for i in range(len(years))])
# Climatologies
vqclim = vqt.mean(axis=0)
vclim  = vt.mean(axis=0)
qclim  = qt.mean(axis=0)
# Prime series
vqprm = vqt - vqclim[np.newaxis,:]
vprm =   vt - vclim[np.newaxis,:]
qprm =   qt - qclim[np.newaxis,:]
# Product of prime series
vprmqprm = vprm*qprm
# Trend of various terms
vqtrend   = np.polyfit(years,vqt,1)[0]
vprmt     = np.polyfit(years,vprm,1)[0]
qprmt     = np.polyfit(years,qprm,1)[0]
vprmqprmt = np.polyfit(years,vprmqprm,1)[0]
#vprmt     = getTrend(years,vprm)
#qprmt     = getTrend(years,qprm)
#vprmqprmt = getTrend(years,vprmqprm)
# Terms
term1 = vprmt*qclim
term2 = qprmt*vclim
term3 = vprmqprmt
term4 = vqtrend
term5 = vqtrend - (term1+term2+term3)

term1 = moving_sum(term1,11)
term2 = moving_sum(term2,11)
term3 = moving_sum(term3,11)
term4 = moving_sum(term4,11)
term5 = moving_sum(term5,11)
lons = lons[5:360-5]

pl.plot(lons,term1,'r',linewidth=1.5,alpha=0.8,label='(Vp)T.Qc')
pl.plot(lons,term2,'b',linewidth=1.5,alpha=0.8,label='(Qp)T.Vc')
pl.plot(lons,term3,'g',linewidth=1.25,alpha=0.6,label='(Vp.Qp)T')
pl.plot(lons,term4,'DarkMagenta',linewidth=2,alpha=1,label='Total')
pl.plot(lons,term5,'k--',linewidth=1.25,alpha=0.6,label='residual')
pl.grid()
pl.xlabel('Longitude')
pl.ylabel('Moisture flux trend')
pl.xlim(0,359)
pl.show()
"""
