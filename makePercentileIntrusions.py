from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from UnPickle import *
from scipy import interpolate,stats
from LambertProjector import *

import matplotlib.pyplot as pl
import numpy as np
import os

def myround(x, base=5):
	return int(base * round(float(x)/base))

def interpolateFluxFile(vq):
        n0,n1 = vq.shape
        print 'Interpolating data to 1 degree resolution...'
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

def moistureSeries(vq,d,years,mio=(11,3),thresh=0):
        F    = {}
        for year in years:
                F[year] = []
        for ii in range(len(d)):
                date = d[ii]
                year,month,day,hour = date
                if (month >= mio[0]) or (month <= mio[1]):
                        if month >= 10: year = year + 1
                        try:
                                flux = np.ma.masked_where(vq[ii]<thresh,vq[ii]).sum()
                                if type(flux) != np.ma.core.MaskedConstant: F[year].append(flux)
                        except:
                                pass
        f = np.array([np.sum(F[key]) for key in F.keys()])
        return f

def seasonalVar(years,Season='NDJFM',blat=80.,Source='CAM4xCO2'):
        # Seasonal mean field
        proj = LambertProjector(boundinglat=blat,resolution=200.)
	latm = np.tile(proj.lat,(len(years),1,1))
        mmds = MMDataServer(Field='Ts',Source=Source)
        s    = proj(np.array([mmds.getSeason(year,Season=Season) for year in years]),mmds.lon,mmds.lat)
	s    = np.ma.masked_where(latm<blat,s)
        s    = s.mean(axis=-1).mean(axis=-1)
        return s

def intrusionVar(Q,D,years,mio=(11,3)):
        INT = {}
        for year in years:
                INT[year] = []
        for ii in range(len(D)):
                date = D[ii][0]
                year,month,day,hour = date
                if (month >=mio[0]) or (month <= mio[1]):
                        if month >= 10: year = year + 1
                        try:
                                INT[year].append(Q[ii])
                        except:
                                pass
        # Mean number, flux and intensity per event
        N = np.array([len(INT[key]) for key in INT.keys()])
        E = np.array([sum([sum([sum(q1) for q1 in q]) for q in INT[year]]) for year in INT.keys()])
        I = E/N
        # Total duration
        dur   = np.array([sum([len(q) for q in INT[year]]) for year in INT.keys()])
        return N,E,I

def testInt(G,D):
	gs = {}
	for i in range(len(D)):
		for j in range(len(D[i])):
			try:
				gs[D[i][j]].append(G[i][j])
			except:
				gs[D[i][j]] = []
				gs[D[i][j]].append(G[i][j])
	Gs = {}
	for i in gs:
		if len(gs[i]) > 1:
			inter = list(set.intersection(*map(set,gs[i])))
			if inter != []:	
				Gs[i] = [[inte for inte in inter] for kk in range(len(gs[i])-1)]
	print len(gs),len(Gs)
	for i in Gs:
		print i,gs[i]

"""
Expt     = '16xCO2'
intfile  = '/mnt/climstorage/cian/intrusions/CAM%s_intrusions.1901-1919.NDJFM.6x6hr.9deg.%s.6dt.20.5.p' % (Expt,420)
G,Q,D    = unpick(intfile)
fluxfile = '/mnt/climstorage/cian/scripts/newfluxfiles/70/CAM%s/CAM%s.moist.%s-%s.0-1000hPa.70N.%s.p' % (Expt,Expt,1901,1919,'NDJFM')
vq,dates = unpick(fluxfile)
vq       = interpolateFluxFile(vq)

pl.plot(vq[dates.index((1907,1,31,0)),:])
pl.plot([112,123],[420,420],'k--')
pl.show()
"""

"""
Expt     = 'ERAInt'
flux     = 200
years    = range(1980,2012+1,1)
fluxfile = '/mnt/climstorage/cian/scripts/newfluxfiles/70/%s/%s.moist.%s-%s.0-1000hPa.70N.%s.p' % (Expt,Expt,years[0],years[-1],'NDJFM')
intfile  = '/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.filtered.p' % (Expt,years[0],years[-1],flux)

vq,dates = unpick(fluxfile)
vq       = interpolateFluxFile(vq)
G,Q,D    = unpick(intfile)
N,E,I    = intrusionVar(Q,D,years,(11,3))
s        = seasonalVar(years,'NDJFM',80,Source='%s' % (Expt))
F        = moistureSeries(vq,dates,years,mio=(11,3),thresh=flux)

print stats.percentileofscore(vq.reshape(-1)[np.where(vq.reshape(-1)>0)],flux)
print np.ma.masked_where(vq<flux,vq).sum()/np.ma.masked_where(vq<0,vq).sum()
print E.sum()/np.ma.masked_where(vq<0,vq).sum()
print np.corrcoef(E,s)[0][1]
print np.corrcoef(F,s)[0][1]
"""

# Input
blat        = 60
Expt        = 'CAM8xCO2'
years       = range(1901,1919+1,1)
LevRange    = (0,1000)
percentiles = range(30,98+1,1)
vq,dates    = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.moist.%s-%s.%s-%shPa.%sN.%s.p' % (blat,Expt,Expt,years[0],years[-1],LevRange[0],LevRange[1],blat,'NDJFM'))
vq          = interpolateFluxFile(vq)
# Postive fluxes
vqpos = vq.reshape(-1)[np.where(vq.reshape(-1)>0)]
# Percentile fluxes (rounded to nearest 5 Tg day**-1 deg**-1)
Flux  = [myround(stats.scoreatpercentile(vqpos,pi),base=1) for pi in percentiles]
Fluxp = [        stats.scoreatpercentile(vqpos,pi)         for pi in percentiles]
# Seasonal variation
s  = seasonalVar(years,'NDJFM',blat,Source='%s' % (Expt))
# Full flux
F     = moistureSeries(vq,dates,years,mio=(11,3),thresh=0)
rE,rR = [],[]
rF,rT = [],[]
mE,sE = [],[]
p,f   = [],[]
for ii in range(len(percentiles)):
	flux    = Flux[ii]
	print flux
	intfile = '/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.%sN.filtered.%sN.p' % (Expt,years[0],years[-1],flux,blat,blat+10)
	G,Q,D   = unpick(intfile)
	print len(G)
	N,E,I   = intrusionVar(Q,D,years,mio=(11,3))
	R       = F - E
	T       = moistureSeries(vq,dates,years,mio=(11,3),thresh=flux)
	slope, intercept, r_value, p_value, std_err = stats.linregress(E,s)

	mE.append(slope)
	sE.append(E.std())
	rE.append(np.corrcoef(E,s)[0][1])
	rR.append(np.corrcoef(R,s)[0][1])
	rT.append(np.corrcoef(T,s)[0][1])
	p.append(100.*E.sum()/F.sum())	
	f.append(flux)

fig,ax1 = pl.subplots(num=1)
ax2     = ax1.twinx()
ax1.plot(percentiles, rE, 'r'  , linewidth=1.55, alpha=0.8, label='E vs. Ts')
ax1.plot(percentiles, rR, 'k'  , linewidth=1.55, alpha=0.8, label='R vs. Ts')
ax1.plot(percentiles, rT, 'g'  , linewidth=1.55, alpha=0.8, label='T vs. Ts')
ax2.plot(percentiles,  p, 'b--', linewidth=1.55, alpha=0.8, label='E/F')
ax1.set_xlabel('Flux percentile [%]')
ax1.set_ylabel('Correlation coefficient')
ax2.set_ylabel('Proportion of northward transport [%]')
pl.title(Expt)
ax1.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax2.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
ax1.set_xlim(30,98)
ax1.set_ylim(-0.8,0.8)
#ax2.set_ylim(0,40)
pl.savefig('/mnt/climstorage/cian/scripts/figs/CAM/%s.percentiles.%s-%s.NDJFM.%sN.pdf' % (Expt,years[0],years[-1],blat), format='pdf')
pl.show()
