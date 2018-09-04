from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from scipy import interpolate

import matplotlib.pyplot as pl
import numpy as np
import glob

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

lats = [70]
lons = range(2,10+1,1)
levs = [250,500,850]

Models = [g[9:] for g in glob.glob('../rcp85/*')]
print Models
bins   = np.linspace(-30,30,60)
Season = 'DJF'
plevs  =  (200,1000)
years1 = range(1981,2005+1,1)
years2 = range(2075,2100+1,1)

"""
# ERAInt PDF
S0 = [np.array([]) for i in range(len(levs))]
d0 = reDataServer(Field='V',LevType='plev',LevRange=plevs,LatRange=(65,80))
for year in years1:
	print year
	s0 = d0.getDataSnaps(Year=year,Season=Season,dailymean=True).squeeze()
	s0 = interpolateND(s0,d0.lev,levs,axis=1)
	s0 = interpolateND(s0,d0.lat,lats,axis=2)
	s0 = interpolateND(s0,d0.lon,lons,axis=3).squeeze().reshape((1,3,-1)).squeeze()
	S0 = np.append(S0,s0,axis=1)
HH0 = [np.histogram(S0[i],bins=bins,range=None,normed=True,weights=None,density=None)[0] for i in range(len(levs))]
"""

# Model PDFS
HH1,HH2 = [],[]
var1    = []
for Model in Models:
	d1 = cmipDataServer(Field='va',Source=Model,ExpType='historical',LevType='surface',DataFreq='day',LevRange=plevs)
	d2 = cmipDataServer(Field='va',Source=Model,ExpType='rcp85',LevType='surface',DataFreq='day',LevRange=plevs)
	S1 = [np.array([]) for i in range(len(levs))]
	for year in years1:
		s1 = d1.getDataSnaps(Year=year,Season=Season).squeeze()
		s1 = interpolateND(s1,d1.lev,levs,axis=1)
		s1 = interpolateND(s1,d1.lat,lats,axis=2)
		s1 = interpolateND(s1,d1.lon,lons,axis=3).squeeze().reshape((1,3,-1)).squeeze()
		S1  = np.append(S1,s1,axis=1)
	S2 = [np.array([]) for i in range(len(levs))]
	for year in years2:
	        s2 = d2.getDataSnaps(Year=year,Season=Season).squeeze()
	        s2 = interpolateND(s2,d2.lev,levs,axis=1)
	        s2 = interpolateND(s2,d2.lat,lats,axis=2)
		s2 = interpolateND(s2,d2.lon,lons,axis=3).squeeze().reshape((1,3,-1)).squeeze()
	        S2  = np.append(S2,s2,axis=1)
	var.append((S2.std()-S1.std())/S1.std())
	H1 = [np.histogram(S1[i],bins=bins,range=None,normed=True,weights=None,density=None)[0] for i in range(len(levs))]
	H2 = [np.histogram(S2[i],bins=bins,range=None,normed=True,weights=None,density=None)[0] for i in range(len(levs))]
	HH1.append(H1)
	HH2.append(H2)
HH1,HH2 = np.array(HH1).mean(axis=0),np.array(HH2).mean(axis=0)
bins = (bins[1:] + bins[0:-1])/2.

# Plot
for i in range(len(levs)):
	pl.plot(bins,HH1[i],'b',linewidth=2,alpha=0.7)
	pl.plot(bins,HH2[i],'r',linewidth=2,alpha=0.7)
	pl.plot(bins,HH0[i],'g',linewidth=2,alpha=0.7)
	pl.title('PDFs of V at %shPa' % (levs[i]))
	pl.xlabel('Velocity [m s$ ^{-1}$]')
	pl.ylabel('Normalised probability')
	pl.xlim(bins[0],bins[-1])
	pl.savefig('figs/pdf/v.%s-%sE.%shPa.pdf' % (lons[0],lons[-1],levs[i]),format='pdf')
	pl.close()
