from cmipDataServer import *
from scipy import interpolate

import matplotlib.pyplot as pl
import numpy as np
import glob

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

lats   = np.arange(65,70+0.25,0.25)
lons   = range(0,360,1)
years1 = range(1981,2005+1,1)
years2 = range(2075,2095+1,1)
notM   = ['FGOALS-g2','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M']
Models = [g[9:] for g in glob.glob('../rcp85/*') if g[9:] not in notM]
lonm   = np.array(list(np.zeros(310)) + list(np.ones(30)) + list(np.zeros(20)))

D1,D2 = [],[]
for Model in Models:
	ds1 = DataServer(Field='prw',Source=Model,ExpType='historical',LevType='surface',DataFreq='mon',LevRange=(50,50))
	ds2 = DataServer(Field='prw',Source=Model,ExpType='rcp85',LevType='surface',DataFreq='mon',LevRange=(50,50))

	s1 = np.array([ds1.getDataSnaps(Year=year,Season='DJF').mean(axis=0) for year in years1]).mean(axis=0)
	s2 = np.array([ds2.getDataSnaps(Year=year,Season='DJF').mean(axis=0) for year in years2]).mean(axis=0)

	s1 = interpolateND(s1,ds1.lat,lats,axis=0)
	s2 = interpolateND(s2,ds2.lat,lats,axis=0)
	s1 = np.append(s1,(s1[:,-1])[:,np.newaxis],axis=1)
	s2 = np.append(s2,(s2[:,-1])[:,np.newaxis],axis=1)
	l1 = np.append(ds1.lon,ds1.lon[0]+360)
	l2 = np.append(ds2.lon,ds2.lon[0]+360)
        s1 = interpolateND(s1,l1,lons,axis=1)
        s2 = interpolateND(s2,l2,lons,axis=1)

	m1,c1 = np.polyfit(lats,s1,1)
	m2,c2 = np.polyfit(lats,s2,1)

	m1,m2 = np.ma.masked_where(lonm==1,m1),np.ma.masked_where(lonm==1,m2)

	pl.plot(lons,m1,'b')
	pl.plot(lons,m2,'r')
	pl.plot([0,359],[0,0],'k--')
	pl.xlim(0,359)
	pl.title(Model)
	pl.show()

"""
	D1.append(m1)
	D2.append(m2)

D1,D2 = np.array(D1).mean(axis=0),np.array(D2).mean(axis=0)

pl.plot(D1,'b')
pl.plot(D2,'r')
pl.plot([0,359],[0,0],'k--')
pl.xlim(0,359)
pl.show()
"""
