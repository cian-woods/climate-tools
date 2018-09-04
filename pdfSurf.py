from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from scipy import stats,interpolate
from LambertProjector import *

import matplotlib.pyplot as pl
import numpy as np
import glob

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

def bivarPDF(x,y,xedges,yedges,norm=True):
	H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=[[xedges[0],xedges[-1]],[yedges[0],yedges[-1]]],normed=norm,weights=None)
	xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
	return xedges,yedges,np.rollaxis(H,1,0)

# Parameters
Modelsno    = ['GFDL-CM3']	# prw daily data not functioning...
Models      = [g[9:] for g in glob.glob('../rcp85/*') if g[9:] not in Modelsno]
bins2       = np.linspace(80,350,36)
bins1       = np.linspace(4,10,35)
Season      = 'DJF'
years1      = range(1982,2005+1,1)
years2      = range(2075,2099+1,1)
proj        = LambertProjector(boundinglat=80,resolution=400.)

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
L1,L2   = [],[]
for Model in ['BNU-ESM']:#Models:
	d1 = cmipDataServer(Field='prw',Source=Model,ExpType='historical',LevType='surface',DataFreq='day')
	d2 = cmipDataServer(Field='prw',Source=Model,ExpType='rcp85',LevType='surface',DataFreq='day')
        d3 = cmipDataServer(Field='rlds' ,Source=Model,ExpType='historical',LevType='surface',DataFreq='day')
        d4 = cmipDataServer(Field='rlds' ,Source=Model,ExpType='rcp85',LevType='surface',DataFreq='day')
	S1,S3 = np.array([]),np.array([])
	for year in years1:
		s1 = proj(d1.getDataSnaps(Year=year,Season=Season),d1.lon,d1.lat).reshape(-1)
		S1 = np.append(S1,s1,axis=0)
		s3 = proj(d3.getDataSnaps(Year=year,Season=Season),d3.lon,d3.lat).reshape(-1)
		S3 = np.append(S3,s3,axis=0)
	d1.closeFiles(d1.Handles)
	d3.closeFiles(d3.Handles)
	S2,S4 = np.array([]),np.array([])
	for year in years2:
	        s2 = proj(d2.getDataSnaps(Year=year,Season=Season),d2.lon,d2.lat).reshape(-1)
	        S2 = np.append(S2,s2,axis=0)
                s4 = proj(d4.getDataSnaps(Year=year,Season=Season),d4.lon,d4.lat).reshape(-1)
                S4 = np.append(S4,s4,axis=0)
	d2.closeFiles(d2.Handles)
	d4.closeFiles(d4.Handles)
	#H1,edges = np.histogram(S1,bins=bins,range=None,normed=True,weights=None,density=None)
	#H2,edges = np.histogram(S2,bins=bins,range=None,normed=True,weights=None,density=None)
	S1t,S2t = [S1[i] for i in range(len(S1)) if bins1[0] <= S1[i] <= bins1[-1]],[S2[i] for i in range(len(S2)) if bins1[0] <= S2[i] <= bins1[-1]]
	S3t,S4t = [S3[i] for i in range(len(S1)) if bins1[0] <= S1[i] <= bins1[-1]],[S4[i] for i in range(len(S2)) if bins1[0] <= S2[i] <= bins1[-1]]
	S1,S2   = S1t[:],S2t[:]
	S3,S4   = S3t[:],S4t[:]
	slope1,intercept1,r_value1,p_value1,std_err1 = stats.linregress(S1,S3)
	slope2,intercept2,r_value2,p_value2,std_err2 = stats.linregress(S2,S4)
	xe,ye,H1 = bivarPDF(S1,S3,bins1,bins2,norm=True)
	xe,ye,H2 = bivarPDF(S2,S4,bins1,bins2,norm=True)
	HH1.append(H1)
	HH2.append(H2)
	L1.append([slope1*x + intercept1 for x in bins1])
	L2.append([slope2*x + intercept2 for x in bins1])
HH1,HH2 = np.array(HH1).mean(axis=0),np.array(HH2).mean(axis=0)
L1,L2   = np.array(L1).mean(axis=0),np.array(L2).mean(axis=0)

cseqa = np.arange(-0.02,0.02+0.004,0.004)

pl.figure(1)
pl.contourf(xe,ye,H1,12,cmap=pl.cm.OrRd,extend='max')
pl.colorbar()
pl.plot(bins1,L1,'b')
pl.plot(bins1,L2,'r')
pl.xlabel('Precipitable water')
pl.ylabel('Downward radiation')
pl.figure(2)
pl.contourf(xe,ye,H2,12,cmap=pl.cm.OrRd,extend='max')
pl.colorbar()
pl.plot(bins1,L1,'b')
pl.plot(bins1,L2,'r')
pl.xlabel('Precipitable water')
pl.ylabel('Downward radiation')
pl.figure(3)
pl.contourf(xe,ye,H2-H1,cseqa,cmap=pl.cm.RdBu_r,extend='both')
pl.colorbar()
pl.plot(bins1,L1,'b')
pl.plot(bins1,L2,'r')
pl.xlabel('Precipitable water')
pl.ylabel('Downward radiation')
pl.show()

"""
# Plot
pl.plot(bins,HH1,'b',linewidth=2,alpha=0.7,label='historical: 1981 - 2005')
pl.plot(bins,HH2,'r',linewidth=2,alpha=0.7,label='RCP8.5: 2075 - 2100')
#pl.plot(bins,HH0,'g',linewidth=2,alpha=0.7)
pl.title('PDFs of %s' % (d1.Field))
pl.xlabel('%s [%s]' % (d1.long_name,d1.units))
pl.ylabel('Normalised probability')
pl.xlim(bins[0],bins[-1])
pl.grid()
lg = pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.savefig('figs/pdf/%s.pdf' % (d1.Field),format='pdf')
pl.show()
#pl.close()
"""
