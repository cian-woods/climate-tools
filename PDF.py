from ReanalysisDataServer import *
from LambertProjector import *
from RegressFlux import *
from scipy import stats

import numpy as np
import matplotlib.pyplot as pl

# Attributes
rf     = RegressFlux(Type='back1')
dpw    = DataServer(Field='pw',LevType='surface_analysis',LatRange=(70,90))
dnl    = DataServer(Field='fls',LevType='surface_forecast',LatRange=(70,90))
proj   = LambertProjector(boundinglat=80,resolution=200.)
#years1 = [1982,1989,1998,1999,2003,2004]#range(1981,1998+1,1)
#years2 = [1985,1990,1991,1995,2000,2002]#range(1999,2016+1,1)
years1 = range(1981,1998+1,1)
years2 = range(1999,2016+1,1)

# Get data
pw1 = dpw.getDataSnaps(Year=years1[0],Season='DJF')
pw2 = dpw.getDataSnaps(Year=years2[0],Season='DJF')
nl1 = dnl.getDataSnaps(Year=years1[0],Season='DJF')
nl2 = dnl.getDataSnaps(Year=years2[0],Season='DJF')
for year in years1[1:]:
	print year
	pw1  = np.append(pw1,dpw.getDataSnaps(Year=year,Season='DJF'),axis=0)
	nl1  = np.append(nl1,dnl.getDataSnaps(Year=year,Season='DJF'),axis=0)
for year in years2[1:]:
	print year
	pw2  = np.append(pw2,dpw.getDataSnaps(Year=year,Season='DJF'),axis=0)
	nl2  = np.append(nl2,dnl.getDataSnaps(Year=year,Season='DJF'),axis=0)
pw1,pw2 = proj(pw1,dpw.lon,dpw.lat).reshape(-1),proj(pw2,dpw.lon,dpw.lat).reshape(-1)
nl1,nl2 = proj(nl1,dnl.lon,dnl.lat).reshape(-1),proj(nl2,dnl.lon,dnl.lat).reshape(-1)

# PDFs
xe     = np.linspace(0,4.5,30+1)
ye     = np.linspace(-80,10,30+1)
x,y,H1 = rf.bivarPDF(pw1,nl1,xe,ye,norm=True)
x,y,H2 = rf.bivarPDF(pw2,nl2,xe,ye,norm=True)
H1,H2  = np.rollaxis(H1,1,0),np.rollaxis(H2,1,0)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(pw1,nl1)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(pw2,nl2)
line1 = [xi*slope1 + intercept1 for xi in xe]
line2 = [xi*slope2 + intercept2 for xi in xe]

cseqf,cseqa = np.arange(0,0.015+0.001,0.001),np.arange(-0.004,0.004+0.0005,0.0005)
pl.figure(1)
cf   = pl.contourf(x,y,H2-H1,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label('Frequency')
pl.plot(xe,line1,'b--',linewidth=1.75,alpha=0.6)
pl.plot(xe,line2,'r--',linewidth=1.75,alpha=0.6)
pl.xlabel('%s [%s]' % (dpw.long_name,dpw.units))
pl.ylabel('%s [%s]' % (dnl.long_name,dnl.units))
pl.title('%s to %s minus %s to %s' % (years1[0],years1[-1],years2[0],years2[-1]))
pl.xlim(x[0],x[-1])
pl.ylim(y[0],y[-1])
pl.savefig('figs/regress/pw_fls_change.pdf',format='pdf')
pl.figure(2)
cf   = pl.contourf(x,y,H1,cseqf,cmap=pl.cm.OrRd,extend='max')
cbar = pl.colorbar(cf)
cbar.set_label('Frequency')
pl.xlabel('%s [%s]' % (dpw.long_name,dpw.units))
pl.ylabel('%s [%s]' % (dnl.long_name,dnl.units))
pl.title('%s to %s' % (years1[0],years1[-1]))
pl.savefig('figs/regress/pw_fls_full.pdf',format='pdf')
pl.show()

"""
h1,e1  = np.histogram(s1,50,normed=False)
h2,e2  = np.histogram(s2,50,normed=False)
pl.plot(e1[0:-1],h1,'b',linewidth=1.5,alpha=0.6)
pl.plot(e2[0:-1],h2,'r',linewidth=1.5,alpha=0.6)
pl.xlabel('%s [%s]' % (d.long_name,d.units))
pl.ylabel('Frequency')
pl.show()
"""
