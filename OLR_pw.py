from ReanalysisDataServer import *
from RegressFlux import *
from scipy import stats
import matplotlib.pyplot as pl

# Attributes
rf     = RegressFlux(Type='back1')
proj   = LambertProjector(boundinglat=85,resolution=100.)
pwds   = DataServer(Field='pw',LevType='surface_analysis',Source='ERAInt',LatRange=(70,90))
trds   = DataServer(Field='flds',LevType='surface_forecast',Source='ERAInt',LatRange=(70,90))
years1 = range(1980,1985+1,1)
years2 = range(2010,2015+1,1)
# Bins
#xedges,yedges  = np.linspace(-0.4,0.4,30),np.linspace(-5,5,40)		# T2
xedges,yedges = np.linspace(-0.4,0.4,40),np.linspace(-10,10,40)		# flds
#xedges,yedges = np.linspace(0,5,40),np.linspace(-200,-135,50)		# ttr
HH1,HH2 = [],[]
xx,yy   = [],[]
xxs     = np.linspace(-0.5,0.5,50)
for year in years1:
	print year
	# Get data
	pw = proj(pwds.getDataSnaps(Year=year,Season='DJF'),pwds.lon,pwds.lat)#.mean(axis=1).mean(axis=1)
	tr = proj(trds.getDataSnaps(Year=year,Season='DJF'),trds.lon,trds.lat)#.mean(axis=1).mean(axis=1)
	pw0 = pw[0:-1,:,:].reshape(-1)
	dpw = np.diff(pw,axis=0).reshape(-1)
	dtr = np.diff(tr,axis=0).reshape(-1)
	xx  = xx + list(dpw)
	yy  = yy + list(dtr)
	#dpw0 = dpw#[np.where(pw0<2)]
	#dtr0 = dtr#[np.where(pw0<2)]
	# Histogram
	#x,y,H1 = rf.bivarPDF(dpw0,dtr0,xedges=xedges,yedges=yedges,norm=True)
	#HH1.append(H1)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(xx,yy)
line1 = [slope1*ix + intercept1 for ix in xxs]
xx,yy = [],[]
for year in years2:
        print year
        # Get data
        pw = proj(pwds.getDataSnaps(Year=year,Season='DJF'),pwds.lon,pwds.lat)#.mean(axis=1).mean(axis=1)
        tr = proj(trds.getDataSnaps(Year=year,Season='DJF'),trds.lon,trds.lat)#.mean(axis=1).mean(axis=1)
        pw0 = pw[0:-1,:,:].reshape(-1)
        dpw = np.diff(pw,axis=0).reshape(-1)
        dtr = np.diff(tr,axis=0).reshape(-1)
        xx  = xx + list(dpw)
        yy  = yy + list(dtr)
        #dpw0 = dpw#[np.where(pw0<2)]
        #dtr0 = dtr#[np.where(pw0<2)]
        # Histogram
        #x,y,H2 = rf.bivarPDF(dpw0,dtr0,xedges=xedges,yedges=yedges,norm=True)
        #HH2.append(H2)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(xx,yy)
line2 = [slope2*ix + intercept2 for ix in xxs]

pl.plot(xxs,line1,'k',linewidth=1.75,alpha=0.4)
pl.plot(xxs,line2,'r',linewidth=1.75,alpha=0.4)
pl.grid()
pl.show()

"""
HH1 = np.array(HH1).sum(axis=0)
HH2 = np.array(HH2).sum(axis=0)
cseqf,cseqa = np.arange(0,5.6+0.4,0.4),np.arange(-1,1+0.2,0.2)
pl.figure(1)
pl.contourf(y,x,HH1,cseqf,cmap=pl.cm.RdBu_r,extend='max')
pl.colorbar()
pl.xlabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (trds.long_name,trds.units))
pl.ylabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (pwds.long_name,pwds.units))
pl.figure(2)
pl.contourf(y,x,HH2,cseqf,cmap=pl.cm.RdBu_r,extend='max')
pl.colorbar()
pl.xlabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (trds.long_name,trds.units))
pl.ylabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (pwds.long_name,pwds.units))
pl.figure(3)
pl.contourf(y,x,HH2-HH1,cseqa,cmap=pl.cm.RdBu_r,extend='max')
pl.colorbar()
pl.xlabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (trds.long_name,trds.units))
pl.ylabel(r'$\Delta$%s [%s {6 hr}$^{-1}$]' % (pwds.long_name,pwds.units))
pl.show()
"""
