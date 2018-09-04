from cmipDataServer import DataServer as cmipDataServer
from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from stipling import stipling
from scipy import stats

import numpy as np
import matplotlib.pyplot as pl
import glob

def za(field):
	return field - field.mean(axis=-1)[...,np.newaxis]

proj   = LambertProjector(boundinglat=55,resolution=200.)
Years0 = range(1981,2016+1,1)
Years1 = range(1981,2016+1,1)
Years2 = range(2060,2095+1,1)
plev   = 500

dERA  = MMDataServer(Field='V')
lx    = np.argmin((dERA.lev-plev)**2)
print lx,dERA.lev[lx]
sERA  = np.array([za(dERA.getSeason(Year,Season='DJF')[lx,:,:]) for Year in Years0]).mean(axis=0)
sERA  = proj(sERA,dERA.lon,dERA.lat)

Models  = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')]
ModelsH = ['MRI-ESM1','bcc-csm1-1','bcc-csm1-1-m','BNU-ESM','CMCC-CESM','IPSL-CM5A-LR','IPSL-CM5A-MR','ACCESS1-3']	# High gradient
ModelsL = [Model for Model in Models if Model not in ModelsH]
SH1,SH2   = [],[]
mh_tmp,ml_tmp = [],[]
for Model in ModelsH:
	try:
		dCMIP = cmipDataServer(Field='va',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='mon')
		lx    = np.argmin((dCMIP.lev-plev)**2)
		print dCMIP.lev
		print lx,dCMIP.lev[lx]
		s1    = np.array([za(dCMIP.getSeason(Year,Season='DJF')[lx,:,:]) for Year in Years1]).mean(axis=0)
		s2    = np.array([za(dCMIP.getSeason(Year,Season='DJF')[lx,:,:]) for Year in Years2]).mean(axis=0)
		s1    = proj(s1,dCMIP.lon,dCMIP.lat)
		s2    = proj(s2,dCMIP.lon,dCMIP.lat)
		SH1.append(s1)
		SH2.append(s2)
		mh_tmp.append(Model)
	except:
		print 'Model %s did not complete...' % (Model)
SL1,SL2 = [],[]
for Model in ModelsL:
        try:
                dCMIP = cmipDataServer(Field='va',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='mon')
                lx    = np.argmin((dCMIP.lev-plev)**2)
		print dCMIP.lev
                print lx,dCMIP.lev[lx]
                s1    = np.array([za(dCMIP.getSeason(Year,Season='DJF')[lx,:,:]) for Year in Years1]).mean(axis=0)
		s2    = np.array([za(dCMIP.getSeason(Year,Season='DJF')[lx,:,:]) for Year in Years2]).mean(axis=0)
                s1    = proj(s1,dCMIP.lon,dCMIP.lat)
		s2    = proj(s2,dCMIP.lon,dCMIP.lat)
                SL1.append(s1)
		SL2.append(s2)
		ml_tmp.append(Model)
        except:
	       print 'Model %s did not complete...' % (Model)

print mh_tmp+ml_tmp

nh1,nl1   = len(SH1),len(SL1)
nh2,nl2   = len(SH2),len(SL2)
SH1,SH2   = np.array(SH1),np.array(SH2)
SL1,SL2   = np.array(SL1),np.array(SL2)
SH1m,SH2m = SH1.mean(axis=0),SH2.mean(axis=0)
SL1m,SL2m = SL1.mean(axis=0),SL2.mean(axis=0)

S1,S2         = np.append(SH1,SL1,axis=0),np.append(SH2,SL2,axis=0)
n1,n2         = len(S1),len(S2)
S1a,S2a       = S1 - sERA[np.newaxis,:,:],S2 - S1
S1m,S2m       = S1.mean(axis=0),S2.mean(axis=0)
stipx1,stipy1 = stipling(S1a,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
stipx2,stipy2 = stipling(S2a,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)

bias  = (S1m - sERA)[::6,::6].reshape(-1)
force = (S2m -  S1m)[::6,::6].reshape(-1)
slope, intercept, r_value, p_value, std_err = stats.linregress(bias,force)
xs   = np.linspace(bias.min(),bias.max(),20)
line = [slope*x + intercept for x in xs]

print len(bias),np.corrcoef(bias,force)[0][1]

pl.figure(1)
pl.plot(bias,force,'k.',alpha=0.5)
pl.plot(xs,line,'k--')
pl.grid()
pl.xlabel('Bias')
pl.ylabel('Forced')
pl.title('slope = %s; p = %s; r = %s' % (round(slope,3),round(p_value/2,4),round(r_value,3)))

#cseqf = np.arange(-6,6+1,1)
#cseqa = np.arange(-2,2+0.4,0.4)
cseqf  = np.arange(-9,9+1.5,1.5)
cseqa  = np.arange(-3,3+0.5,0.5)
"""
pl.figure(1)
cf   = pl.contourf(proj.x,proj.y,S1m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (1981-2005): %s hPa (%s models)' % (plev,n1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.hist.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(2)
cf   = pl.contourf(proj.x,proj.y,S2m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (2070-2095): %s hPa (%s models)' % (plev,n2))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.rcp85.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(3)
cf   = pl.contourf(proj.x,proj.y,SH1m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (1981-2005): %s hPa (%s models)' % (plev,nh1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high.hist.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(4)
cf   = pl.contourf(proj.x,proj.y,SL1m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (1981-2005): %s hPa (%s models)' % (plev,nl1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/low.hist.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(5)
cf   = pl.contourf(proj.x,proj.y,SH2m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (1981-2005): %s hPa (%s models)' % (plev,nh2))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high.rcp85.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(6)
cf   = pl.contourf(proj.x,proj.y,SL2m,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5 (2075-2095): %s hPa (%s models)' % (plev,nl2))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/low.rcp85.full.%shPa.pdf' % (plev),format='pdf')

pl.figure(7)
cf   = pl.contourf(proj.x,proj.y,sERA,cseqf,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('ERAInt (1981-2005): %s hPa' % (plev))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/ERAInt.full.%shPa.pdf' % (plev),format='pdf')
"""

pl.figure(8)
cf   = pl.contourf(proj.x,proj.y,S1m-sERA,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cl   = pl.contour(proj.x,proj.y,sERA,levels=[-6.0,-4.5,-3.0,-1.5,1.5,3.0,4.5,6.0],colors='k',linewidths=0.65,alpha=0.55)
cl0  = pl.contour(proj.x,proj.y,sERA,levels=[0],colors='k',linewidths=2,alpha=0.75)
#pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
#pl.clabel(cl0,fmt='%0.1f',colors='k',fontsize=9)
pl.plot(stipx1,stipy1,'k.',alpha=0.5)
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.5,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5-ERAInt: %s hPa (%s models)' % (plev,n1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.bias.%shPa.pdf' % (plev),format='pdf')

"""
pl.figure(9)
cf   = pl.contourf(proj.x,proj.y,SH1m-sERA,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5-ERAInt: %s hPa (%s models)' % (plev,nh1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high.bias.%shPa.pdf' % (plev),format='pdf')

pl.figure(10)
cf   = pl.contourf(proj.x,proj.y,SL1m-sERA,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('CMIP5-ERAInt: %s hPa (%s models)' % (plev,nl1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/low.bias.%shPa.pdf' % (plev),format='pdf')
"""

pl.figure(11)
cf   = pl.contourf(proj.x,proj.y,S2m-S1m,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cl   = pl.contour(proj.x,proj.y,S2m,levels=[-6.0,-4.5,-3.0,-1.5,1.5,3.0,4.5,6.0],colors='k',linewidths=0.65,alpha=0.55)
cl0  = pl.contour(proj.x,proj.y,S2m,levels=[0],colors='k',linewidths=2,alpha=0.75)
pl.plot(stipx2,stipy2,'k.',alpha=0.5)
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('RCP8.5 - HIST: %s hPa (%s - %s models)' % (plev,n2,n1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/all.change.%shPa.pdf' % (plev),format='pdf')

"""
pl.figure(12)
cf   = pl.contourf(proj.x,proj.y,SH2m-SH1m,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('RCP8.5 - HIST: %s hPa (%s - %s models)' % (plev,nh2,nh1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high.change.%shPa.pdf' % (plev),format='pdf')

pl.figure(13)
cf   = pl.contourf(proj.x,proj.y,SL2m-SL1m,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('RCP8.5 - HIST: %s hPa (%s - %s models)' % (plev,nl2,nl1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/low.change.%shPa.pdf' % (plev),format='pdf')

pl.figure(14)
cf   = pl.contourf(proj.x,proj.y,SH1m-SL1m,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('%s High - %s Low: %s hPa' % (nh1,nl1,plev))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high-low.hist.%shPa.pdf' % (plev),format='pdf')

pl.figure(15)
cf   = pl.contourf(proj.x,proj.y,SH2m-SL2m,cseqa,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
cbar.set_label(r'$\Delta{\overline{v}}$ [m s$^{-1}$]')
drawCoastlinesNoRivers(proj.m,linewidth=0.6,color='0.2')
proj.m.drawparallels([70,80],latmax=90)
pl.title('%s High - %s Low: %s hPa' % (nh2,nl2,plev))
pl.savefig('/mnt/climstorage/cian/scripts/figs/eddy/stationary/high-low.rcp85.%shPa.pdf' % (plev),format='pdf')
"""

pl.show()
