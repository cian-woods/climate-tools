from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from LambertProjector import *
from drawCoastlinesNoRivers import drawCoastlinesNoRivers
from netCDF4 import Dataset
from UnPickle import *
from TrackDataServer import DataServer as trDataServer
from scipy import interpolate

import matplotlib.pyplot as pl
import numpy as np

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
	vq   = interpolateFluxFile(vq)
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
        proj = LambertProjector(boundinglat=blat,resolution=80.)
        mmds = MMDataServer(Field='Ts',Source=Source)
        s    = np.array([mmds.getSeason(year,Season=Season) for year in years])
        s    = proj(s,mmds.lon,mmds.lat).mean(axis=-1).mean(axis=-1)
        return s,proj.x,proj.y

def plotClim(years,Field,Season,blat):
	proj = LambertProjector(boundinglat=blat,resolution=80.)
	# 4xCO2
	mmds = MMDataServer(Field=Field,Source='CAM4xCO2')
	s4   = proj(np.array([mmds.getSeason(year,Season=Season) for year in years]).mean(axis=0),mmds.lon,mmds.lat)
        # 8xCO2
        mmds = MMDataServer(Field=Field,Source='CAM8xCO2')
        s8   = proj(np.array([mmds.getSeason(year,Season=Season) for year in years]).mean(axis=0),mmds.lon,mmds.lat)
        # 16xCO2
        mmds = MMDataServer(Field=Field,Source='CAM16xCO2')
        s16  = proj(np.array([mmds.getSeason(year,Season=Season) for year in years]).mean(axis=0),mmds.lon,mmds.lat)

	pl.figure(1)
	cf   = pl.contourf(proj.x,proj.y,s4,14,cmap=pl.cm.OrRd,extend='both')
	cbar = pl.colorbar(cf)
	proj.m.drawparallels([70,80],latmax=90)

        pl.figure(2)
        cf   = pl.contourf(proj.x,proj.y,s8-s4,np.arange(2,20+2,2),cmap=pl.cm.OrRd,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawparallels([70,80],latmax=90)

        pl.figure(3)
        cf   = pl.contourf(proj.x,proj.y,s16-s8,np.arange(2,20+2,2),cmap=pl.cm.OrRd,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawparallels([70,80],latmax=90)
        pl.show()

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
        E = np.array([sum([sum([sum(q1) for q1 in q]) for q in INT[year]])/4. for year in INT.keys()])
        I = E/N
        # Total duration
        dur   = np.array([sum([len(q) for q in INT[year]]) for year in INT.keys()])
	return N,E,I

def plotMap(x,y,field,cseq,cmap,extend,cbarlabel,title,proj,savename=None,P=None,sxy=None,Tang=None,stip=None):
        cf = pl.contourf(x,y,field,cseq,cmap=cmap,extend=extend)
        if P != None:
                cl = pl.contour(x,y,P,levels=[0.20],colors='k',linewidths=2,alpha=0.8)
        if sxy != None:
                sx,sy = sxy
                pl.plot(sx,sy,'k.',alpha=0.3)
        if Tang != None:
                sx,sy,u,v = proj.x,proj.y,Tang[:,:,0],Tang[:,:,1]
                Q = pl.quiver(sx,sy,u,v,units='inches',scale=2,\
                        scale_units='inches',headwidth=3,headlength=5,headaxislength=4.5,pivot='tail')
                qk = pl.quiverkey(Q, 0.2, 1.02, 1, '%s%s' % (100,'%'), labelpos='W',fontproperties={'weight': 'bold'})
        if stip != None:
                pl.plot(stip[0],stip[1],'g.',markersize=8,alpha=0.5)
        cbar   = pl.colorbar(cf)
        cbar.set_label(cbarlabel)
        pl.title(title)  
        #drawCoastlinesNoRivers(proj.m,color='0.4')
        proj.m.drawparallels([70,80],latmax=90)
        if savename != None:
                pl.savefig(savename)
                pl.close()
        else:
                pl.show()

#Expt,flux  = 'CAM16xCO2',420
YearRange  = (1901,1919)
years      = range(YearRange[0],YearRange[1]+1,1)
#proj       = LambertProjector(boundinglat=70,resolution=350.)
#d          = trDataServer(Source='%s' % (Expt),Type='fwrd',steps=21)
#G,Q,D      = unpick('/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.NDJFM.6x6hr.9deg.%s.6dt.20.5.filtered.p' % (Expt,YearRange[0],YearRange[1],flux))

plotClim(years,'pw','NDJFM',40)

"""
LON,LAT,P        = d.getIntrusionTrajectories(G,D)
nn,Ntot,Ttot,x,y = d.density(LON,LAT,proj)
Ntot,x,y         = d.interp2d(Ntot,x,y,6,kind='linear')

plotMap(x,y,Ntot,cseq=14,cmap=pl.cm.OrRd,extend='both',cbarlabel='',title='',proj=proj,savename=None,P=None,sxy=None,Tang=Ttot,stip=None)
"""
