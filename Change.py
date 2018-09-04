from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from stipling import stipling
from drawCoastlinesNoRivers import *

import sys
import glob
import matplotlib.pyplot as pl

Field    = str(sys.argv[1])
LevType  = str(sys.argv[2])
c1,c2,dc = float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])
cseq     = np.arange(c1,c2+dc,dc)
Models   = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')]
print Models
#Models   = [	'ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2',\
#		'CCSM4', 'CMCC-CESM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2',\
#		'GFDL-CM3', 'GFDL-ESM2G', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR',\
#		'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M'	]
proj    = LambertProjector(boundinglat=65,resolution=150.)
x70,y70 = proj.m(np.linspace(0,360,1000), np.zeros(1000)+70 )
x80,y80 = proj.m(np.linspace(0,360,1000), np.zeros(1000)+80 )
years1   = range(2070,2095+1,1)
years2   = range(1981,2005+1,1)

S,G1,G2,M = [],[],[],[]
for Source in Models:
	try:
		ds1 = cmipDataServer(Field=Field,Source=Source,ExpType='rcp85',LevType=LevType,DataFreq='mon',LevRange=(250,250))
		s1  = np.array([proj(ds1.getDataSnaps(Year=year,Season='NDJF').squeeze().mean(axis=0),ds1.lon,ds1.lat) for year in years1]).mean(axis=0)/100.
        	ds2 = cmipDataServer(Field=Field,Source=Source,ExpType='rcp85',LevType=LevType,DataFreq='mon',LevRange=(250,250))
        	s2  = np.array([proj(ds2.getDataSnaps(Year=year,Season='NDJF').squeeze().mean(axis=0),ds2.lon,ds2.lat) for year in years2]).mean(axis=0)/100.
		S.append(s1-s2)
		G1.append(s1)
		G2.append(s2)
		cf   = pl.contourf(proj.x,proj.y,s1-s2,cseq,cmap=pl.cm.RdBu_r,extend='both')
		cbar = pl.colorbar(cf)
		proj.m.drawcoastlines()
		proj.m.drawparallels([70,80],latmax=90)
		pl.title(Source)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/bias/%s/%s.diff.pdf' % (Field,Source),format='pdf')
		pl.close()
		cf   = pl.contourf(proj.x,proj.y,s2,20,cmap=pl.cm.RdBu_r,extend='both')
		cbar = pl.colorbar(cf)
		proj.m.drawcoastlines()
		proj.m.drawparallels([70,80],latmax=90)
		pl.title(Source)
		pl.savefig('/mnt/climstorage/cian/scripts/figs/bias/%s/%s.full.pdf' % (Field,Source),format='pdf')
		pl.close()
		M.append(Source)
	except:
		print '\n%s did not finish...\n' % (Source)

S  = np.array(S)
G1 = np.array(G1)
G2 = np.array(G2)
NS,NG1,NG2 = len(S),len(G1),len(G2)
stipxs,stipys = stipling(S,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
stipxg1,stipyg1 = stipling(G1,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
stipxg2,stipyg2 = stipling(G2,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
S = S.mean(axis=0)
G1 = G1.mean(axis=0)
G2 = G2.mean(axis=0)

cf    = pl.contourf(proj.x,proj.y,S,cseq,cmap=pl.cm.coolwarm,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Sea-level pressure change [hPa]')
plt.plot(stipxs,stipys,'k.',alpha=0.5)
pl.plot(x70,y70,'k-',linewidth=1,alpha=0.65,zorder=9)
pl.plot(x80,y80,'k-',linewidth=1,alpha=0.65,zorder=9)
drawCoastlinesNoRivers(proj.m,linewidth=0.5,color='0.25')
proj.m.drawparallels([70,80],latmax=90)
pl.title('all %s' % (NS))
pl.savefig('/mnt/climstorage/cian/scripts/figs/bias/%s/%s.diff.pdf' % (Field,'all'),format='pdf')
pl.close()

cf    = pl.contourf(proj.x,proj.y,G1,20,cmap=pl.cm.RdBu_r,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Sea-level pressure [hPa]')
plt.plot(stipxg1[::3,::3],stipyg1[::3,::3],'k.',alpha=0.5)
proj.m.drawcoastlines()
proj.m.drawparallels([70,80],latmax=90)
pl.title('all %s' % (NG1))
pl.savefig('/mnt/climstorage/cian/scripts/figs/bias/%s/%s.full.rcp85.pdf' % (Field,'all'),format='pdf')
pl.close()

cf    = pl.contourf(proj.x,proj.y,G2,20,cmap=pl.cm.RdBu_r,extend='both')
cbar  = pl.colorbar(cf)
cbar.set_label('Sea-level pressure [hPa]')
plt.plot(stipxg2[::3,::3],stipyg2[::3,::3],'k.',alpha=0.5)
proj.m.drawcoastlines()
proj.m.drawparallels([70,80],latmax=90)
pl.title('all %s' % (NG2))
pl.savefig('/mnt/climstorage/cian/scripts/figs/bias/%s/%s.full.hist.pdf' % (Field,'all'),format='pdf')
pl.close()
