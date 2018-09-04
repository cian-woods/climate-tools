from cmipDataServer import DataServer as cmipDataServer
from LambertProjector import *
from UnPickle import *
from toPick import *
from stipling import stipling

import matplotlib.pyplot as pl
import numpy as np
import glob,os

PLEV   = int(sys.argv[1])
yearsr = range(2075,2100+1,1)
yearsh = range(1981,2005+1,1)
proj   = LambertProjector(boundinglat=10,resolution=80.)

Models = [g[9:] for g in glob.glob('../rcp85/*')]

KEH,KER = [],[]
for Source in Models:
	fnameh = 'eddyfiles/%s/%s.EKE.historical.%s-%s.%shPa.DJF.p' % (PLEV,Source,yearsh[0],yearsh[-1],PLEV)
	fnamer = 'eddyfiles/%s/%s.EKE.rcp85.%s-%s.%shPa.DJF.p' % (PLEV,Source,yearsr[0],yearsr[-1],PLEV)
	print fnamer
	if os.path.isfile(fnameh) == False:
		udsh   = cmipDataServer(Field='ua',Source=Source,ExpType='historical',LevType='plev',DataFreq='day',LevRange=(PLEV,PLEV))
		vdsh   = cmipDataServer(Field='va',Source=Source,ExpType='historical',LevType='plev',DataFreq='day',LevRange=(PLEV,PLEV))
		KEh    = []
		for year in yearsh:
			uh,vh = udsh.getDataSnaps(Year=year,Season='DJF').squeeze(),vdsh.getDataSnaps(Year=year,Season='DJF').squeeze()
			uh,vh = uh-uh.mean(axis=-1)[:,:,np.newaxis],vh-vh.mean(axis=-1)[:,:,np.newaxis]
			KEh.append(0.5*(uh**2 + vh**2).mean(axis=0))
		KEh,lonh,lath = np.array(KEh).mean(axis=0),udsh.lon,udsh.lat
		toPick([KEh,lonh,lath],fnameh)
	else:
		KEh,lonh,lath = unpick(fnameh)
	if os.path.isfile(fnamer) == False:
		udsr   = cmipDataServer(Field='ua',Source=Source,ExpType='rcp85',LevType='plev',DataFreq='day',LevRange=(PLEV,PLEV))
		vdsr   = cmipDataServer(Field='va',Source=Source,ExpType='rcp85',LevType='plev',DataFreq='day',LevRange=(PLEV,PLEV))
		KEr    = []
		for year in yearsr:
			ur,vr = udsr.getDataSnaps(Year=year,Season='DJF').squeeze(),vdsr.getDataSnaps(Year=year,Season='DJF').squeeze()
			ur,vr = ur-ur.mean(axis=-1)[:,:,np.newaxis],vr-vr.mean(axis=-1)[:,:,np.newaxis]
			KEr.append(0.5*(ur**2 + vr**2).mean(axis=0))
		KEr,lonr,latr = np.array(KEr).mean(axis=0),udsr.lon,udsr.lat
		toPick([KEr,lonr,latr],fnamer)
	else:
	        KEr,lonr,latr = unpick(fnamer)

	KEh,KEr = proj(KEh,lonh,lath),proj(KEr,lonr,latr)
	KEH.append(KEh)
	KER.append(KEr)
	pl.figure(1)
	cseq = np.arange(0,100+5,5)
#	cseq = 20
	cf   = pl.contourf(proj.x,proj.y,KEh,cseq,cmap=pl.cm.RdBu_r,extend='both')
	cbar = pl.colorbar(cf)
	proj.m.drawcoastlines()
	proj.m.drawparallels([70,80],latmax=90)
	pl.savefig('figs/eddy/full/%s.pdf' % (Source),format='pdf')
	pl.close()
	pl.figure(2)
	cseq = np.arange(-8,8+1,1)
#	cseq = 20
	cf   = pl.contourf(proj.x,proj.y,KEr-KEh,cseq,cmap=pl.cm.RdBu_r,extend='both')
	cbar = pl.colorbar(cf)
	proj.m.drawcoastlines()
	proj.m.drawparallels([70,80],latmax=90)
	pl.savefig('figs/eddy/diff/%s.pdf' % (Source),format='pdf')
	pl.close()

KEH,KER     = np.array(KEH),np.array(KER)
stipx,stipy = stipling(KER-KEH,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
KEH,KER     = KEH.mean(axis=0),KER.mean(axis=0)

cseq = np.arange(-8,8+1,1)
#cseq = 20
cf   = pl.contourf(proj.x,proj.y,KER-KEH,cseq,cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
pl.plot(stipx[::3,::3],stipy[::3,::3],'k.',alpha=0.5)
proj.m.drawcoastlines()
proj.m.drawparallels([70,80],latmax=90)
pl.savefig('figs/eddy/diff/%s.pdf' % ('all'),format='pdf')
pl.close()




