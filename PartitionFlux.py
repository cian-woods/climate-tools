from scipy import interpolate,stats
from LambertProjector import *
from UnPickle import *
from matplotlib.colors import LinearSegmentedColormap

import glob
import sys
import matplotlib.pyplot as pl

def moving_average(a,n=3) :
	ret = np.cumsum(a, dtype=float)
	ret[n:] = ret[n:] - ret[:-n]
	return ret[n - 1:]/n

def moving_sum(a,n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:]

def interp(vq):
	n0,n1 = vq.shape
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

def corr():
	# Lambert projector and km per deg longitude at blats
	Re        = 6371
	proj      = LambertProjector(boundinglat=50,resolution=80.)
	xres      = 2*np.pi*Re*np.cos(70*np.pi/180)/360.
	sf        = 1e09/(24.*60*60*xres*1000)
	RH,RR     = [],[]
	SH,SR     = [],[]
	PqH,PqR   = [],[]
	PvH,PvR   = [],[]
	PvqH,PvqR = [],[]
	# Loop through models
	Models = [i[9:] for i in glob.glob('../rcp85/*')]
	for Source in Models:
		print Source
		# Open historical flux files
		hdir         = '/qnap/cian/cmip/scripts'
		vqht,datesh  = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))
		vht,datesh   = unpick('%s/newfluxfiles/%s/%s/%s.mass.1981-2005.600-1000hPa.%sN.DJF.p'  % (hdir,70,Source,Source,70))
		qht          = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))[0]/(40000/9.8)
		vqht,vht,qht = sf*interp(vqht),sf*interp(vht),interp(qht)
		# Open rcp85 flux files
		rdir         = '/mnt/climstorage/cian/scripts'
		vqrt,datesr  = unpick('%s/newfluxfiles/%s/%s/%s.moist.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))
		vrt,datesr   = unpick('%s/newfluxfiles/%s/%s/%s.mass.2075-2100.600-1000hPa.%sN.DJF.p'  % (rdir,70,Source,Source,70))
		qrt          = unpick('%s/newfluxfiles/%s/%s/%s.vapor.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))[0]/(40000/9.8)
		vqrt,vrt,qrt = sf*interp(vqrt),sf*interp(vrt),interp(qrt)
		Sh,Sr,Rh,Rr  = [],[],[],[]
		Pvh,Pvr      = [],[]
		Pqh,Pqr      = [],[]
		Pvh,Pvr      = [],[]
		Pvqh,Pvqr      = [],[]
		print vqht.max()
		for i in range(360):
			vqh = np.array([vqht[j,i] for j in range(len(vqht[:,i])) if vqht[j,i] > 0])
			vh  = np.array([vht[j,i]  for j in range(len(vht[:,i]))  if vqht[j,i] > 0])
			qh  = np.array([qht[j,i]  for j in range(len(qht[:,i]))  if vqht[j,i] > 0])
			vqr = np.array([vqrt[j,i] for j in range(len(vqrt[:,i])) if vqrt[j,i] > 0])
                        vr  = np.array([vrt[j,i]  for j in range(len(vrt[:,i]))  if vqrt[j,i] > 0])
                        qr  = np.array([qrt[j,i]  for j in range(len(qrt[:,i]))  if vqrt[j,i] > 0])		
			slopeh,intercepth,rh,p_valueh,std_errh = stats.linregress(vh,qh)
			sloper,interceptr,rr,p_valuer,std_errr = stats.linregress(vr,qr)
			Sh.append(slopeh)
			Sr.append(sloper)
			Rh.append(rh)
			Rr.append(rr)
#			Pvh.append(np.ma.masked_where(vh<stats.scoreatpercentile(vh,90),vh).sum()/vh.sum())
#			Pvr.append(np.ma.masked_where(vr<stats.scoreatpercentile(vr,90),vr).sum()/vr.sum())
#			Pqh.append(np.ma.masked_where(qh<stats.scoreatpercentile(qh,90),qh).sum()/qh.sum())
#			Pqr.append(np.ma.masked_where(qr<stats.scoreatpercentile(qr,90),qr).sum()/qr.sum())
#			Pvqh.append(np.ma.masked_where(vqh<stats.scoreatpercentile(vqh,90),vqh).sum()/vqh.sum())
#			Pvqr.append(np.ma.masked_where(vqr<stats.scoreatpercentile(vqr,90),vqr).sum()/vqr.sum())
#			Pvh.append(np.ma.masked_where(vh<vh.mean()+2*vh.std(),vh).sum()/vh.sum())
#			Pvr.append(np.ma.masked_where(vr<vr.mean()+2*vr.std(),vr).sum()/vr.sum())
#			Pqh.append(np.ma.masked_where(qh<qh.mean()+2*qh.std(),qh).sum()/qh.sum())
#			Pqr.append(np.ma.masked_where(qr<qr.mean()+2*qr.std(),qr).sum()/qr.sum())
#			Pvqh.append(np.ma.masked_where(vqh<vqh.mean()+2*vqh.std(),vqh).sum()/vqh.sum())
#			Pvqr.append(np.ma.masked_where(vqr<vqr.mean()+2*vqr.std(),vqr).sum()/vqr.sum())
		RH.append(Rh)
		RR.append(Rr)
		SH.append(Sh)
		SR.append(Sr)
#		PvH.append(Pvh)
#		PvR.append(Pvr)
#		PqH.append(Pqh)
#		PqR.append(Pqr)
#		PvqH.append(Pvqh)
#		PvqR.append(Pvqr)
	lons            = np.arange(0,360,1)
	RH,RR,SH,SR     = np.array(RH),np.array(RR),np.array(SH),np.array(SR)
	RHm,RRm,SHm,SRm = RH.mean(axis=0),RR.mean(axis=0),SH.mean(axis=0),SR.mean(axis=0)
	RHs,RRs,SHs,SRs = RH.std(axis=0),RR.std(axis=0),SH.std(axis=0),SR.std(axis=0)
#	PvH,PvR,PqH,PqR = np.array(PvH),np.array(PvR),np.array(PqH),np.array(PqR)
#	PvqH,PvqR       = np.array(PvqH),np.array(PvqR)
#	PvHm,PvRm       = PvH.mean(axis=0),PvR.mean(axis=0)
#	PqHm,PqRm       = PqH.mean(axis=0),PqR.mean(axis=0)
#	PvqHm,PvqRm     = PvqH.mean(axis=0),PvqR.mean(axis=0)
#	PvHs,PvRs       = PvH.std(axis=0),PvR.std(axis=0)
#	PqHs,PqRs       = PqH.std(axis=0),PqR.std(axis=0)
#	PvqHs,PvqRs     = PvqH.std(axis=0),PvqR.std(axis=0)
#	PvHm,PvRm       = moving_average(PvHm,11),moving_average(PvRm,11)
#	PvHs,PvRs       = moving_average(PvHs,11),moving_average(PvRs,11)
#       PqHm,PqRm       = moving_average(PqHm,11),moving_average(PqRm,11)
#       PqHs,PqRs       = moving_average(PqHs,11),moving_average(PqRs,11)
#       PvqHm,PvqRm     = moving_average(PvqHm,11),moving_average(PvqRm,11)
#       PvqHs,PvqRs     = moving_average(PvqHs,11),moving_average(PvqRs,11)
#       lons = lons[5:360-5]

	"""
	pl.figure(1)
	pl.plot(lons,PvHm,'b')
	pl.fill_between(lons,PvHm-PvHs,PvHm+PvHs,color='b',alpha=0.1)
	pl.plot(lons,PvRm,'r')
	pl.fill_between(lons,PvRm-PvRs,PvRm+PvRs,color='r',alpha=0.1)
	pl.title('Mass: hist %s rcp85 %s' % (round(PvHm.mean(),5),round(PvRm.mean(),5)))
	pl.xlim(0,359)
        pl.figure(2)
        pl.plot(lons,PqHm,'b')
        pl.fill_between(lons,PqHm-PqHs,PqHm+PqHs,color='b',alpha=0.1)
        pl.plot(lons,PqRm,'r')
        pl.fill_between(lons,PqRm-PqRs,PqRm+PqRs,color='r',alpha=0.1)
        pl.title('Vapor:  %s rcp85 %s' % (round(PqHm.mean(),5),round(PqRm.mean(),5)))
        pl.xlim(0,359)
        pl.figure(3)
        pl.plot(lons,PvqHm,'b')
        pl.fill_between(lons,PvqHm-PvqHs,PvqHm+PvqHs,color='b',alpha=0.1)
        pl.plot(lons,PvqRm,'r')
        pl.fill_between(lons,PvqRm-PvqRs,PvqRm+PvqRs,color='r',alpha=0.1)
        pl.title('Moist: %s rcp85 %s' % (round(PvqHm.mean(),5),round(PvqRm.mean(),5)))
        pl.xlim(0,359)
	pl.show()
	"""

	pl.figure(1)
	pl.plot(lons,RHm,'b')
	pl.fill_between(lons,RHm-RHs,RHm+RHs,color='b',alpha=0.1)
        pl.plot(lons,RRm,'r')
        pl.fill_between(lons,RRm-RRs,RRm+RRs,color='r',alpha=0.1)
	pl.title('correlation: hist %s rcp85 %s' % (round(RHm.mean(),5),round(RRm.mean(),5)))
	pl.xlim(0,359)
	pl.xlabel('Longitude')
	pl.ylabel('Correlation coefficient')
        pl.figure(2)
        pl.plot(lons,SHm,'b')
        pl.fill_between(lons,SHm-SHs,SHm+SHs,color='b',alpha=0.1)
        pl.plot(lons,SRm,'r')
        pl.fill_between(lons,SRm-SRs,SRm+SRs,color='r',alpha=0.1)
	pl.title('regression')
	pl.xlim(0,359)
	pl.xlabel('Longitude')
	pl.ylabel('Regression coefficient')
	pl.show()

def pdf(norm=False):
	if norm == False: char = 'full'
	if norm == True:  char = 'norm'
        # Lambert projector and km per deg longitude at blats
        Re        = 6371
        proj      = LambertProjector(boundinglat=50,resolution=80.)
        xres      = 2*np.pi*np.cos(70*np.pi/180)*Re/360.
	sf        = 1e09/(24.*60*60*xres*1000)
	Ha1,Ha2   = [],[]
	Rh,Rr     = [],[]
	Ir,Ih     = [],[]
	QH,QR,QRT = [],[],[]
	VH,VR     = [],[]
	VQH,VQR   = [],[]
	# Loop through models
	Models    = [i[9:] for i in glob.glob('../rcp85/*')]
	lons      = range(0,360)
	lon1,lon2 = 0,359
	indices   = range(lon1,lon2+1)
	tvar      = []
	for Source in Models:
		print Source
		# Open historical flux files
		hdir        = '/qnap/cian/cmip/scripts'
		vqht,datesh = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.30-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))
		vht,datesh  = unpick('%s/newfluxfiles/%s/%s/%s.mass.1981-2005.30-1000hPa.%sN.DJF.p'  % (hdir,70,Source,Source,70))
		qht         = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1981-2005.30-1000hPa.%sN.DJF.p' % (hdir,70,Source,Source,70))[0]/(97000/9.8)
		vhft        = vqht/qht
		# Interpolate
#		print np.ma.masked_where(vqht<195,vqht).sum()/np.ma.masked_where(vqht<0,vqht).sum()
		vqht,vht,qht,vhft = sf*interp(vqht),sf*interp(vht)/1000,interp(qht),sf*interp(vhft)/1000
		vqht,vht,qht,vhft = vqht.take(indices,mode='wrap',axis=1),vht.take(indices,mode='wrap',axis=1),\
				    qht.take(indices,mode='wrap',axis=1),vhft.take(indices,mode='wrap',axis=1)
		vqht,vht,qht,vhft = vqht.reshape(-1),vht.reshape(-1),qht.reshape(-1),vhft.reshape(-1)
                # Open rcp85 flux files
		rdir        = '/mnt/climstorage/cian/scripts'
		vqrt,datesr = unpick('%s/newfluxfiles/%s/%s/%s.moist.2075-2100.30-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))
		vrt,datesr  = unpick('%s/newfluxfiles/%s/%s/%s.mass.2075-2100.30-1000hPa.%sN.DJF.p'  % (rdir,70,Source,Source,70))
		qrt         = unpick('%s/newfluxfiles/%s/%s/%s.vapor.2075-2100.30-1000hPa.%sN.DJF.p' % (rdir,70,Source,Source,70))[0]/(97000/9.8)
		vrft        = vqrt/qrt
		# Interpolate
		vqrt,vrt,qrt,vrft = sf*interp(vqrt),sf*interp(vrt)/1000,interp(qrt),sf*interp(vrft)/1000
		vqrt,vrt,qrt,vrft = vqrt.take(indices,mode='wrap',axis=1),vrt.take(indices,mode='wrap',axis=1),\
				    qrt.take(indices,mode='wrap',axis=1),vrft.take(indices,mode='wrap',axis=1)
		vqrt,vrt,qrt,vrft = vqrt.reshape(-1),vrt.reshape(-1),qrt.reshape(-1),vrft.reshape(-1)
		Xh,Yh,Zh,Xr,Yr,Zr = [],[],[],[],[],[]
		for i in range(len(vht)):
			if vht[i]>0:#-10e30:
				Xh.append(qht[i])
				Yh.append(vhft[i])
				Zh.append(vqht[i])
		for i in range(len(vrt)):
			if vrt[i]>0:#-10e30:
                                Xr.append(qrt[i])
                                Yr.append(vrft[i])
                                Zr.append(vqrt[i])
		qht,vht,vqht = np.array(Xh),np.array(Yh),np.array(Zh)
		qrt,vrt,vqrt = np.array(Xr),np.array(Yr),np.array(Zr)
		tvar.append(100.*(vrt.std()-vht.std())/vht.std())
#		qht,vht,vqht = (qht-qht.mean())/qht.std(),(vht-vht.mean())/vht.std(),(vqht-vqht.mean())/vqht.std()
#		qrt,vrt,vqrt = (qrt-qrt.mean())/qrt.std(),(vrt-vrt.mean())/vrt.std(),(vqrt-vqrt.mean())/vqrt.std()
		rat = qrt.mean()/qht.mean()
		Qh,edgesq   = np.histogram(qht,bins=np.linspace(0,0.007,75),range=None,normed=True,weights=None,density=None)
		Qr,edgesq   = np.histogram(qrt,bins=np.linspace(0,0.007,75),range=None,normed=True,weights=None,density=None)
		Qrt,edgesq  = np.histogram(rat*qht,bins=np.linspace(0,0.007,75),range=None,normed=True,weights=None,density=None)
		Vh,edgesv   = np.histogram(vht,bins=np.linspace(0,200,75),range=None,normed=True,weights=None,density=None)
		Vr,edgesv   = np.histogram(vrt,bins=np.linspace(0,200,75),range=None,normed=True,weights=None,density=None)
		VQh,edgesvq = np.histogram(vqht,bins=np.linspace(0,400,75),range=None,normed=True,weights=None,density=None)
		VQr,edgesvq = np.histogram(vqrt,bins=np.linspace(0,400,75),range=None,normed=True,weights=None,density=None)
#		Qh,edgesq   = np.histogram(qht,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		Qr,edgesq   = np.histogram(qrt,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		Qrt,edgesq  = np.histogram(rat*qht,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		Vh,edgesv   = np.histogram(vht,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		Vr,edgesv   = np.histogram(vrt,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		VQh,edgesvq = np.histogram(vqht,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
#		VQr,edgesvq = np.histogram(vqrt,bins=np.linspace(-4,4,75),range=None,normed=True,weights=None,density=None)
		QH.append(Qh)
		QR.append(Qr)
		QRT.append(Qrt)
                VH.append(Vh)
                VR.append(Vr)
                VQH.append(VQh)
                VQR.append(VQr)

		# Standardise data
#		qht,vht,vqht = (qht-qht.mean())/qht.std(),(vht-vht.mean())/vht.std(),(vqht-vqht.mean())/vqht.std()
#		qrt,vrt,vqrt = (qrt-qrt.mean())/qrt.std(),(vrt-vrt.mean())/vrt.std(),(vqrt-vqrt.mean())/vqrt.std()
		# Correlation coefficients
		slopeh,intercepth,rh,p_valueh,std_errh = stats.linregress(vht,qht)
		sloper,interceptr,rr,p_valuer,std_errr = stats.linregress(vrt,qrt)
		Rh.append(slopeh)
		Rr.append(sloper)
                Ih.append(intercepth)
                Ir.append(interceptr)
		# Bivariate histogram
#		yedges,xedges    = np.arange(-1.2e05,1.2e05+5000,5000),np.linspace(0,0.004,31)
		yedges,xedges    = np.linspace(0,220,50),np.linspace(0,0.002,25)
#		yedges,xedges    = 25,25
		Hh,xedges,yedges = np.histogram2d(qht,vht,bins=[xedges,yedges],range=None,normed=norm,weights=None)
		Hr,xedges,yedges = np.histogram2d(qrt,vrt,bins=[xedges,yedges],range=None,normed=norm,weights=None)
		xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
		if norm == False:
			Hh = Hh*xedges[:,np.newaxis]
			Hh = Hh*yedges[np.newaxis,:]/len(datesh)
			Hr = Hr*xedges[:,np.newaxis]
			Hr = Hr*yedges[np.newaxis,:]/len(datesr)
		Ha1.append(Hh)
		Ha2.append(Hr-Hh)

	"""
	print np.array(tvar).mean()
	fig,ax1   = pl.subplots(num=1)
	ax2       = ax1.twinx()
	QH,QR,QRT = np.array(QH).mean(axis=0),np.array(QR).mean(axis=0),np.array(QRT).mean(axis=0)
	ax1.plot(edgesq[0:-1],QH,'b',linewidth=2)
	ax1.plot(edgesq[0:-1],QR,'r',linewidth=2)
#	ax1.plot(edgesq[0:-1],QRT,'r--',linewidth=1)
	ax2.plot(edgesq[0:-1],QR/QH,'k--',linewidth=1)
	ax1.set_xlim(edgesq[0],edgesq[-2])
	ax2.set_ylim(0,500)
	ax1.set_xlabel('Standard deviation')
	pl.title('Humidity: %s-%sE' % (lon1,lon2))
	pl.savefig('figs/pdf/q.pdf.%s-%sE.pdf' % (lon1,lon2))

	fig,ax1 = pl.subplots(num=2)
	ax2     = ax1.twinx()
	VH,VR   = np.array(VH).mean(axis=0),np.array(VR).mean(axis=0)
	ax1.plot(edgesv[0:-1],VH,'b',linewidth=2)
	ax1.plot(edgesv[0:-1],VR,'r',linewidth=2)
	ax2.plot(edgesv[0:-1],VR/VH,'k--',linewidth=1)
	pl.xlim(edgesv[0],edgesv[-2])
	ax2.set_ylim(0,2)
	ax1.set_xlabel('Standard deviation')
	pl.title('Mass flux: %s-%sE' % (lon1,lon2))
	pl.savefig('figs/pdf/v.pdf.%s-%sE.pdf' % (lon1,lon2))

	fig,ax1 = pl.subplots(num=3)
	ax2     = ax1.twinx()
	VQH,VQR = np.array(VQH).mean(axis=0),np.array(VQR).mean(axis=0)
	ax1.plot(edgesvq[0:-1],VQH,'b',linewidth=2)
	ax1.plot(edgesvq[0:-1],VQR,'r',linewidth=2)
	ax2.plot(edgesvq[0:-1],VQR/VQH,'k--',linewidth=1)
	pl.xlim(edgesvq[0],edgesvq[-2])
	ax2.set_ylim(0,200)
	ax1.set_xlabel('Standard deviation')
	pl.title('Moisture flux: %s-%sE' % (lon1,lon2))
	pl.savefig('figs/pdf/vq.pdf.%s-%sE.pdf' % (lon1,lon2))
	pl.show()
	"""

	Ha1 = np.array(Ha1).mean(axis=0)*1000
	Ha2 = np.array(Ha2).mean(axis=0)*1000
        Rh  = np.array(Rh).mean(axis=0)
        Rr  = np.array(Rr).mean(axis=0)
        Ih  = np.array(Ih).mean(axis=0)
        Ir  = np.array(Ir).mean(axis=0)

	#cseq1 = np.arange(0,12+1,1)
	#cseq2 = np.arange(-8,8+1,1)
	cseq1,cseq2 = 15,15
	pl.figure(1)
	cmap   = pl.cm.coolwarm
	colors = cmap(np.linspace(0.5, 1, cmap.N/2))
	cmap   = LinearSegmentedColormap.from_list('Upper Half', colors)
	cf     = pl.contourf(yedges,xedges,Ha1,cseq1,cmap=cmap,extend='max')
	cbar   = pl.colorbar(cf)
	cbar.set_label('Moisture flux [kg m$^{-1} s^{-1}$]')
	pl.xlabel('Mass flux [Mg m$^{-1}$ s$^{-1}$]')
	pl.ylabel('Specific humidity [kg kg$^{-1}$]')
#	pl.plot(yedges,[Rh*yedges[i] + Ih for i in range(len(yedges))],'k--')
	pl.plot(np.linspace(yedges[0],yedges[-1],100),0.073/np.linspace(yedges[0],yedges[-1],100),'k')
	pl.xlim(yedges[0],yedges[-1])
	pl.ylim(xedges[0],xedges[-1])
	pl.title(Rh)
	pl.savefig('figs/bivar.hist.%s.%s-%sE.pdf' % (char,lons[lon1],lons[lon2]))

	pl.figure(2)
	cf   = pl.contourf(yedges,xedges,Ha2,cseq2,cmap=pl.cm.coolwarm,extend='both')
	cl   = pl.contour(yedges,xedges,Ha2,levels=[0],colors='k',linewidths=2)
	cbar = pl.colorbar(cf)
	cbar.set_label('Moisture flux [kg m$^{-1} s^{-1}$]')
	pl.xlabel('Mass flux [Mg m$^{-1}$ s$^{-1}$]')
	pl.ylabel('Specific humidity [kg kg$^{-1}$]')
#	pl.plot(yedges,[Rr*yedges[i] + Ir for i in range(len(yedges))],'k--')
	pl.plot(np.linspace(yedges[0],yedges[-1],100),0.073/np.linspace(yedges[0],yedges[-1],100),'k')
	pl.xlim(yedges[0],yedges[-1])
	pl.ylim(xedges[0],xedges[-1])
	pl.title(Rr)
	pl.savefig('figs/bivar.%s.%s-%sE.pdf' % (char,lons[lon1],lons[lon2]))
	pl.show()


def reg():

	# Open historical flux files
	hdir        = '/qnap/cian/cmip/scripts'
	blat        = 70
	Models      = ['ERAInt'] + [g[27:] for g in glob.glob('/qnap/cian/cmip/historical/*')]
	S = []
	for Source in Models:
		vqft,datesf = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.30-1000hPa.%sN.DJF.p' % (hdir,blat,Source,Source,blat))
		vqst,datess = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,blat,Source,Source,blat))
		# Interpolate
		vqft,vqst = interp(vqft).reshape(-1),interp(vqst).reshape(-1)
		# Regress
		slope,intercept,r,p_value,std_err = stats.linregress(vqft,vqst)
		S.append(slope)
		# Plot
#		xs   = np.linspace(-1500,2000,20)
#		line = [x*slope + intercept for x in xs]
#		pl.plot(vqft,vqst,'k+')
#		pl.plot(xs,line,'k--')
#		pl.title('%s: y = %sx + %s' % (Source,slope,r))
#		pl.show()
	S = np.array(S)
	print S[0],S[1:].mean(),S[1:].std()

def split(blats,case='+ve',thresh1=0,thresh2=100):

	# Lambert projector and km per deg longitude at blats
	Re   = 6371
	proj = LambertProjector(boundinglat=50,resolution=80.)
	xres = 2*np.pi*np.cos(blats*np.pi/180)*Re/360.

	F1,F2,F3 = [],[],[]
	F4,F5,F6 = [],[],[]
	for blat in blats:
		# Loop through models
		Models = [i[9:] for i in glob.glob('../rcp85/*')]
		VQ,V,Q = [],[],[]	
		for Source in Models:
			print Source
			# Open historical flux files
			hdir        = '/qnap/cian/cmip/scripts'
			vqht,datesh = unpick('%s/newfluxfiles/%s/%s/%s.moist.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,blat,Source,Source,blat))
			vht,datesh  = unpick('%s/newfluxfiles/%s/%s/%s.mass.1981-2005.600-1000hPa.%sN.DJF.p'  % (hdir,blat,Source,Source,blat))
			qht         = unpick('%s/newfluxfiles/%s/%s/%s.vapor.1981-2005.600-1000hPa.%sN.DJF.p' % (hdir,blat,Source,Source,blat))[0]/(40000/9.8)
			# Interpolate
			vqht,vht,qht   = interp(vqht),interp(vht),interp(qht)

			# Open rcp85 flux files
			rdir        = '/mnt/climstorage/cian/scripts'
			vqrt,datesr = unpick('%s/newfluxfiles/%s/%s/%s.moist.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,blat,Source,Source,blat))
			vrt,datesr  = unpick('%s/newfluxfiles/%s/%s/%s.mass.2075-2100.600-1000hPa.%sN.DJF.p'  % (rdir,blat,Source,Source,blat))
			qrt         = unpick('%s/newfluxfiles/%s/%s/%s.vapor.2075-2100.600-1000hPa.%sN.DJF.p' % (rdir,blat,Source,Source,blat))[0]/(40000/9.8)
			# Interpolate
			vqrt,vrt,qrt   = interp(vqrt),interp(vrt),interp(qrt)
			# Northward clim
			if case == '+ve':
				# Historical
				vqh  = np.ma.masked_where(vqht<thresh1,vqht)
				vqh  = np.ma.masked_where(vqht>thresh2,vqh).sum(axis=0)/len(datesh)
				vh   = np.ma.masked_where(vqht<thresh1,vht)
				vh   = np.ma.masked_where(vqht>thresh2,vh).sum(axis=0)/len(datesh)
				qh   = np.ma.masked_where(vqht<thresh1,qht)
				qh   = np.ma.masked_where(vqht>thresh2,qh).mean(axis=0)
				# RCP85
				vqr  = np.ma.masked_where(vqrt<thresh1,vqrt)
				vqr  = np.ma.masked_where(vqrt>thresh2,vqr).sum(axis=0)/len(datesr)
				vr   = np.ma.masked_where(vqrt<thresh1,vrt)
				vr   = np.ma.masked_where(vqrt>thresh2,vr).sum(axis=0)/len(datesr)
				qr   = np.ma.masked_where(vqrt<thresh1,qrt)
				qr   = np.ma.masked_where(vqrt>thresh2,qr).mean(axis=0)
			vq = (vqr.data - vqh.data)
			v  = qh*(vr.data - vh.data)
			q  = vh*(qr.data - qh.data)

			"""
			pl.plot(vqh,'DarkMagenta',linewidth=1,alpha=1)
			pl.plot(1.0*vh*qh,'r',linewidth=1,alpha=1)
			pl.title('%s historical' % (Source))
			pl.plot(vqr,'DarkMagenta',linestyle='dashed',linewidth=1,alpha=1)
			pl.plot(1.0*vr*qr,'r',linestyle='dashed',linewidth=1,alpha=1)
			pl.title('%s rcp85' % (Source))
			pl.show()
			"""

			VQ.append(vq)
			V.append(v)
			Q.append(q)
		VQ,V,Q = np.array(VQ),np.array(V),np.array(Q)
		# CMIP5 model mean bias
		vm,qm,vqm   = V.mean(axis=0),Q.mean(axis=0),VQ.mean(axis=0)
		vs,qs,vqs   = V.std(axis=0),Q.std(axis=0),VQ.std(axis=0)
		vqmax,vqmin = V.max(axis=0),V.min(axis=0)
		F1.append(vqm)
		F2.append(vm)
		F3.append(qm)
                F4.append(vqs)
                F5.append(vs)
                F6.append(qs)
	F1,F2,F3 = np.array(F1),np.array(F2),np.array(F3)
	F4,F5,F6 = np.array(F4),np.array(F5),np.array(F6)
	F1,F2,F3 = F1.squeeze(),F2.squeeze(),F3.squeeze()
	F4,F5,F6 = F4.squeeze(),F5.squeeze(),F6.squeeze()

	# Shift grid
	lons = np.arange(0,360+1,1)
	l0   = np.argmin((lons-(360-80))**2)
	vm   = np.array(list(F2[l0+1:]) + list(F2[:l0+1]))
	qm   = np.array(list(F3[l0+1:]) + list(F3[:l0+1]))
	vqm  = np.array(list(F1[l0+1:]) + list(F1[:l0+1]))
        vs   = np.array(list(F5[l0+1:]) + list(F5[:l0+1]))
        qs   = np.array(list(F6[l0+1:]) + list(F6[:l0+1]))
        vqs  = np.array(list(F4[l0+1:]) + list(F4[:l0+1]))
	lons = np.arange(-179,180+1,1)
	# Smooth
#	vm,qm,vqm   = moving_sum(vm,11),moving_sum(qm,11),moving_sum(vqm,11)
#	vs,qs,vqs   = moving_sum(vs,11),moving_sum(qs,11),moving_sum(vqs,11)
#	vqmin,vqmax = moving_sum(vqmin,11),moving_sum(vqmax,11)
	vm,qm,vqm   = moving_average(vm,11),moving_average(qm,11),moving_average(vqm,11)
	vs,qs,vqs   = moving_average(vs,11),moving_average(qs,11),moving_average(vqs,11)
	vqmin,vqmax = moving_average(vqmin,11),moving_average(vqmax,11)
	lons = lons[5:360-5]
	return lons,vqm,vqs,vm,qm

def plot(lons,vqm,vqs,vm,qm):

	pl.figure(1)
	PP,MM = [],[]
	pp, = pl.plot(lons,qm,'b',linewidth=1.5,alpha=1)
	PP.append(pp)
	MM.append(r'$\bar{q}^{\prime}\bar{v}^{R}$')

	pp, = pl.plot(lons,vm,'r',linewidth=1.5,alpha=1)
	PP.append(pp)
	MM.append(r'$\bar{v}^{\prime}\bar{q}^{R}$')

	pp, = pl.plot(lons,vqm,'DarkMagenta',linewidth=1.5,alpha=1)
	PP.append(pp)
	MM.append(r'$F^{\prime}$')

	pp, = pl.plot(lons,vqm-(vm+qm),'k',linewidth=1,alpha=0.5)
	PP.append(pp)
	MM.append(r'$\epsilon$')
	pl.fill_between(lons,vqm-vqs,vqm+vqs,color='DarkMagenta',alpha=0.1)

	pl.plot([lons[0],lons[-1]],[0,0],'k--')
	pl.xlim(-179,180)
	pl.ylim(-10,50)
	pl.xlabel('Longitude')
	pl.ylabel('Moisture transport bias [Tg day$^{-1}$ deg$^{-1}$]')
	pl.xticks([-160,-100,-40,20,80,140],['-60','0','60','120','180','-120'])
#	pl.yticks([-50,-40,-30,-20,-10,0,10,20,30],[-50,-40,-30,-20,-10,0,10,20,30])
	lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
	pl.setp(lg.get_title(),fontsize=10)
	pl.savefig('figs/testa.pdf',format='pdf')

#reg()
pdf(norm=True)
#corr()

#blats   = np.arange(70,70+1,1)
#case    = str(sys.argv[1])
#thresh1 = int(sys.argv[2])
#thresh2 = int(sys.argv[3])
#lons,vqm,vqs,vm,qm = split(blats=blats,case=case,thresh1=thresh1,thresh2=thresh2)
#plot(lons,vqm,vqs,vm,qm)

"""
lons,vqm1,vqs1,vm1,qm1 = split(blats=np.arange(70,70+1,1),case='+ve',thresh1=0,thresh2=240)
lons,vqm2,vqs2,vm2,qm2 = split(blats=np.arange(70,70+1,1),case='+ve',thresh1=240,thresh2=3000)

vqm,vm,qm = vqm1+vqm2,vm1+vm2,qm1+qm2
plot(lons,vqm,vqs1,vm,qm)
"""
