from scipy.signal import butter, lfilter
from stipling import stipling

import glob

def butter_bandpass(lowcut, highcut, fs, order=5):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	b, a = butter(order, [low, high], btype='band')
	return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
	b, a = butter_bandpass(lowcut, highcut, fs, order=order)
	y = lfilter(b, a, data, axis=0)
	return y

def all(lwl,hwl):
	# Lambert
	proj = LambertProjector(boundinglat=40,resolution=80.,Projection='nplaea')
	# CMIP5
	syn1,syn2 = [],[]
	Models = [g[9:] for g in glob.glob('../rcp85/*') if g[9:]!='FGOALS-g2']
	for Model in Models:
		#S1,lons1,lats1 = unpick('synopfiles/%s.msl.%s-%sdays-pass.std.p' % (Model,lwl,hwl))
		S1,lons1,lats1 = unpick('/qnap/cian/cmip/scripts/synopfiles/%s.msl.%s-%sdays-pass.std.rcp85.p' % (Model,lwl,hwl))
		S1 = proj(S1,lons1,lats1)
		S2,lons2,lats2 = unpick('/qnap/cian/cmip/scripts/synopfiles/%s.msl.%s-%sdays-pass.std.historical.p' % (Model,lwl,hwl))
		S2 = proj(S2,lons2,lats2)
		# Plot std
		cseq = np.arange(-200,200+25,25)
		plt.figure(1)
		cf   = proj.m.contourf(proj.x,proj.y,S1-S2,cseq,cmap=plt.cm.RdBu_r,extend='both')
		cbar = plt.colorbar(cf)
		cbar.set_label('Standard deviation [hPa]')
		proj.m.drawcoastlines(color='0.4')
		proj.m.drawparallels([70,80],latmax=90)
		plt.savefig('figs/bias/synop/%s.%s-%sdays.pdf' % (Model,lwl,hwl),format='pdf')
		plt.close()
		syn1.append(S1)
		syn2.append(S2)
	syn1 = np.array(syn1)/10.
	syn2 = np.array(syn2)/10.
	syn  = syn1 - syn2
#	syn  = syn1/syn1.max() - syn2/syn2.max()
	stipx,stipy = stipling(syn,xx=proj.x,yy=proj.y,x=None,y=None,thresh=0.8)
	syn = syn.mean(axis=0)

	# Plot std
	cseq = 17
	plt.figure(1)
	cf   = plt.contourf(proj.x,proj.y,syn,cseq,cmap=plt.cm.RdBu_r,extend='both')
	plt.plot(stipx[::3,::3],stipy[::3,::3],'k.',alpha=0.5)
	cbar = plt.colorbar(cf)
	cbar.set_label('Standard deviation [10$^{-1}$ hPa]')
	proj.m.drawcoastlines(color='0.4')
	proj.m.drawparallels([70,80],latmax=90)
	plt.savefig('figs/bias/synop/%s.%s-%sdays.pdf' % ('all',lwl,hwl),format='pdf')
	plt.close()

        plt.figure(1)
        cf   = plt.contourf(proj.x,proj.y,syn1.mean(axis=0),cseq,cmap=plt.cm.RdBu_r,extend='both')
        cbar = plt.colorbar(cf)
        cbar.set_label('Standard deviation [10$^{-1}$ hPa]')
        proj.m.drawcoastlines(color='0.4')
        proj.m.drawparallels([70,80],latmax=90)
        plt.savefig('figs/bias/synop/%s.%s-%sdays.rcp85.pdf' % ('all',lwl,hwl),format='pdf')
        plt.close()

        plt.figure(1)
        cf   = plt.contourf(proj.x,proj.y,syn2.mean(axis=0),cseq,cmap=plt.cm.RdBu_r,extend='both')
        cbar = plt.colorbar(cf)
        cbar.set_label('Standard deviation [10$^{-1}$ hPa]')
        proj.m.drawcoastlines(color='0.4')
        proj.m.drawparallels([70,80],latmax=90)
        plt.savefig('figs/bias/synop/%s.%s-%sdays.hist.pdf' % ('all',lwl,hwl),format='pdf')
        plt.close()

def makeFile(Source,lwl,hwl,Season):
        # Make file
        years = range(1980,2012+1,1)
	fname = 'synopfiles/%s.msl.%s-%sdays-pass.std.%s-%s.%s.p' % (Source,lwl,hwl,years[0],years[-1],Season)
	print fname
        if os.path.isfile(fname) == False:
                # MSL DataServer CMIP
                ds = reDataServer(Field='slp',Source=Source,LevType='surface_analysis')
                S1 = []
                for year in years:
                        x1 = ds.getDataSnaps(year,Season=Season)
			print year,x1.shape
                        x1 = x1-x1.mean(axis=0)
                        if (lwl!=0) and (hwl!=0):
                                # Sample rate and desired cutoff frequencies (days per DJF season).
                                fs1 = len(x1)
                                lowcut,highcut = 1.*fs1/(4*hwl),1.*fs1/(4*lwl)
                                # Apply filter
                                y1 = butter_bandpass_filter(x1, lowcut, highcut, fs1, order=6)
                        else:
                                y1 = x1
                        S1.append(y1.std(axis=0))
                S1 = np.array(S1)
                lons1,lats1 = ds.lon,ds.lat
                toPick([S1,lons1,lats1],fname)
        else:
                S1,lons1,lats1 = unpick(fname)

if __name__ == "__main__":
	import numpy as np
	import matplotlib.pyplot as plt
	import os,sys
	from ReanalysisDataServer import DataServer as reDataServer
	from cmipDataServer import DataServer as cmipDataServer
	from LambertProjector import *
	from toPick import *
	from UnPickle import *

	Source  = str(sys.argv[1])
	lwl,hwl = int(sys.argv[2]),int(sys.argv[3])	

	all(lwl,hwl)
	sys.exit()
 
	# Make file
	years = range(1901,1919+1,1)
	if os.path.isfile(fname) == False:
		# MSL DataServer CMIP
		ds = cmipDataServer(Field='psl',Source=Source,ExpType='rcp85',LevType='surface',DataFreq='day')
		S1 = []
		for year in years:
			print year
			x1 = ds.getDataSnaps(year,Season='DJF',step=24.)
			x1 = x1-x1.mean(axis=0)
			if (lwl!=0) and (hwl!=0):
				# Sample rate and desired cutoff frequencies (days per DJF season).
				fs1 = len(x1)
				lowcut,highcut = 1.*fs1/hwl,1.*fs1/lwl
				# Apply filter
				y1 = butter_bandpass_filter(x1, lowcut, highcut, fs1, order=6)
			else:
				y1 = x1
			S1.append(y1.std(axis=0))
		S1 = np.array(S1).mean(axis=0)
		lons1,lats1 = ds.lon,ds.lat
		toPick([S1,lons1,lats1],fname)
	else:
		S1,lons1,lats1 = unpick(fname)

	"""
	proj = LambertProjector(boundinglat=30,resolution=80.)
	S1,lons1,lats1 = unpick('synopfiles/%s.msl.2-6days-pass.std.p' % (Source))
	S1   = proj(S1,lons1,lats1)

	cseq = 20
	cf   = plt.contourf(proj.x,proj.y,S1,cseq,cmap=plt.cm.RdBu_r,extend='both')
	cbar = plt.colorbar(cf)
	cbar.set_label('Standard deviation [hPa]')
	proj.m.drawcoastlines(color='0.4')
	proj.m.drawparallels([70,80],latmax=90)
	plt.show()
	"""
