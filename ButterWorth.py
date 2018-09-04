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
