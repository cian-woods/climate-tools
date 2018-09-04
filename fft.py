from cmipDataServer import *
from scipy import interpolate,fftpack
from scipy.fftpack import fft, rfft, fftfreq

import numpy as np
import matplotlib.pyplot as pl
import sys
import glob

def interpolateND(field,xold,xnew,axis):
	# field in (time)x(lev)x(lon)
	f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=True)
	field = f(xnew)
	return field

years1 = range(1981,2005+1,1)
years2 = range(2075,2099+1,1)
lats   = range(68,72+1,1)
lons   = range(0,10+1,1)

SS1,F1 = [],[]
SS2,F2 = [],[]
Models = [g[9:] for g in glob.glob('../rcp85/*') if (g[9:]!='GFDL-CM3') and (g[9:]!='GFDL-ESM2G') and (g[9:]!='GFDL-ESM2M')]
for Model in Models:
	ds1    = DataServer(Field='psl',Source=Model,ExpType='historical',LevType='surface',DataFreq='day')
	ds2    = DataServer(Field='psl',Source=Model,ExpType='rcp85',LevType='surface',DataFreq='day')
	dlons1 = np.append(ds1.lon,[360+ds1.lon[0]],axis=0)
	dlons2 = np.append(ds2.lon,[360+ds2.lon[0]],axis=0)
	S1     = np.zeros((0,len(lats),len(lons)))
	A1,B1  = [],[]
        for year in years1:
                s1    = ds1.getDataSnaps(Year=year,Season='ONDJ')
                send1 = s1[:,:,0]
                s1    = np.append(s1,send1[:,:,np.newaxis],axis=2)
                s1    = interpolateND(s1,ds1.lat,lats,axis=1)
                s1    = interpolateND(s1,dlons1,lons,axis=2).squeeze()
                S1    = np.append(S1,s1,axis=0)

                s0        = s1.mean(axis=1).mean(axis=1)
                x         = range(len(s0))
                sl,inc    = np.polyfit(x,s0,1)
                line      = np.array([sl*xi + inc for xi in x])
                s0        = s0 - line

                W2        = fftfreq(len(s0), d=86400)
                n0        = np.argmin(W2)-1
                f_signal1 = fft(s0)
                A1.append(np.abs(f_signal1))
                B1.append(s0)
	S2  = np.zeros((0,len(lats),len(lons)))
	A2,B2 = [],[]
	for year in years2:
		s2    = ds2.getDataSnaps(Year=year,Season='ONDJ')
		send2 = s2[:,:,0]
		s2    = np.append(s2,send2[:,:,np.newaxis],axis=2)
		s2    = interpolateND(s2,ds2.lat,lats,axis=1)
		s2    = interpolateND(s2,dlons2,lons,axis=2).squeeze()
		S2    = np.append(S2,s2,axis=0)

		s0        = s2.mean(axis=1).mean(axis=1)
		x         = range(len(s0))
		sl,inc    = np.polyfit(x,s0,1)
		line      = np.array([sl*xi + inc for xi in x])
		s0        = s0 - line

		W2        = fftfreq(len(s0), d=86400)
		n0        = np.argmin(W2)-1
		f_signal2 = fft(s0)
		A2.append(np.abs(f_signal2))
		B2.append(s0)
	A1,A2 = np.array(A1).mean(axis=0),np.array(A2).mean(axis=0)
	B1,B2 = np.array(B1).mean(axis=0),np.array(B2).mean(axis=0)
	SS2.append(B2)
	SS1.append(B1)
	F2.append(A2)
	F1.append(A1)
SS1,F1 = np.array(SS1).mean(axis=0),np.array(F1).mean(axis=0)
SS2,F2 = np.array(SS2).mean(axis=0),np.array(F2).mean(axis=0)
x      = range(len(SS2))

pl.subplot(2,1,1)
pl.plot(x,SS2,'k')
pl.subplot(2,1,2)
pl.plot(W2[1:n0],F2[1:n0]**2,'r')
pl.plot(W2[1:n0],F1[1:n0]**2,'b')
pl.show()

"""
	S1 = S1[0:9125,:,:].mean(axis=1).mean(axis=1)
	S1 = S1 - S1.mean()
        S2 = S2[0:9125,:,:].mean(axis=1).mean(axis=1)
        S2 = S2 - S2.mean()
	N1 = len(S1)
	N2 = len(S2)
	print N1,N2
	W1 = fftfreq(N1, d=86400)
	W2 = fftfreq(N2, d=86400)
	f_signal1 = fft(S1)
	f_signal2 = fft(S2)
	SS1.append(S1)
	SS2.append(S2)
	F1.append(np.abs(f_signal1))
	F2.append(np.abs(f_signal2))

	pl.subplot(2,1,1)
	pl.plot(range(N2),SS2[-1])
	pl.subplot(2,1,2)
	pl.plot(W2,F2[-1])
	pl.xlabel('Frequency (Hz)')
	pl.xscale('log')
	pl.show()
"""

SS1,F1 = np.array(SS1).mean(axis=0),np.array(F1).mean(axis=0)
SS2,F2 = np.array(SS2).mean(axis=0),np.array(F2).mean(axis=0)
x      = range(len(SS1))

pl.subplot(2,1,1)
pl.plot(x,SS2-SS1)
pl.subplot(2,1,2)
pl.plot(W2,F2-F1)
pl.xlabel('Frequency (Hz)')
pl.xscale('log')
pl.show()
