from UnPickle import *
from scipy import interpolate

import matplotlib.pyplot as pl
import numpy as np
import glob

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

thresh = 0
Dir    = '/qnap/cian/cmip/scripts/newfluxfiles/70'

# ERAInt
Model = 'ERAInt'
# Open file
vqERAt,datesERA = unpick('%s/%s/%s.moist.1981-2005.600-1000hPa.70N.DJF.p' % (Dir,Model,Model))
vERAt,datesERA  = unpick('%s/%s/%s.mass.1981-2005.600-1000hPa.70N.DJF.p'  % (Dir,Model,Model))
qERAt,datesERA  = unpick('%s/%s/%s.vapor.1981-2005.600-1000hPa.70N.DJF.p' % (Dir,Model,Model))
# Interpolate
vqERAt,vERAt,qERAt = interp(vqERAt),interp(vERAt),interp(qERAt)/(40000/9.8)
# Mask
vqERA = np.ma.masked_where(vqERAt<thresh,vqERAt).sum(axis=0)/len(datesERA)
vERA  = np.ma.masked_where(vqERAt<thresh,vERAt).sum(axis=0)/len(datesERA)
qERA  = np.ma.masked_where(vqERAt<thresh,qERAt).mean(axis=0)

# CMIP5
VQ,VQ_,V,Q = [],[],[],[]
Models = [g[27:] for g in glob.glob('/qnap/cian/cmip/historical/*')]
for Model in Models:
	print Model
	# Open files
	vqt,dates = unpick('%s/%s/%s.moist.1981-2005.600-1000hPa.70N.DJF.p' % (Dir,Model,Model))
	vt,dates  = unpick('%s/%s/%s.mass.1981-2005.600-1000hPa.70N.DJF.p'  % (Dir,Model,Model))
	qt,dates  = unpick('%s/%s/%s.vapor.1981-2005.600-1000hPa.70N.DJF.p' % (Dir,Model,Model))
	# Interpolate
	vqt,vt,qt = interp(vqt),interp(vt),interp(qt)/(40000/9.8)
	# Mask
	vq = np.ma.masked_where(vqt<thresh,vqt).sum(axis=0)/len(dates)
	v  = np.ma.masked_where(vqt<thresh,vt).sum(axis=0)/len(dates)
	q  = np.ma.masked_where(vqt<thresh,qt).mean(axis=0)
	VQ.append(vq)
	VQ_.append(v*q)
	V.append(v)
	Q.append(q)
vqm,vqm_,vm,qm = np.array(VQ).mean(axis=0),np.array(VQ_).mean(axis=0),np.array(V).mean(axis=0),np.array(Q).mean(axis=0)
# Plot
pl.figure(1)
pl.plot(vqm,'DarkMagenta')
pl.plot(vqm_,'DarkMagenta',linestyle='dashed')
pl.show()
