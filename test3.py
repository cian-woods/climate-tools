from UnPickle import *

import numpy as np
import matplotlib.pyplot as pl

def sortflux(flux,dates):
	E = {}
	for i in range(len(flux)):
		year,month,day,hour = dates[i]
		if month == 12:
			year = year + 1
		try:
			E[year].append(np.ma.masked_where(flux[i]<0,flux[i]).sum())
		except:
			E[year] = []
			E[year].append(np.ma.masked_where(flux[i]<0,flux[i]).sum())
	return E

G,Q,D = unpick('298events.p')
timei = sum([len(D[i]) for i in range(len(D))])
q,d1  = [ np.array([np.sum(Q[i][j]) for j in range(len(Q[i]))]) for i in range(len(Q))],[D[i][0] for i in range(len(D))]
p,d2  = unpick('/qnap/cian/cmip/scripts/fluxfiles/ERAInt.0.1000.1980.2012.DJF.p')
timea = d2.index((2010,2,28,18))-d2.index((1989,12,1,0))
print 81987./(timea*360)

p     = np.array(p)
print np.ma.masked_where(p.reshape(-1)<0,p.reshape(-1)).mean()

E1    = sortflux(q,d1)
E2    = sortflux(p,d2)

years = range(1990,2010+1,1)
E1 = [E1[year] for year in years]
E2 = [E2[year] for year in years]
E1 = np.array([np.sum(E1[i]) for i in range(len(E1))])
E2 = np.array([np.sum(E2[i]) for i in range(len(E2))])
R  = E2-E1
R  = (R-R.mean())/R.std()

pl.plot(years,E1,'k',linewidth=2)
pl.plot(years,E2,'r',linewidth=2)
pl.ylim(-2.2,2.8)
pl.xlim(1990,2010)
pl.show()
pl.savefig('figs/resid.pdf',format='pdf')

