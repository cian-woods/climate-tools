from matplotlib import pyplot as mp
import numpy as np

def gaussian(x, mu, sig):
	return (1/(sig*np.sqrt(2*np.pi)))*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def gaussian2(x, amp, cen, wid):
	return amp * exp(-(x-cen)**2 /wid)

N1,N2     = 200,6
mu1,mu2   = 100,110
sig1,sig2 = 15,15
cut       = 125

xs = np.linspace(-5*sig1+mu1,5*sig1+mu1,350)
g1 = N1*gaussian(xs, mu1, sig1)
g2 = N2*gaussian(xs, mu2, sig2)

ixs  = np.argmin((xs-cut)**2)
int1 = g1[ixs:].sum()
int2 = g2[ixs:].sum()

PP,MM = [],[]
fig,ax1 = mp.subplots(num=1)
ax2     = ax1.twinx()
ax2.plot(xs,g2/g1,'k')
ax2.plot(xs[ixs],(g2/g1)[ixs],'ko')
ax2.set_yscale('log')
pp, = ax1.plot(xs,g1,'b')
PP.append(pp)
MM.append('$\mu = %s; s = %s$' % (mu1,sig1))
pp, = ax1.plot(xs,g2,'r')
PP.append(pp)
MM.append('$\mu = %s; s = %s$' % (mu2,sig2))
ax1.set_ylabel('Frequency')
ax2.set_ylabel('Ratio [red:blue]')
mp.xlim(xs[0],xs[-1])
lg = mp.legend(PP,MM,loc=2,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
mp.title('%s' % round(int2/int1,3))
mp.savefig('test.png',format='png')
mp.show()
