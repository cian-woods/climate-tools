from netCDF4 import Dataset
from stipling import stipling
from scipy import interpolate
from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from UnPickle import *
from scipy import stats

import numpy as np
import matplotlib.pyplot as pl
import glob
import sys

def interp2d(field,x,y,n,kind='cubic'):
        # Assumes axes are increasing
        nx,ny   = int(len(x)*n),int(len(y)*n)
        new_x   = np.linspace(x[0],x[-1],nx)
        new_y   = np.linspace(y[0],y[-1],ny)
        # 2d interp function
        f = interpolate.interp2d(x,y,field,kind=kind,bounds_error=True)
        # Interpolate to new grid
        field = f(new_x,new_y)
        return field,new_x,new_y

def bivarPDF(x,y,xedges,yedges,norm=True):
        H,xedges,yedges = np.histogram2d(x,y,bins=[xedges,yedges],range=None,normed=norm,weights=None)
        xedges,yedges    = 0.5*(xedges[1:]+xedges[0:-1]),0.5*(yedges[1:]+yedges[0:-1])
        return xedges,yedges,H

def bivar2d(x,y,z,xedges,yedges):
        H,N = np.zeros((len(xedges),len(yedges))),np.zeros((len(xedges),len(yedges)))
        for i in range(len(x)):
                ix,iy = np.argmin((xedges-x[i])**2),np.argmin((yedges-y[i])**2)
                H[iy,ix] = H[iy,ix] + z[i]
                N[iy,ix] = N[iy,ix] + 1
        H,N = np.ma.masked_where(N==0,H),np.ma.masked_where(N==0,N)
        return xedges,yedges,H/N,N/N.sum(),N

def decomp(q_Amp,v_Amp,q_Phi,v_Phi):
	r        = np.cos(q_Phi-v_Phi)
	qm,vm,rm = q_Amp.mean(),v_Amp.mean(),r.mean()
	qa,va,ra = q_Amp - qm ,v_Amp - vm   ,r - rm
	total = ((q_Amp*v_Amp)*r).mean()/2.
	term1 = qm*vm*rm/2.
	term2 = rm*((qa*va).mean())/2.
	term3 = qm*((va*ra).mean())/2.
	term4 = vm*((qa*ra).mean())/2.
	term5 = (ra*qa*va).mean()/2.
	tot   = term1 + term2 + term3 + term4 + term5
	"""
	print total
	print '%s = %s + %s + %s + %s + %s' % (tot,term1,term2,term3,term4,term5)
	pl.bar(np.arange(6),[tot,term1,term2,term3,term4,term5],align='center')
	pl.xticks(np.arange(6),['total',r'$\bar{v} \bar{q} \bar{r}$',r'$\bar{r} \overline{v^{\prime} q^{\prime}}$',\
                                                             r'$\bar{q} \overline{v^{\prime} r^{\prime}}$',
                                                             r'$\bar{v} \overline{q^{\prime} r^{\prime}}$',
                                                             r'$\overline{v^{\prime} q^{\prime} r^{\prime}}$'])
	pl.show()
	"""
	return tot,term1,term2,term3,term4,term5

# Initial file for axes, plus other attributes
masktype    = None
years       = range(1980,2005+1,1)
dayinds     = range(60) + range(335,364+1,1)
latx        = 75
levx        = range(9,9+1,1)
wvns        = np.array([1])
Phi_edges   = np.linspace(-np.pi,np.pi,11)
qAmp_edges = np.linspace(0,12.5e-04,11)
vAmp_edges = np.linspace(0,7,11)		# all F
#qAmp_edges  = np.linspace(0,17.5e-04,11)	# F + sigma
#vAmp_edges  = np.linspace(0,9,11)
vqAmp_edges = np.linspace(0,4e-03,11)
cos_edges   = np.linspace(-1,1,11)
dv          = Dataset('../v_decomp/v_decomp_1980.nc','r')
time        = dv.variables['time'][:]
lat         = dv.variables['lat'][:]
lev         = dv.variables['lev'][:]
wvn         = dv.variables['wavenumber'][:]
dv.close()
# ERAInt DataServer
dsERA = reDataServer(Field='T2',LevType='surface_analysis')
# Mass of girdbox on each level
dP_ = []
dP  = np.diff(lev)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665
# ERAInt PDF
q_Amp,q_Phi = np.zeros(0),np.zeros(0)
v_Amp,v_Phi = np.zeros(0),np.zeros(0)
datelist    = []
print 'ERAInt'
for year in years:
        # Data files
        dq = Dataset('../q_decomp/q_decomp_%s.nc' % (year),'r')
        dv = Dataset('../v_decomp/v_decomp_%s.nc' % (year),'r')
        # Extract variables
        qAmp = dq.variables['Amp'][dayinds,levx,latx,wvns].reshape(-1)
        qPhi = dq.variables['Phi'][dayinds,levx,latx,wvns].reshape(-1)
        vAmp = dv.variables['Amp'][dayinds,levx,latx,wvns].reshape(-1)
        vPhi = dv.variables['Phi'][dayinds,levx,latx,wvns].reshape(-1)
	# Add dates to datelist
	datelist = datelist + [dsERA.getDate(24*day)[0:3] for day in dq.variables['time'][dayinds]]
        # Close files
        dq.close()
        dv.close()
        # Append data
        q_Amp = np.append(q_Amp,qAmp,axis=0)
        q_Phi = np.append(q_Phi,qPhi,axis=0)
        v_Amp = np.append(v_Amp,vAmp,axis=0)
        v_Phi = np.append(v_Phi,vPhi,axis=0)
vq_Amp  = v_Amp*q_Amp/2.
vq_dphi = np.cos(v_Phi-q_Phi)
# Mask for high flux events, injections existing etc.
if masktype == 'sigma':
	F     = (q_Amp*v_Amp)*np.cos(q_Phi-v_Phi)/2.
	Fmean = F.mean()
	Fstd  = F.std()
	xs   = np.where(F>Fmean+Fstd)[0]
elif masktype == 'injection':
	G,Q,D = unpick('/qnap/cian/cmip/intrusions/ERAInt_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p')
	E     = np.array([sum([sum(q1) for q1 in q]) for q in Q])
	D     = [D[i] for i in np.where(E>=stats.scoreatpercentile(E,75))[0]]
	D     = [D[i][j][0:3] for i in range(len(D)) for j in range(len(D[i]))]
	xs_   = [datelist.index(injdate) for injdate in D if datelist.count(injdate)>0]
	xs    = []
	for i in range(len(xs_)):
		if xs.count(xs_[i]) == 0:
			xs.append(xs_[i])
	#xs = [i for i in range(len(q_Amp)) if i not in xs]
elif masktype == None:
	xs = range(len(q_Amp))
q_Amp,v_Amp,q_Phi,v_Phi = q_Amp[xs],v_Amp[xs],q_Phi[xs],v_Phi[xs]
vq_Amp,vq_dphi          = vq_Amp[xs],vq_dphi[xs]
print len(datelist),len(xs),1.*len(xs)/len(datelist)
#totE,t1E,t2E,t3E,t4E,t5E = decomp(q_Amp,v_Amp,q_Phi,v_Phi)

# Bivariate PDFs
xeq,yeq,HqERA = bivarPDF( q_Amp,  q_Phi, qAmp_edges,  Phi_edges, norm=False)
xev,yev,HvERA = bivarPDF( v_Amp,  v_Phi, vAmp_edges,  Phi_edges, norm=False)
xe0,ye0,H0ERA = bivarPDF( v_Amp,  q_Amp, vAmp_edges, qAmp_edges, norm=False)
xe1,ye1,H1ERA = bivarPDF( v_Phi,  q_Phi,  Phi_edges,  Phi_edges, norm=False)
xe2,ye2,H2ERA = bivarPDF(vq_Amp,vq_dphi, vqAmp_edges, cos_edges, norm=False)

# Cycle through Models
Models = [g[9:] for g in glob.glob('../rcp85/*')]
HvCMIP,HqCMIP = [],[]
HH0,HH1,HH2   = [],[],[]
#Terms        = []
for Model in Models[:1]:
	print Model
	# CMIP5 DataServer
	dsCMIP = cmipDataServer(Field='va',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='day')
	Dir = '/mnt/climstorage/cian/rcp85/%s/day/plev' % (Model)
	q_Amp,q_Phi = np.zeros(0),np.zeros(0)
	v_Amp,v_Phi = np.zeros(0),np.zeros(0)
	datelist    = []
	for year in years:
		# Data files
		dq = Dataset('%s/q_decomp/q_decomp_%s.nc' % (Dir,year),'r')
		dv = Dataset('%s/v_decomp/v_decomp_%s.nc' % (Dir,year),'r')
		# Extract variables
		qAmp = dq.variables['Amp'][dayinds,levx,latx,wvns].reshape(-1)#*(10000/9.8)
		qPhi = dq.variables['Phi'][dayinds,levx,latx,wvns].reshape(-1)
		vAmp = dv.variables['Amp'][dayinds,levx,latx,wvns].reshape(-1)
		vPhi = dv.variables['Phi'][dayinds,levx,latx,wvns].reshape(-1)
		# Add dates to datelist
		datelist = datelist + [dsCMIP.getDate(day)[0:3] for day in dq.variables['time'][dayinds]]
		# Close files
		dq.close()
		dv.close()
		# Append data
		q_Amp = np.append(q_Amp,qAmp,axis=0)
		q_Phi = np.append(q_Phi,qPhi,axis=0)
		v_Amp = np.append(v_Amp,vAmp,axis=0)
		v_Phi = np.append(v_Phi,vPhi,axis=0)
	vq_Amp  = v_Amp*q_Amp/2.
	vq_dphi = np.cos(v_Phi-q_Phi)
	# Mask
	if masktype == 'sigma':
		F     = (q_Amp*v_Amp)*np.cos(q_Phi-v_Phi)/2.
		Fmean = F.mean()
		Fstd  = F.std()
		xs   = np.where(F>Fmean+Fstd)[0]
	elif masktype == 'injection':
		G,Q,D = unpick('/qnap/cian/cmip/intrusions/%s_intrusions.DJF.2x24hr.9deg.240.24dt.20.15.p' % (Model))
		E     = np.array([sum([sum(q1) for q1 in q]) for q in Q])
		D     = [D[i] for i in np.where(E>=stats.scoreatpercentile(E,75))[0]]
		D     = [D[i][j][0:3] for i in range(len(D)) for j in range(len(D[i]))]
		xs_   = [datelist.index(injdate) for injdate in D if datelist.count(injdate)>0]
		xs    = []
		for i in range(len(xs_)):
		        if xs.count(xs_[i]) == 0:
		                xs.append(xs_[i])
		#xs = [i for i in range(len(q_Amp)) if i not in xs]
	elif masktype == None:
        	xs = range(len(q_Amp))
	q_Amp,v_Amp,q_Phi,v_Phi = q_Amp[xs],v_Amp[xs],q_Phi[xs],v_Phi[xs]
	vq_Amp,vq_dphi          = vq_Amp[xs],vq_dphi[xs]
	print len(datelist),len(xs),1.*len(xs)/len(datelist)
	#totC,t1C,t2C,t3C,t4C,t5C = decomp(q_Amp,v_Amp,q_Phi,v_Phi)
	#Terms.append([totC,t1C,t2C,t3C,t4C,t5C])
	# Bivariate PDFs
	xeq,yeq,Hq = bivarPDF(q_Amp,q_Phi,qAmp_edges,Phi_edges,norm=False)
	xev,yev,Hv = bivarPDF(v_Amp,v_Phi,vAmp_edges,Phi_edges,norm=False)
	xe0,ye0,H0 = bivarPDF(v_Amp,q_Amp,vAmp_edges,qAmp_edges,norm=False)
	xe1,ye1,H1 = bivarPDF(v_Phi,q_Phi, Phi_edges, Phi_edges,norm=False)
	xe2,ye2,H2 = bivarPDF(vq_Amp,vq_dphi, vqAmp_edges, cos_edges, norm=False)
	HvCMIP.append(Hv)
	HqCMIP.append(Hq)
	HH0.append(H0)
	HH1.append(H1)
	HH2.append(H2)

"""
Terms         = np.array(Terms)
Terms_m       = Terms.mean(axis=0)
Terms_s       = Terms.std(axis=0)

totC,t1C,t2C,t3C,t4C,t5C = Terms_m[:]
dtot,dt1,dt2,dt3,dt4,dt5 = totC-totE,t1C-t1E,t2C-t2E,t3C-t3E,t4C-t4E,t5C-t5E
pl.figure(1)
pl.bar(np.arange(6),[dtot,dt1,dt2,dt3,dt4,dt5],align='center')
pl.xticks(np.arange(6),['total',r'$\bar{v} \bar{q} \bar{r}$',r'$\bar{r} \overline{v^{\prime} q^{\prime}}$',\
                                                     r'$\bar{q} \overline{v^{\prime} r^{\prime}}$',
                                                     r'$\bar{v} \overline{q^{\prime} r^{\prime}}$',
                                                     r'$\overline{v^{\prime} q^{\prime} r^{\prime}}$'])
pl.grid()

pl.figure(2)
pl.bar(np.arange(6),[totE,t1E,t2E,t3E,t4E,t5E],align='center')
pl.xticks(np.arange(6),['total',r'$\bar{v} \bar{q} \bar{r}$',r'$\bar{r} \overline{v^{\prime} q^{\prime}}$',\
                                                     r'$\bar{q} \overline{v^{\prime} r^{\prime}}$',
                                                     r'$\bar{v} \overline{q^{\prime} r^{\prime}}$',
                                                     r'$\overline{v^{\prime} q^{\prime} r^{\prime}}$'])
pl.grid()
pl.show()
"""

HvCMIP,HqCMIP  = np.array(HvCMIP),np.array(HqCMIP)
HH0,HH1,HH2    = np.array(HH0),np.array(HH1),np.array(HH2)
dHv,dHq        = HvCMIP-HvERA[np.newaxis,:,:],HqCMIP-HqERA[np.newaxis,:,:]
dH0,dH1,dH2    = HH0 - H0ERA[np.newaxis,:,:],HH1 - H1ERA[np.newaxis,:,:],HH2 - H2ERA[np.newaxis,:,:]
dHvm,dHqm      = dHv.mean(axis=0),dHq.mean(axis=0)
dH0m,dH1m,dH2m = dH0.mean(axis=0),dH1.mean(axis=0),dH2.mean(axis=0)
HvCMIP,HqCMIP  = HvCMIP.mean(axis=0),HqCMIP.mean(axis=0)
HH0,HH1,HH2    = HH0.mean(axis=0),HH1.mean(axis=0),HH2.mean(axis=0)

dHqm,xq,yq   = interp2d(  dHqm,yeq,xeq,3,'linear')
dHvm,xv,yv   = interp2d(  dHvm,yev,xev,3,'linear')
HqCMIP,xq,yq = interp2d(HqCMIP,yeq,xeq,3,'linear')
HvCMIP,xv,yv = interp2d(HvCMIP,yev,xev,3,'linear')
HqERA,xq,yq  = interp2d( HqERA,yeq,xeq,3,'linear')
HvERA,xv,yv  = interp2d( HvERA,yev,xev,3,'linear')
H0ERA,x0,y0  = interp2d( H0ERA,ye0,xe0,3,'linear')
H1ERA,x1,y1  = interp2d( H1ERA,ye1,xe1,3,'linear')
H2ERA,x2,y2  = interp2d( H2ERA,ye2,xe2,3,'linear')
dH0m,x0,y0   = interp2d(  dH0m,ye0,xe0,3,'linear')
dH1m,x1,y1   = interp2d(  dH1m,ye1,xe1,3,'linear')
dH2m,x2,y2   = interp2d(  dH2m,ye2,xe2,3,'linear')
HH0,x0,y0    = interp2d(   HH0,ye0,xe0,3,'linear')
HH1,x1,y1    = interp2d(   HH1,ye1,xe1,3,'linear')
HH2,x2,y2    = interp2d(   HH2,ye2,xe2,3,'linear')

dH2_ = np.array([interp2d(dH2_0,ye2,xe2,3,'linear')[0] for dH2_0 in dH2])
stipx,stipy = stipling(dH2_,x=x2,y=y2,xx=None,yy=None,thresh=0.8)

# Plot
cseq1f  = np.arange(-25,25+5,5)
cseq1l  = np.append(np.arange(-25,0,5),np.arange(5,25+5,5),axis=0)
cseq2   = np.arange(10,100+10,10)

"""
fig,ax1 = pl.subplots(num=1)
ax2     = ax1.twinx()
cf      = ax1.contourf(xq,yq,dHqm,cseq1f,cmap=pl.cm.RdBu_r,extend='both')
cl      = ax2.contour(xv,yv,dHvm,cseq1l,colors='0.65',linewidths=1.2)
cl0      = ax2.contour(xv,yv,dHvm,[0],colors='0.60',linewidths=2)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax1.set_xlabel('Phase [rad]')
ax1.set_ylabel('Amplitude [kg/kg]')
ax2.set_ylabel('Amplitude [m s$^{-1}$]')
ax1.set_title('CMIP5 - ERAInt; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/all.bias.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

fig,ax1 = pl.subplots(num=2)
ax2     = ax1.twinx()
cf      = ax1.contourf(xq,yq,HqCMIP,cseq2,cmap=pl.cm.OrRd,extend='both')
cl      = ax2.contour(xv,yv,HvCMIP,cseq2,colors='0.65',linewidths=1.2)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax1.set_xlabel('Phase [rad]')
ax1.set_ylabel('Amplitude [kg/kg]')
ax2.set_ylabel('Amplitude [m s$^{-1}$]')
ax1.set_title('CMIP5; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/all.full.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

fig,ax1 = pl.subplots(num=3)
ax2     = ax1.twinx()
cf      = ax1.contourf(xq,yq,HqERA,cseq2,cmap=pl.cm.OrRd,extend='both')
cl      = ax2.contour(xv,yv,HvERA,cseq2,colors='0.65',linewidths=1.2)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax1.set_xlabel('Phase [rad]')
ax1.set_ylabel('Amplitude [kg/kg]')
ax2.set_ylabel('Amplitude [m s$^{-1}$]')
ax1.set_title('ERAInt; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/ERAInt.full.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')
"""

fig = pl.figure(4)
ax  = fig.add_subplot(111)
axR = ax.twinx()
axT = ax.twiny()
cf  = ax.contourf(x0,y0,dH0m,cseq1f,cmap=pl.cm.RdBu_r,extend='both')
cl  = ax.contour(x0,y0,dH1m,cseq1l,colors='0.65',linewidths=1.2)
cl0 = ax.contour(x0,y0,dH1m,[0],colors='0.60',linewidths=2)
ax.plot(x0,y0+np.pi/2,'g--',linewidth=1.2,alpha=0.5)
ax.plot(x0,y0-np.pi/2,'g--',linewidth=1.2,alpha=0.5)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.78)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax.set_xlabel('Amplitude [kg/kg]')
ax.set_ylabel('Amplitude [m s$^{-1}$]')
axR.set_ylabel('Phase of v [rad]')
axT.set_xlabel('Phase of q [rad]')
axR.set_yticks(y1)
axR.set_yticks(np.arange(-2.8,2.8+0.8,0.8))
axT.set_xticks(x1)
axT.set_xticks(np.arange(-2.8,2.8+0.8,0.8))
ax.set_ylim(yv[0],yv[-1])
ax.set_xlim(x0[0],x0[-1])
pl.savefig('figs/bivar/all.bias.inv.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

fig = pl.figure(5)
ax  = fig.add_subplot(111)
axR = ax.twinx()
axT = ax.twiny()
cf  = ax.contourf(x0,y0,HH0,11,cmap=pl.cm.OrRd,extend='both')
cl  = ax.contour(x0,y0,HH1,11,colors='0.65',linewidths=1.2)
ax.plot(x0,y0+np.pi/2,'g--',linewidth=1.2,alpha=0.5)
ax.plot(x0,y0-np.pi/2,'g--',linewidth=1.2,alpha=0.5)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.78)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax.set_xlabel('Amplitude [kg/kg]')
ax.set_ylabel('Amplitude [m s$^{-1}$]')
axR.set_ylabel('Phase of v [rad]')
axT.set_xlabel('Phase of q [rad]')
axR.set_yticks(y1)
axR.set_yticks(np.arange(-2.8,2.8+0.8,0.8))
axT.set_xticks(x1)
axT.set_xticks(np.arange(-2.8,2.8+0.8,0.8))
ax.set_ylim(yv[0],yv[-1])
ax.set_xlim(x0[0],x0[-1])
pl.savefig('figs/bivar/all.full.inv.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

fig = pl.figure(6)
ax  = fig.add_subplot(111)
axR = ax.twinx()
axT = ax.twiny()
cf  = ax.contourf(x0,y0,H0ERA,11,cmap=pl.cm.OrRd,extend='both')
cl  = ax.contour(x0,y0,H1ERA,11,colors='0.65',linewidths=1.2)
ax.plot(x0,y0+np.pi/2,'g--',linewidth=1.2,alpha=0.5)
ax.plot(x0,y0-np.pi/2,'g--',linewidth=1.2,alpha=0.5)
pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
fig.subplots_adjust(right=0.78)
cbar_ax = fig.add_axes([0.88, 0.1, 0.03, 0.8])
cbar    = fig.colorbar(cf, cax=cbar_ax)
cbar.set_label('Frequency')
ax.set_xlabel('Amplitude [kg/kg]')
ax.set_ylabel('Amplitude [m s$^{-1}$]')
axR.set_ylabel('Phase of v [rad]')
axT.set_xlabel('Phase of q [rad]')
axR.set_yticks(y1)
axR.set_yticks(np.arange(-2.8,2.8+0.8,0.8))
axT.set_xticks(x1)
axT.set_xticks(np.arange(-2.8,2.8+0.8,0.8))
ax.set_ylim(yv[0],yv[-1])
ax.set_xlim(x0[0],x0[-1])
pl.savefig('figs/bivar/ERAInt.full.inv.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

"""
xx,yy = np.meshgrid(x2,y2)
zz    = xx*yy
zz    = np.append(zz[:,:],zz[:,:][:,::-1],axis=1)
dH2m  = (10000/9.8)*zz*np.append(dH2m[:,:],dH2m[:,:][:,::-1],axis=1)/(2*90.*len(years))
H2ERA = np.append(H2ERA[:,:],H2ERA[:,:][:,::-1],axis=1)/2.
HH2   = np.append(HH2[:,:],HH2[:,:][:,::-1],axis=1)/2.
x2    = range(1,60+1,1)

print (2260e03/1e06)*((10000/9.8)*(zz*H2ERA)).sum()/H2ERA.sum()
print (2260e03/1e06)*((10000/9.8)*(  zz*HH2)).sum()/HH2.sum()

print (2260e03/1e06)*((10000/9.8)*(zz*H2ERA)).sum()/(90.*len(years))
print (2260e03/1e06)*((10000/9.8)*(  zz*HH2)).sum()/(90.*len(years))

pl.figure(7)
print (2260e03/1e06)*dH2m.sum()
cf   = pl.contourf(x2,y2,dH2m,np.arange(-0.035,0.035,0.005),cmap=pl.cm.RdBu_r,extend='both')
cbar = pl.colorbar(cf)
#cl = pl.contour(x2,y2,dH2m,np.append(np.arange(-70,-10+10,10),np.arange(10,70+10,10),axis=0),colors='0.65',linewidths=1.2)
#pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
cl0 = pl.contour(x2,y2,dH2m,[0],colors='0.60',linewidths=2)
#pl.plot(stipx,stipy,'k.',alpha=0.5)
pl.plot([14.5,14.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.plot([46.5,46.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
cbar.set_label('Zonal mean moisture flux [kg s$^{-1}$ m$^{-1}$]')
pl.xlabel(r'$cos(d\phi)$')
pl.xticks([1,14.5,30.5,46.5,60],[-1,0,1,0,-1])
pl.ylabel(r'$\frac{vq}{2}$ [m s$^{-1}$ kg kg$^{-1}$]')
pl.title('CMIP5 - ERAInt; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/all.bias.inv2.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

pl.figure(8)
#cl  = pl.contour(x2,y2,H2ERA,np.arange(15,225+15,15),colors='0.65',linewidths=1.2)
#pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
cf   = pl.contourf(x2,y2,H2ERA,12,cmap=pl.cm.OrRd,extend='max')
cbar = pl.colorbar(cf)
cbar.set_label('Frequency')
pl.plot([14.5,14.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.plot([46.5,46.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.xlabel(r'$cos(d\phi)$')
pl.xticks([1,14.5,30.5,46.5,60],[-1,0,1,0,-1])
pl.ylabel(r'$\frac{vq}{2}$ [m s$^{-1}$ kg kg$^{-1}$]')
pl.title('ERAInt; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/ERAInt.full.inv2.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

pl.figure(9)
#cl  = pl.contour(x2,y2,HH2,np.arange(15,225+15,15),colors='0.65',linewidths=1.2)
#pl.clabel(cl,fmt='%0.1f',colors='k',fontsize=9)
cf   = pl.contourf(x2,y2,HH2,12,cmap=pl.cm.OrRd,extend='max')
cbar = pl.colorbar(cf)
cbar.set_label('Frequency')
pl.plot([14.5,14.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.plot([46.5,46.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.xlabel(r'$cos(d\phi)$')
pl.xticks([1,14.5,30.5,46.5,60],[-1,0,1,0,-1])
pl.ylabel(r'$\frac{vq}{2}$ [m s$^{-1}$ kg kg$^{-1}$]')
pl.title('CMIP5; k = %s; level = %s-%shPa; lat = %sN' % (wvns+1,lev[levx[0]],lev[levx[-1]],lat[latx]))
pl.savefig('figs/bivar/all.full.inv2.%s-%sk.%d-%dhPa.%sN.pdf' % (wvns[0]+1,wvns[-1]+1,lev[levx[0]],lev[levx[-1]],lat[latx]),format='pdf')

pl.figure(10)
pl.contourf(x2,y2,zz,14,cmap=pl.cm.RdBu_r,extend='both')
pl.colorbar()
pl.plot([14.5,14.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.plot([46.5,46.5],[y2[0],y2[-1]],'g--',linewidth=1.2,alpha=0.5)
pl.xlabel(r'$cos(d\phi)$')
pl.xticks([1,14.5,30.5,46.5,60],[-1,0,1,0,-1])
pl.ylabel(r'$\frac{vq}{2}$ [m s$^{-1}$ kg kg$^{-1}$]')
pl.savefig('figs/bivar/zm.pdf',format='pdf')
"""

pl.show()
