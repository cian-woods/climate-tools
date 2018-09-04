from ReanalysisMonthlyMeanDataServer import *
from LambertProjector import *
from scipy import interpolate

import numpy as np
import matplotlib.pyplot as pl

def seasonalVar(years,Season='NDJFM',blat=80.,Source='CAM4xCO2'):
	# Seasonal mean field
	proj = LambertProjector(boundinglat=blat,resolution=80.)
	mmds = DataServer(Field='Ts',Source=Source)
	s    = np.array([mmds.getSeason(year,Season=Season) for year in years])
	s    = proj(s,mmds.lon,mmds.lat).mean(axis=-1).mean(axis=-1)
	return s

def interpolateND(field,xold,xnew,axis):
        # field in (time)x(lev)x(lon)
        f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=False,fill_value='extrapolate')
        field = f(xnew)
        return field

def heatTransport(Source,years,lats):
	# Parameters
	Re    = 6370e03
	dy    = 2*np.pi*Re/360. # m per degree latitude
	# DataServers
	dNLR = DataServer(Field='ttr',Source=Source)
	dNSR = DataServer(Field='tsr',Source=Source)
	# Heat flux from time mean TOA radiation
	nlr   = interpolateND(np.array([dNLR.getSeason(Year=year,Season='Annual') for year in years]).mean(axis=-1),dNLR.lat,lats,axis=1)*2*np.pi*Re*np.cos(lats[np.newaxis,:]*np.pi/180.)
	nsr   = interpolateND(np.array([dNSR.getSeason(Year=year,Season='Annual') for year in years]).mean(axis=-1),dNSR.lat,lats,axis=1)*2*np.pi*Re*np.cos(lats[np.newaxis,:]*np.pi/180.)
	nor   = nsr-nlr
	iy    = np.argmin((lats-0)**2)
	# TOA imbalance for each year in years in W m**-2
	dF    = ((nor[:,0:-1] + nor[:,1:])/2).mean(axis=1)*180*dy/(4*np.pi*(Re**2))
	# TOA imbalance for each year in zonal-mean format i.e. W m**-1
	dFlat = np.array([dF[i]*2*np.pi*Re*np.cos(lats*np.pi/180.) for i in range(len(dF))])
	# Correct annual mean TOA flux
	nor   = nor - dFlat
	E     = np.array([dy*(lats[i]-lats[0])*nor[:,0:i+1].mean(axis=1) for i in range(0,len(lats),1)])
	ES    = np.array([dy*(lats[i]-lats[0])*nor[:,0:i+1].mean(axis=1) for i in range(1,iy+1,1)])
	EN    = np.array([dy*(lats[::-1][i]-lats[::-1][0])*nor[:,::-1][:,0:i+1].mean(axis=1) for i in range(1,iy+1,1)])[::-1,:]
	# Time averages
	dF   = dF.mean()
	E    = E.mean(axis=1)
	ES   = ES.mean(axis=1)
	EN   = EN.mean(axis=1)
	return E,ES,EN,dF

# Attributes
case  = '_high_continent'
if case == '': cstr = 'low'
if case != '': cstr = 'high'
years = range(1900,1919+1,1)
lats  = np.arange(-90,90+1,1)
lats0 = (lats[0:-1]+lats[1:])/2.
i0    = np.argmin((lats-0)**2)
# Energy fluxes
E4 ,  ES4,  EN4,  dF4 = heatTransport( 'CAM4xCO2%s' % (case),years,lats)
E8 ,  ES8,  EN8,  dF8 = heatTransport( 'CAM8xCO2%s' % (case),years,lats)
E16, ES16, EN16, dF16 = heatTransport('CAM16xCO2%s' % (case),years,lats)

print dF4,dF8,dF16

# Plot
pl.figure(1)
pl.plot(lats0[0:90],  ES4/1e15, 'b', linewidth=1., alpha=0.75, label= r'4xCO2%s; $dF=%s$ W m$^{2}$' % (case,round( dF4,3)))
pl.plot( lats0[90:],  EN4/1e15, 'b', linewidth=1., alpha=0.75)
pl.plot(lats0[0:90],  ES8/1e15, 'g', linewidth=1., alpha=0.75, label= r'8xCO2%s; $dF=%s$ W m$^{2}$' % (case,round( dF8,3)))
pl.plot( lats0[90:],  EN8/1e15, 'g', linewidth=1., alpha=0.75)
pl.plot(lats0[0:90], ES16/1e15, 'r', linewidth=1., alpha=0.75, label=r'16xCO2%s; $dF=%s$ W m$^{2}$' % (case,round(dF16,3)))
pl.plot( lats0[90:], EN16/1e15, 'r', linewidth=1., alpha=0.75)
pl.ylabel('Heat transport [PW]')
pl.xlabel('Latitude')
pl.title('Integrating from both poles')
pl.grid()
pl.xlim(lats[0],lats[-1])
pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.savefig('/mnt/climstorage/cian/scripts/figs/TOA/%s.both_poles.full.pdf' % (cstr),format='pdf')

pl.figure(2)
pl.plot(lats0[0:90],  (ES8-ES4)/1e15, 'k', linewidth=1.75, alpha=0.75, label='8xCO2%s - 4xCO2%s' % (case,case))
pl.plot( lats0[90:],  (EN8-EN4)/1e15, 'k', linewidth=1.75, alpha=0.75)
pl.plot(lats0[0:90], (ES16-ES4)/1e15, 'r', linewidth=1.75, alpha=0.75, label='16xCO2%s - 4xCO2%s' % (case,case))
pl.plot( lats0[90:], (EN16-EN4)/1e15, 'r', linewidth=1.75, alpha=0.75)
pl.ylabel('Heat transport [PW]')
pl.xlabel('Latitude')
pl.title('Integrating from both poles')
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.xlim(lats[0],lats[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/TOA/%s.both_poles.change.pdf' % (cstr),format='pdf')

pl.figure(3)
pl.plot(lats,  (E8-E4)/1e15, 'k', linewidth=1.75, alpha=0.75, label='8xCO2%s - 4xCO2%s' % (case,case))
pl.plot(lats, (E16-E4)/1e15, 'r', linewidth=1.75, alpha=0.75, label='16xCO2%s - 4xCO2%s' % (case,case))
pl.ylabel('Heat transport [PW]')
pl.xlabel('Latitude')
pl.title('Integrating from south pole')
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.xlim(lats[0],lats[-1])
pl.savefig('/mnt/climstorage/cian/scripts/figs/TOA/%s.south_pole.change.pdf' % (cstr),format='pdf')

pl.figure(4)
pl.plot(lats, E4/1e15,'b',linewidth=1.,alpha=0.75, label='4xCO2%s; $dF=%s$ W m$^{2}$' % (case,round(  dF4,3)))
pl.plot(lats, E8/1e15,'g',linewidth=1.,alpha=0.75, label='8xCO2%s; $dF=%s$ W m$^{2}$' % (case,round(  dF8,3)))
pl.plot(lats,E16/1e15,'r',linewidth=1.,alpha=0.75,label='16xCO2%s; $dF=%s$ W m$^{2}$' % (case,round( dF16,3)))
pl.ylabel('Heat transport [PW]')
pl.xlabel('Latitude')
pl.grid()
pl.legend(loc=0,frameon=False,prop={'size':12},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
pl.xlim(lats[0],lats[-1])
pl.title('Integrating from south pole')
pl.savefig('/mnt/climstorage/cian/scripts/figs/TOA/%s.south_pole.full.pdf' % (cstr),format='pdf')

pl.show()
