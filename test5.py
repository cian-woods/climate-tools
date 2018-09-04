from ReanalysisMonthlyMeanDataServer import *
from LambertProjector import *

import matplotlib.pyplot as pl

from netCDF4 import *

def plotField(x,y,Field,title='',sname='default',cseq=np.arange(0,8+1,1),cmap=pl.cm.RdBu_r,extend='both',clabel='',stipxy=None,Field0=None,cseq0=None,cmap0=None,l0=0,l1=1):

        fig, ax = pl.subplots(1,1)
        cf      = ax.pcolormesh(x-0.5,y,Field,cmap=cmap,vmin=cseq[0],vmax=cseq[-1])
        if Field0 != None:
                if cmap0 != 'red/blue':
                        if list(cseq0).count(0) > 0: cseq0 = np.delete(cseq0,list(cseq0).index(0))
                        cl = ax.contour(x,y,Field0,cseq0,linewidths=1.5,cmap=cmap_trunc(cmap0,l0,l1),alpha=0.85)
                        ax.clabel(cl,fmt='%1.1f',inline=1,colors='k',fontsize=9)
                else:
                        i0 = np.argmin((cseq0 - 0)**2)
                        clb = ax.contour(x,y,Field0,cseq0[:i0],linewidths=1.5,colors='b',alpha=0.85)
                        clr = ax.contour(x,y,Field0,cseq0[i0+1:],linewidths=1.5,colors='r',alpha=0.85)
                        ax.clabel(clr,fmt='%1.1f',inline=1,colors='k',fontsize=9)
                        ax.clabel(clb,fmt='%1.1f',inline=1,colors='k',fontsize=9)
        if stipxy != None: ax.plot(stipxy[0][::6,:],stipxy[1][::6,:],'k.',markersize=5,alpha=0.65)
        cbar = pl.colorbar(cf)
        cbar.set_label(clabel)
	cbar.set_ticks(cseq)
	cbar.set_ticklabels(cseq)
        cseq  = [1000,2000,3000,4000,5000,6500,10000,20000]
        poss  = [[11.25,74.0],[11.25,56.0],[11.25,32.5],[9.9,7.5],[7.9,7.5],[6.1,7.5],[3.95,7.5],[2,7.5]]
        xis   = [75,60,40,20,16,14,12,10]
        rots  = [-12.90,-24.75,-47.50,-82.5,-85,-88,-90,-90]
        for i in range(len(cseq)):
                line_string = r'%s km' % (format(cseq[i],","))
                ks          = 2*np.pi*np.cos(np.pi*y/180.)*6370/cseq[i]
                pp,         = ax.plot(ks,y,'k--',linewidth=0.6,alpha=0.25)
                pos         = poss[i]
                ltex        = pl.text(pos[0], pos[1], line_string, size=9, rotation=rots[i], color = pp.get_color(),\
                                       ha="center", va="center")#,bbox = dict(ec='1',fc='1'))
        ax.set_xlabel('Zonal wavenumber')
        ax.set_ylabel('Latitude')
        pl.title('%s' % (title))  
	ax.set_xticks(np.arange(0,12+1,1))
	ax.set_xticklabels(np.arange(0,12+1,1))
        ax.set_xlim(x[0]-0.5,x[-1]-0.5)
        ax.set_ylim(y[0],y[-1])
        pl.savefig('%s' % (sname),format='pdf')
        pl.close()

YearRange0 = (1981,2016)        # ERAInt
dERA       = Dataset('/mnt/climstorage/cian/%s/%s.mon.mean.new.nc' % ('vq','vq'))
sERA       = dERA.variables['vq'][:,:,:]*2260e03/1e06
months     = np.tile(np.arange(1,12+1,1),sERA.shape[0]/12.)
xs         = np.where((months==1)|(months==2)|(months==12)==True)[0]
sERA       = sERA[xs,:,1:12+1+1]
lats       = dERA.variables['lat'][:]
wvns       = dERA.variables['wavenumber'][1:12+1+1]
sERA       = sERA[2:-1,:,:].reshape((-1,3,len(lats),len(wvns))).mean(axis=1)
years      = range(1980,2016+1,1)
x1,x2      = years.index(YearRange0[0]),years.index(YearRange0[1])
sERAclim   = sERA[x1:x2+1,:,:].mean(axis=0)

Model      = 'BNU-ESM'
dcmip      = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/%s/%s.mon.mean.new.nc' % (Model,'vq','vq'))
scmip      = dcmip.variables['vq'][:,:,:]*2260e03/1e06
months     = np.tile(np.arange(1,12+1,1),scmip.shape[0]/12.)
xs         = np.where((months==1)|(months==2)|(months==12)==True)[0]
scmip      = scmip[xs,:,1:12+1+1]
scmip      = scmip[2:-1,:,:].reshape((-1,3,len(lats),len(wvns))).mean(axis=1)
scmipclim  = scmip[x1:x2+1,:,:].mean(axis=0)

# Total climatology
cseq   = 1*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % ('latent')
plotField(x=wvns,y=lats,Field=sERAclim,title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.clim.pdf' % ('vq'),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
cseq   = 1*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % ('latent')
plotField(x=wvns,y=lats,Field=scmipclim-sERAclim,title='%s - ERAInt' % (Model),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/%s.bias.pdf' % ('vq',Model),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='max',clabel=clabel,stipxy=None)
