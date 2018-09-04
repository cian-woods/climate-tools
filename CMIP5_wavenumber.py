def scatter_plot(xs,ys,Models,xlabel,ylabel,sname,xlims=None,ylims=None,showERA=False):
	fig,ax = pl.subplots(num=1)
	if len(xs) - len(Models) == 1:
		slope, intercept, r_value, p_value, std_err = stats.linregress(xs[1:],ys[1:])   # Exclude ERAInt in regression
		# Plot CMIP5 models
		for i in range(len(Models)):
		        ax.plot(xs[i+1],ys[i+1],marker=markers[i],markerfacecolor="None",markeredgecolor=scattermap(norm(ixs[i])),markersize=7,mew=2,label=Models[i],linestyle="None")
		# Plot CMIP5 multi model mean and ERAInt
		ax.plot(np.mean(xs[1:]),np.mean(ys[1:]),marker='o',markeredgecolor="None",markerfacecolor='grey',markersize=10,mew=2,label='Multi-model mean',linestyle="None")
		if showERA: ax.plot(xs[0],ys[0],marker='o',markeredgecolor="None",markerfacecolor='k',markersize=7,mew=2,label='ERA-Interim',linestyle="None")
	else:
		slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
		# Plot CMIP5 models     
		for i in range(len(Models)):
		        ax.plot(xs[i],ys[i],marker=markers[i],markerfacecolor="None",markeredgecolor=scattermap(norm(ixs[i])),markersize=7,mew=2,label=Models[i],linestyle="None")
		# Plot CMIP5 multi model mean
		ax.plot(np.mean(xs[:]),np.mean(ys[:]),marker='o',markeredgecolor="None",markerfacecolor='grey',markersize=10,mew=2,label='Multi-model mean',linestyle="None")
	if xlims == None:
		xmin,xmax = ax.get_xlim()
	else:
		xmin,xmax = xlims
	if ylims == None:
		ymin,ymax = ax.get_ylim()
	else:
		ymin,ymax = ylims
	xs         = np.linspace(xmin,xmax,10)
	line       = [slope*x + intercept for x in xs]
	ax.plot(xs,line,'k--',linewidth=1.25,alpha=0.75)
	#ax.plot([0,0],[ymin,ymax],'k--',linewidth=0.75,alpha=0.5)
	#ax.plot([xmin,xmax],[0,0],'k--',linewidth=0.75,alpha=0.5)
	if p_value <= 0.05: color = '0.00'
	if p_value >  0.05: color = '0.75'
	ax.text(5.0,-2.75,'$r^{2}=%.2f$' % (round(r_value**2,2)),fontweight='bold',fontsize=15,color=color)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	# Shrink current axis by 20% and add legend
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	ax.legend(loc='center left',frameon=True,prop={'size':8},ncol=1,columnspacing=1,handletextpad=0.2,title=None,bbox_to_anchor=(1, 0.5),numpoints=1)
	pl.savefig(sname+'pdf',format='pdf')
	pl.close()

def plot_lines(wvns,total,stat,trans,sname):
	pl.figure(1)
	pl.plot(wvns,total,'k.-',linewidth=2.0,alpha=0.80,label='Total')
	pl.plot(wvns, stat,'r.-',linewidth=1.3,alpha=0.65,label='Stationary')
	pl.plot(wvns,trans,'b.-',linewidth=1.3,alpha=0.65,label='Transient')
	pl.ylabel('Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % ('latent'))
	pl.xlabel('Zonal wavenumber')
	pl.xticks(wvns,[int(i) for i in wvns])
	pl.xlim(wvns[0],wvns[-1])
	pl.grid()
	pl.legend(loc=0,frameon=True,prop={'size':11},ncol=1,columnspacing=1,handletextpad=0.2,title=None,bbox_to_anchor=(1, 0.5),numpoints=1)
	pl.savefig(sname,format='pdf')
	pl.close()

def getRegression(x,Y):
	n0,n1,n2    = Y.shape
	r_val,p_val = np.zeros((n1,n2)),np.zeros((n1,n2))
	m_val       = np.zeros((n1,n2))
	for i in range(n1):
	       for j in range(n2):
	               slope, intercept, r_value, p_value, std_err = stats.linregress(x,Y[:,i,j])
	               m_val[i,j] = slope
	               r_val[i,j] = r_value**2
	               p_val[i,j] = p_value/2
	return m_val,r_val,p_val

# Parameters
case = str(sys.argv[1])
wvN0 = int(sys.argv[2])
wvN  = int(sys.argv[3])
Models   = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')]
if case == 'vq': c_sf,sf,orr = 1, 2260e03/1e06, 'latent'
if case == 'vT': c_sf,sf,orr = 8, 1005./1e06  , 'sensible'
if case == 'vv': c_sf,sf,orr = 8, 0.5/1e06    , 'kinetic'
if case == 'vu': c_sf,sf,orr = 8, 0.5/1e06    , 'kinetic'
YearRange0 = (1981,2016)	# ERAInt
YearRange1 = (1981,2016)	# CMIP5 hist
YearRange2 = (2060,2095)	# CMIP5 future
LonRange   = (0,359)
LatRange   = (-15,15)
# Basemap and LSM
bmap = Basemap(projection='merc',llcrnrlat=LatRange[0],urcrnrlat=LatRange[1],llcrnrlon=LonRange[0],urcrnrlon=LonRange[1],lat_ts=20)
lsm  = lsm_interp(np.arange(LatRange[0],LatRange[1]+1,1),np.arange(LonRange[0],LonRange[1]+1,1))

# ERAInt data
dERA      = Dataset('/mnt/climstorage/cian/%s/%s.mon.mean.new.nc' % (case,case))
dERA_stat = Dataset('/mnt/climstorage/cian/%s/%s.stat.1981-2016.DJF.mon.mean.nc' % (case,case))
sERA      = dERA.variables[case][:,:,:]*sf
sERA_stat = dERA_stat.variables[case][:,:,:].squeeze()*sf

months    = np.tile(np.arange(1,12+1,1),sERA.shape[0]/12.)
xs        = np.where((months==1)|(months==2)|(months==12)==True)[0]
sERA      = sERA[xs,:,wvN0:wvN+1]
sERA_stat = sERA_stat[:,wvN0:wvN+1]
lats      = dERA.variables['lat'][:]
wvns      = dERA.variables['wavenumber'][wvN0:wvN+1]
sERA      = sERA[2:-1,:,:].reshape((-1,3,len(lats),len(wvns))).mean(axis=1)
dERA.close()

# Gradient
zonal_anom = False
latband1   = (-15,15)
latband2   =  (70,90)
latband3   = (-15,15)
# Flux
lat1,lat2   = 60,70
wvn1,wvn2   = int(sys.argv[4]),int(sys.argv[5])
latx1,latx2 = np.argmin((lats-lat1)**2),np.argmin((lats-lat2)**2)
wvnx1,wvnx2 = np.argmin((wvns-wvn1)**2),np.argmin((wvns-wvn2)**2)

# Get ERAInt trend
years                      = range(1980,2016+1,1)
x1,x2                      = years.index(YearRange0[0]),years.index(YearRange0[1])
sERAclim_y                 = sERA[x1:x2+1,:,:]
sERAclim                   = sERA[x1:x2+1,:,:].mean(axis=0)
trendERA,PvalERA           = getTrend(years[x1:x2+1],sERA[x1:x2+1,:,:])
sERA_trans                 = sERAclim - sERA_stat
sERA_stat_yearly           = np.array([Dataset('/mnt/climstorage/cian/%s/%s.stat.%s-%s.DJF.mon.mean.nc'\
                             % (case,case,year,year)).variables[case][:,:,wvN0:wvN+1].squeeze()*sf for year in range(YearRange0[0],YearRange0[1]+1,1)])
trendERAstat ,PvalERAstat  = getTrend(years[x1:x2+1],sERA_stat_yearly)
trendERAtrans              = trendERA - trendERAstat
#trendERAtrans,PvalERAtrans = getTrend(years[x1:x2+1],np.array([Dataset('/mnt/climstorage/cian/%s/%s.trans.%s-%s.DJF.mon.mean.nc' % (case,case,year,year)).variables[case][:] for year in range(YearRange0[0],YearRange0[1]+1,1)]))

# Wavelength in km in lat x wavenumber space (not used anymore)
#wavel = 2*6370*np.cos(np.pi*lats/180.)
latsl = np.linspace(lats[0],lats[-1],100)
wvnsl = np.linspace(wvns[0],wvns[-1],100)
wavel = 2*6370*np.cos(np.pi*latsl/180.)
wavel = np.tile(wavel[:,np.newaxis],(1,len(wvnsl)))/wvnsl[np.newaxis,:]
# Plot
# Trend
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendERA,title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.trend.pdf' % (case),\
	  cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Trend (stationary)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendERAstat,title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.stat.trend.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Trend (transient)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendERAtrans,title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.trans.trend.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)#,Field0=trendERAstat,cseq0=cseq,cmap0=pl.cm.RdBu_r,l0=0,l1=1)
# Total climatology
cseq  = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=sERAclim,title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.clim.pdf' % (case),\
	  cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
# Total convergence climatology
cseq   = c_sf*np.arange(-4,4+0.5,0.5)
clabel = 'Zonal-mean meridional %s heat flux convergence \n[W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(sERAclim),1,2),title='ERAInt',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.conv.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Stationary component
cseq  = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=sERA_stat,title='ERAInt: stationary',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.stat.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
# Stationary residual from yearly files
cseq  = c_sf*np.arange(-8,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=sERA_stat - sERA_stat_yearly.mean(axis=0),title='ERAInt: stationary residual',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.stat.clim.residual.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Transient component
cseq  = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=sERA_trans,title='ERAInt: transient',sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.trans.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None,Field0=sERA_stat,cseq0=c_sf*np.arange(1,5+1,1),cmap0=pl.cm.gist_heat_r,l0=0.15,l1=0.75)
# Line plots
plot_lines(wvns,sERAclim[60:70,:].mean(axis=0),sERA_stat[60:70,:].mean(axis=0),sERA_trans[60:70,:].mean(axis=0),\
	   sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/ERAInt.%s-%sN.clim.pdf' % (case,lat1,lat2))
# Cycle through CMIP5 models
trendCMIP      = []
trendCMIPSTAT  = []
trendCMIPTRANS = []
trendCMIP2     = []
fullCMIP       = []
fullCMIP2      = []
fullCMIPSTAT   = []
fullCMIPTRANS  = []
fullCMIPSTAT2  = []
fullCMIPTRANS2 = []
tropoTrend     = []
wvnSection     = []
tropoClim      = []
wvnClim        = []
wvnClimSTAT    = []
wvnClimTRANS   = []
tropTrend      = []
tropClim       = []
arcClim        = []
Corr           = []
Corr2          = []
wvnClim2       = []
wvnClim2STAT   = []
wvnClim2TRANS  = []
tropoTrend2    = []
tropoClim2     = []
tropTrend2     = []
tropClim2      = []
arcClim2       = []
TasClim        = []
TasClim2       = []
TasReg_y       = []
PatCorr        = []

tas_climE_y,tas_climE,lats_tas,lons_tas = getSurfaceClim('ERAInt',YearRange0,zonal_anom=True,LonRange=LonRange,LatRange=LatRange)
m_valE,r_valE,p_valE                    = getRegression(sERA_stat_yearly[:,latx1:latx2+1,wvnx1:wvnx2+1].mean(axis=1).mean(axis=1),tas_climE_y)
corr_series                             = sERAclim_y[:,latx1:latx2+1,wvnx1:wvnx2+1].mean(axis=1).mean(axis=1)
clim0,trend0,troptrend0,trop0,arc0,corr_seriesE,corr0 = getTropoTrend('ERAInt',years[x1:x2+1],latband1,latband2,latband3,LonRange=LonRange,zonal_anom=zonal_anom,corr_series=corr_series)
TasReg_y.append(m_valE)
tropoTrend.append(trend0)
wvnSection.append(trendERA[latx1:latx2+1,wvnx1:wvnx2+1].mean())
tropoClim.append(clim0)
wvnClim.append(sERAclim[latx1:latx2+1,wvnx1:wvnx2+1].mean())
wvnClimSTAT.append(sERA_stat[latx1:latx2+1,wvnx1:wvnx2+1].mean())
wvnClimTRANS.append(sERA_trans[latx1:latx2+1,wvnx1:wvnx2+1].mean())
tropTrend.append(troptrend0)
tropClim.append(trop0)
arcClim.append(arc0)
Corr.append(corr0)
PatCorr.append(1)

models = []
for Model in Models:
	try:
		# CMIP5 data
		dcmip       = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/%s/%s.mon.mean.nc' % (Model,case,case))
		dcmip_stat  = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/%s/%s.stat.1981-2016.DJF.mon.mean.nc' % (Model,case,case))
		dcmip_stat2 = Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/%s/%s.stat.2060-2095.DJF.mon.mean.nc' % (Model,case,case))
		scmip       = dcmip.variables[case][:,:,:]*sf
		scmip_stat  = dcmip_stat.variables[case][:,:,:].squeeze()*sf
                scmip_stat2 = dcmip_stat2.variables[case][:,:,:].squeeze()*sf
		months      = np.tile(np.arange(1,12+1,1),scmip.shape[0]/12.)
		xs          = np.where((months==1)|(months==2)|(months==12)==True)[0]
		scmip       = scmip[xs,:,wvN0:wvN+1]
		scmip_stat  = scmip_stat[:,wvN0:wvN+1]
		scmip_stat2 = scmip_stat2[:,wvN0:wvN+1]
		dcmip.close()
		scmip = scmip[2:-1,:,:].reshape((-1,3,len(lats),len(wvns))).mean(axis=1)
		# Surface temperature climatology
		tas_clim_y ,tas_clim,lats_tas,lons_tas  = getSurfaceClim(Model,YearRange1,zonal_anom=True,LonRange=LonRange,LatRange=LatRange)
		tas_clim2_y,tas_clim2,lats_tas,lons_tas = getSurfaceClim(Model,YearRange2,zonal_anom=True,LonRange=LonRange,LatRange=LatRange)
		TasClim.append(tas_clim)
		TasClim2.append(tas_clim2)

		# Get CMIP trend and tropo trend (hist)
		years              = range(1981,2100+1,1)
		x1h,x2h            = years.index(YearRange1[0]),years.index(YearRange1[1])
		scmipclim_y        = scmip[x1h:x2h+1,:,:]
		scmipclim          = scmip[x1h:x2h+1,:,:].mean(axis=0)
		scmip_trans        = scmipclim - scmip_stat
		trendcmip,Pvalcmip = getTrend(years[x1h:x2h+1],scmip[x1h:x2h+1,:,:])
                # Stationary and transient trends
                scmip_stat_yearly           = np.array([Dataset('/mnt/climstorage/cian/rcp85/%s/mon/surface/%s/%s.stat.%s-%s.DJF.mon.mean.nc'\
                                              % (Model,case,case,year,year)).variables[case][:,:,wvN0:wvN+1].squeeze()*sf for year in range(YearRange0[0],YearRange0[1]+1,1)])
                trendcmipstat ,Pvalcmipstat = getTrend(years[x1:x2+1],scmip_stat_yearly)
                trendcmiptrans              = trendcmip - trendcmipstat
		# For YearRange2 (future)
                x1f,x2f      = years.index(YearRange2[0]),years.index(YearRange2[1])
		scmipclim2_y = scmip[x1f:x2f+1,:,:]
                scmipclim2   = scmip[x1f:x2f+1,:,:].mean(axis=0)
		scmip_trans2 = scmipclim2 - scmip_stat2

		m_val,r_val,p_val = getRegression(scmip_stat_yearly[:,latx1:latx2+1,wvnx1:wvnx2+1].mean(axis=1).mean(axis=1),tas_clim_y)
		TasReg_y.append(m_val)
		PatCorr.append(np.corrcoef(np.ma.masked_where(lsm==1,m_valE).reshape(-1),np.ma.masked_where(lsm==1,m_val).reshape(-1))[0][1])

		# Append all lists at end 
                fullCMIP.append(scmipclim)
                fullCMIPSTAT.append(scmip_stat)
                fullCMIPTRANS.append(scmip_trans)
                trendCMIP.append(trendcmip)
                trendCMIPSTAT.append(trendcmipstat)
                trendCMIPTRANS.append(trendcmiptrans)
                fullCMIP2.append(scmipclim2)
                fullCMIPSTAT2.append(scmip_stat2)
                fullCMIPTRANS2.append(scmip_trans2)
                trendcmip2,Pvalcmip2 = getTrend(years[x1f:x2f+1],scmip[x1f:x2f+1,:,:])
                trendCMIP2.append(trendcmip2)

		# Hist
		corr_series = scmipclim_y[:,latx1:latx2+1,wvnx1:wvnx2+1].mean(axis=1).mean(axis=1)
		clim0,trend0,troptrend0,trop0,arc0,corr_series0,corr0 = getTropoTrend(Model,years[x1h:x2h+1],latband1,latband2,latband3,LonRange=LonRange,zonal_anom=zonal_anom,corr_series=corr_series)
		tropoTrend.append(trend0)
		wvnSection.append(trendcmip[latx1:latx2+1,wvnx1:wvnx2+1].mean())
		tropoClim.append(clim0)
		wvnClim.append(scmipclim[latx1:latx2+1,wvnx1:wvnx2+1].mean())
		wvnClimSTAT.append(scmip_stat[latx1:latx2+1,wvnx1:wvnx2+1].mean())
		wvnClimTRANS.append(scmip_trans[latx1:latx2+1,wvnx1:wvnx2+1].mean())
		tropTrend.append(troptrend0)
		tropClim.append(trop0)
		arcClim.append(arc0)
		Corr.append(corr0)

		# Future
		corr_series = scmipclim2_y[:,latx1:latx2+1,wvnx1:wvnx2+1].mean(axis=1).mean(axis=1)
                clim0,trend0,troptrend0,trop0,arc0,corr_series0,corr0 = getTropoTrend(Model,years[x1f:x2f+1],latband1,latband2,latband3,LonRange=LonRange,zonal_anom=zonal_anom,corr_series=corr_series)
                tropoTrend2.append(trend0)
                tropoClim2.append(clim0)
                wvnClim2.append(scmipclim2[latx1:latx2+1,wvnx1:wvnx2+1].mean())
                wvnClim2STAT.append(scmip_stat2[latx1:latx2+1,wvnx1:wvnx2+1].mean())
                wvnClim2TRANS.append(scmip_trans2[latx1:latx2+1,wvnx1:wvnx2+1].mean())
                tropTrend2.append(troptrend0)
                tropClim2.append(trop0)
                arcClim2.append(arc0)

		# Multi-model mean regression of interannual wvk flux and surface tropical temperature (model full)	
		m_val = np.ma.masked_where(lsm==1,m_val)
		j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
		cseq  = np.arange(-0.2,0.2+0.025,0.025)
		cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
		bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
		bmap.fillcontinents(color='0.85',zorder=10)
		j0,k0 = bmap(  0,0)
		j1,k1 = bmap(359,0)
		pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
		cbar = pl.colorbar(cf,orientation='horizontal')
		cbar.set_label('Reg. coeff. [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
		pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/models/reg_y.%s.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (Model,case,YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2))
		pl.close()

		print Model
		models.append(Model)
	except:
		print 'Model %s failed to complete ...' % (Model)


Models     = models
scattermap = pl.get_cmap('nipy_spectral')
norm       = colors.Normalize(vmin=0, vmax=1)
markers    = select_markers(len(Models))
ixs        = np.linspace(0.03,0.97,len(Models))

# Express historical climatologies as biases
corrBias     =         np.array(Corr)# -         Corr[0]
tropBias     =     np.array(tropClim) -     tropClim[0]
tropoBias    =    np.array(tropoClim) -    tropoClim[0]
wvnBias      =      np.array(wvnClim) -      wvnClim[0]
wvnBiasSTAT  =  np.array(wvnClimSTAT) -  wvnClimSTAT[0]
wvnBiasTRANS = np.array(wvnClimTRANS) - wvnClimTRANS[0]

xlabel = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.total.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropoBias,wvnBias,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.stat.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropoBias,wvnBiasSTAT,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.trans.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropoBias,wvnBiasTRANS,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel = 'Corr. coeff.'
ylabel = 'Upper tropospheric DJF temperature difference [K]'
sname  = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/corr.%s.total.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
scatter_plot(corrBias,tropBias,Models,xlabel,ylabel,sname,xlims=None,ylims=None,showERA=True)

xlabel = 'Pattern corr. coeff.'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname  = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/pattern_corr.tas.%s-%s.DJF.flux_%s-%sN_%s-%sk.' % (YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2)
scatter_plot(PatCorr,wvnBias,Models,xlabel,ylabel,sname,xlims=None,ylims=None,showERA=True)

xlabel = 'Upper tropospheric DJF temperature difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname  = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.total.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropBias,wvnBias,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel = 'Upper tropospheric DJF temperature difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname  = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.stat.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropBias,wvnBiasSTAT,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel = 'Upper tropospheric DJF temperature difference [K]'
ylabel = '%s energy flux DJF bias [10$^{6}$ W m$^{-1}$]' % (orr)
sname  = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/clim.%s.trans.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange1[0],YearRange1[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
scatter_plot(tropBias,wvnBiasTRANS,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.total.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropClim2) - np.array(tropClim[1:])
dwvnClim  = np.array(wvnClim2)  - np.array(wvnClim[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.stat.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropClim2) - np.array(tropClim[1:])
dwvnClim  = np.array(wvnClim2STAT)  - np.array(wvnClimSTAT[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.trans.%s-%s.DJF.full_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband3[0],latband3[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropClim2) - np.array(tropClim[1:])
dwvnClim  = np.array(wvnClim2TRANS)  - np.array(wvnClimTRANS[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.total.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropoClim2) - np.array(tropoClim[1:])
dwvnClim  = np.array(wvnClim2)  - np.array(wvnClim[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.stat.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropoClim2) - np.array(tropoClim[1:])
dwvnClim  = np.array(wvnClim2STAT)  - np.array(wvnClimSTAT[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

xlabel    = 'Upper tropospheric DJF temperature gradient difference [K]'
ylabel    = '%s energy flux DJF change [10$^{6}$ W m$^{-1}$]' % (orr)
sname     = '/mnt/climstorage/cian/scripts/figs/regress/wavenumber/change.%s.trans.%s-%s.DJF.gradient_%s-%sN_%s-%sN.flux_%s-%sN_%s-%sk.' % (case,YearRange2[0],YearRange2[1],latband1[0],latband1[1],latband2[0],latband2[1],lat1,lat2,wvn1,wvn2)
dtropClim = np.array(tropoClim2) - np.array(tropoClim[1:])
dwvnClim  = np.array(wvnClim2TRANS)  - np.array(wvnClimTRANS[1:])
scatter_plot(dtropClim,dwvnClim,Models,xlabel,ylabel,sname,xlims=(-5.0,7.5),ylims=(-3.0,6.0))

# Trend, multi-model mean and differences with ERAInt
nmodels         = len(trendCMIP)
TasClim         = np.array(TasClim)
TasClim2        = np.array(TasClim2)
TasReg_y        = np.array(TasReg_y)
dTasClim        = TasClim2 - TasClim
trendCMIP       = np.array(trendCMIP)
trendCMIPSTAT   = np.array(trendCMIPSTAT)
trendCMIPTRANS  = np.array(trendCMIPTRANS)
trendCMIP2      = np.array(trendCMIP2)
fullCMIP        = np.array(fullCMIP)
fullCMIP2       = np.array(fullCMIP2)
fullCMIPSTAT    = np.array(fullCMIPSTAT)
fullCMIPTRANS   = np.array(fullCMIPTRANS)
fullCMIPSTAT2   = np.array(fullCMIPSTAT2)
fullCMIPTRANS2  = np.array(fullCMIPTRANS2)
biasCMIP        = fullCMIP - sERAclim[np.newaxis,:,:]
biasCMIPSTAT    = fullCMIPSTAT - sERA_stat[np.newaxis,:,:]
biasCMIPTRANS   = fullCMIPTRANS - sERA_trans[np.newaxis,:,:]
biasCMIPTOTAL   = biasCMIPTRANS + biasCMIPSTAT
changeCMIP      = fullCMIP2 - fullCMIP
changeCMIPSTAT  = fullCMIPSTAT2 - fullCMIPSTAT
changeCMIPTRANS = fullCMIPTRANS2 - fullCMIPTRANS
dtrend          = trendCMIP[:,:,:] - trendERA[np.newaxis,:,:]
stipxdiff,stipydiff               = stipling(biasCMIP       ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxdiffstat,stipydiffstat       = stipling(biasCMIPSTAT   ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxdifftrans,stipydifftrans     = stipling(biasCMIPTRANS  ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxdifftotal,stipydifftotal     = stipling(biasCMIPTOTAL  ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxchange,stipychange           = stipling(changeCMIP     ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxchangestat,stipychangestat   = stipling(changeCMIPSTAT ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxchangetrans,stipychangetrans = stipling(changeCMIPTRANS,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxfull,stipyfull               = stipling(trendCMIP      ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxfull2,stipyfull2             = stipling(trendCMIP2     ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)
stipxdiff,stipydiff               = stipling(dtrend         ,xx=None,yy=None,x=wvns,y=lats,thresh=0.8)

# Multi-model regression with tropical temperature
m_val,r_val,p_val = getRegression(tropClim[1:],biasCMIP)
cseq   = np.arange(0,0.4+0.05,0.05)
clabel = 'r$^{2}$'
plotField_mesh(x=wvns,y=lats,Field=r_val,title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.clim.rsq.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
cseq   = np.arange(-1,1+0.1,0.1)
clabel = 'Reg. coeff. [W m$^{-1}$ K$^{-1}$]'
regx,regy         = np.meshgrid(wvns,lats)
stipregx,stipregy = np.ma.masked_where(p_val>0.05,regx),np.ma.masked_where(p_val>0.05,regy)
plotField_mesh(x=wvns,y=lats,Field=m_val,title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.clim.reg.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipregx,stipregy))

# Multi-model regression of wvk flux and surface tropical temperature biases
m_val,r_val,p_val = getRegression(wvnBiasSTAT[1:],TasClim)
m_val,r_val,p_val = np.ma.masked_where(lsm==1,m_val),np.ma.masked_where(lsm==1,r_val),np.ma.masked_where(lsm==1,p_val)
j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
j0,k0 = bmap(  0,0)
j1,k1 = bmap(359,0)
sx,sy = np.ma.masked_where(p_val>0.01,j),np.ma.masked_where(p_val>0.01,k)
cseq  = np.arange(-0.7,0.7+0.1,0.1)
cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
pl.plot(sx[::6,::6],sy[::6,::6],'k.',alpha=0.5)
bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
bmap.fillcontinents(color='0.85',zorder=10)
pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
cbar = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Reg. coeff. [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/reg.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (case,YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2))
pl.close()
# Multi-model regression of wvk flux and surface tropical temperature changes
dwvnClim          = np.array(wvnClim2STAT)  - np.array(wvnClimSTAT[1:])
m_val,r_val,p_val = getRegression(dwvnClim,dTasClim)
m_val,r_val,p_val = np.ma.masked_where(lsm==1,m_val),np.ma.masked_where(lsm==1,r_val),np.ma.masked_where(lsm==1,p_val)
j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
sx,sy = np.ma.masked_where(p_val>0.01,j),np.ma.masked_where(p_val>0.01,k)
cseq  = np.arange(-0.7,0.7+0.1,0.1)
cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
pl.plot(sx[::6,::6],sy[::6,::6],'k.',alpha=0.5)
bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
bmap.fillcontinents(color='0.85',zorder=10)
pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
cbar = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Reg. coeff. [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/reg.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (case,YearRange2[0],YearRange2[1],lat1,lat2,wvn1,wvn2))
pl.close()
# Multi-model mean regression of interannual wvk flux and surface tropical temperature (ERAInt full)
m_val = np.ma.masked_where(lsm==1,TasReg_y[0])
j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
cseq  = np.arange(-0.2,0.2+0.025,0.025)
cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
bmap.fillcontinents(color='0.85',zorder=10)
pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
cbar = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Reg. coeff. [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/reg_y.ERAInt.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (case,YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2))
pl.close()
# Multi-model mean regression of interannual wvk flux and surface tropical temperature (Multi-model mean full)
m_val = np.ma.masked_where(lsm==1,TasReg_y[1:,:,:].mean(axis=0))
j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
sx,sy = stipling(TasReg_y[1:,:,:],xx=j,yy=k,x=None,y=None,thresh=0.8)
sx,sy = np.ma.masked_where(lsm==1,sx),np.ma.masked_where(lsm==1,sy)
cseq  = np.arange(-0.2,0.2+0.025,0.025)
cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
pl.plot(sx[::6,::6],sy[::6,::6],'k.',alpha=0.5)
bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
bmap.fillcontinents(color='0.85',zorder=10)
pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
cbar = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Reg. coeff. bias [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/reg_y.full.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (case,YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2))
pl.close()
# Multi-model mean regression of interannual wvk flux and surface tropical temperature (Multi-model mean bias)
m_val = TasReg_y[1:,:,:] - TasReg_y[0,:,:]
j,k   = bmap(*np.meshgrid(lons_tas,lats_tas))
sx,sy = stipling(m_val,xx=j,yy=k,x=None,y=None,thresh=0.8)
sx,sy = np.ma.masked_where(lsm==1,sx),np.ma.masked_where(lsm==1,sy)
m_val = np.ma.masked_where(lsm==1,m_val.mean(axis=0))
cseq  = np.arange(-0.2,0.2+0.025,0.025)
cf    = bmap.contourf(j,k,m_val,cseq,cmap=pl.cm.RdBu_r,extend='both')
pl.plot(sx[::6,::6],sy[::6,::6],'k.',alpha=0.5)
bmap.drawcoastlines(color='0.85',linewidth=1,zorder=11)
bmap.fillcontinents(color='0.85',zorder=10)
pl.plot([j0,j1],[k0,k1],'k--',linewidth=1.25,alpha=0.75,zorder=12)
cbar = pl.colorbar(cf,orientation='horizontal')
cbar.set_label('Reg. coeff. bias [K {10$^{6}$ W m$^{-1}$}$^{-1}$]')
pl.savefig('/mnt/climstorage/cian/scripts/figs/regress/wavenumber/reg_y.bias.%s.stat.%s-%s.DJF.full_tas.flux_%s-%sN_%s-%sk.pdf' % (case,YearRange1[0],YearRange1[1],lat1,lat2,wvn1,wvn2))
pl.close()

# Plot model mean full trend and differences with ERAInt
# Multi-model mean
cseq   = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=fullCMIP.mean(axis=0) ,title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.clim.pdf' % (case),\
	  cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
# Multi-model mean convergence
cseq   = c_sf*np.arange(-4,4+0.5,0.5)
clabel = 'Zonal-mean meridional %s heat flux convergence [W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(fullCMIP.mean(axis=0)),1,2),title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.conv.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Multi-model mean stationary
cseq  = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=fullCMIPSTAT.mean(axis=0),title='CMIP5 stationary: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.stat.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None)
# Multi-model mean transient
cseq  = c_sf*np.arange(0,8+1,1)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=fullCMIPTRANS.mean(axis=0),title='CMIP5 transient: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.trans.clim.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.OrRd,extend='max',clabel=clabel,stipxy=None,Field0=fullCMIPSTAT.mean(axis=0),cseq0=c_sf*np.arange(1,5+1,1),cmap0=pl.cm.gist_heat_r,l0=0.15,l1=0.75)
# Multi-model mean bias
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=biasCMIP.mean(axis=0) ,title='CMIP5 - ERAInt: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.bias.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxdiff,stipydiff))
# Multi-model mean convergence bias
cseq   = c_sf*np.arange(-1,1+0.2,0.2)
clabel = 'Zonal-mean meridional %s heat flux convergence [W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(fullCMIP.mean(axis=0)),1,2) - n_point_smooth(conv(sERAclim),1,2),title='CMIP5 - ERAInt: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.conv.bias.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Multi-model mean stationary bias
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=biasCMIPSTAT.mean(axis=0),title='CMIP5 stationary: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.stat.bias.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxdiffstat,stipydiffstat))
# Multi-model mean transient bias
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=biasCMIPTRANS.mean(axis=0),title='CMIP5 transient: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.trans.bias.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxdifftrans,stipydifftrans))#,Field0=biasCMIPSTAT.mean(axis=0),cseq0=c_sf*np.arange(-2,2+0.5,0.5),cmap0=pl.cm.RdBu_r,l0=0,l1=1)
# Multi-model mean transient + stationary bias (should be same as total bias in all.bias.pdf)
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux\n[10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=biasCMIPTOTAL.mean(axis=0),title='CMIP5 transient + stationary: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.total.bias.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxdifftotal,stipydifftotal))
# Multi-model mean change (total)
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=changeCMIP.mean(axis=0) ,title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxchange,stipychange))
# Multi-model mean convergence change (total)
cseq   = c_sf*np.arange(-1,1+0.2,0.2)
clabel = 'Zonal-mean meridional %s heat flux convergence [W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(fullCMIP2.mean(axis=0)),1,2) - n_point_smooth(conv(fullCMIP.mean(axis=0)),1,2),title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.conv.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel)
# Multi-model mean convergence change (stationary)
cseq   = c_sf*np.arange(-1,1+0.2,0.2)
clabel = 'Zonal-mean meridional %s heat flux convergence [W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(fullCMIPSTAT2.mean(axis=0)),1,2) - n_point_smooth(conv(fullCMIPSTAT.mean(axis=0)),1,2),title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.conv.stat.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel)
# Multi-model mean convergence change (transient)
cseq   = c_sf*np.arange(-1,1+0.2,0.2)
clabel = 'Zonal-mean meridional %s heat flux convergence [W m$^{-2}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=n_point_smooth(conv(fullCMIPTRANS2.mean(axis=0)),1,2) - n_point_smooth(conv(fullCMIPTRANS.mean(axis=0)),1,2),title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.conv.trans.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel)
# Multi-model mean change (stationary)
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=changeCMIPSTAT.mean(axis=0) ,title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.stat.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxchangestat,stipychangestat))
# Multi-model mean change (transient)
cseq   = c_sf*np.arange(-2,2+0.25,0.25)
clabel = 'Zonal-mean meridional %s heat flux [10$^{6}$ W m$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=changeCMIPTRANS.mean(axis=0) ,title='CMIP5 - CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.trans.change.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxchangetrans,stipychangetrans),Field0=changeCMIPSTAT.mean(axis=0),cseq0=c_sf*np.arange(-2,2+0.5,0.5),cmap0=pl.cm.RdBu_r,l0=0.65,l1=1)
# Multi-model mean trend (total, historical)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendCMIP.mean(axis=0),title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.trend.hist.pdf' % (case),\
	  cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxfull,stipyfull))
# Multi-model mean trend (stationary, historical)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendCMIPSTAT.mean(axis=0),title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.trend.stat.hist.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Multi-model mean trend (transient, historical)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendCMIPTRANS.mean(axis=0),title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.trend.trans.hist.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=None)
# Multi-model mean trend (future)
cseq   = c_sf*np.arange(-0.35,0.35+0.05,0.05)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=trendCMIP2.mean(axis=0),title='CMIP5: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.trend.future.pdf' % (case),\
          cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxfull2,stipyfull2))
# Multi-model mean trend bias
cseq   = c_sf*np.arange(-1,1+0.2,0.2)
clabel = 'Zonal-mean meridional %s heat flux trend\n[10$^{6}$ W m$^{-1}$ decade$^{-1}$]' % (orr)
plotField_mesh(x=wvns,y=lats,Field=dtrend.mean(axis=0)   ,title='CMIP5 - ERAInt: %s models' % (nmodels),sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.dtrend.pdf' % (case),\
	  cseq=cseq,cmap=pl.cm.RdBu_r,extend='both',clabel=clabel,stipxy=(stipxdiff,stipydiff))
# Line plots
plot_lines(wvns,fullCMIP[:,60:70,:].mean(axis=0).mean(axis=0),fullCMIPSTAT[:,60:70,:].mean(axis=0).mean(axis=0),fullCMIPTRANS[:,60:70,:].mean(axis=0).mean(axis=0),\
           sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/full/all.%s-%sN.clim.pdf' % (case,lat1,lat2))
plot_lines(wvns,biasCMIP[:,60:70,:].mean(axis=0).mean(axis=0),biasCMIPSTAT[:,60:70,:].mean(axis=0).mean(axis=0),biasCMIPTRANS[:,60:70,:].mean(axis=0).mean(axis=0),\
           sname='/mnt/climstorage/cian/scripts/figs/wavenumber/%s/diff/all.%s-%sN.bias.pdf' % (case,lat1,lat2))
