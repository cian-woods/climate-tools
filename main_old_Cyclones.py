        #methods  = ['M02','M03','M06','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        #for method in methods:
        #       CS = CycloneServer(method=method)
        #       d,lon,lat,FC,GC,LC = CS.makeClim()

        #CS = CycloneServer(method='M13',resolution=80.)
        #CS.ice()
        #h,lon,lat,FC,GC,LC = CS.makeClim()
        #print FC.shape,GC.shape,LC.shape
        #FC,GC,LC = CS.proj(FC.mean(axis=0),lon[0,:],lat[:,0]),CS.proj(GC.mean(axis=0),lon[0,:],lat[:,0]),CS.proj(LC.mean(axis=0),lon[0,:],lat[:,0])
        #CS.plotCounts(FC,LC,savename='',num=1,cseq=20)

        """
        # Cyclone densities in sector and time range
        case        = 'warm'
        when        = 'all'
        lags,sector = [3,7],[[-180,180],[80,90]]
        datelist    = unpick('dates.%s.50.ERAInt.p' % (case))[:50]
        #methods     = ['M02','M03','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        methods    = ['M13']
        LAT0,LON0   = [],[]
        LAT1,LON1   = [],[]
        LAT2,LON2   = [],[]
        weights     = []
        FC,GC,LC    = [],[],[]
        Months      = [11,12,1,2,3]
        for method in methods:
                CS                             = CycloneServer(method=method)
                h,lon,lat,FCclim,GCclim,LCclim = CS.makeClim()
                xs                             = [i for i in range(len(CS.datelist)) if CS.datelist[i][1] in Months]
                Fclim,Gclim,Lclim              = FCclim[xs,:,:].mean(axis=0),GCclim[xs,:,:].mean(axis=0),LCclim[xs,:,:].mean(axis=0)
                Cyclones0                      = CS.cyclonesInDateListAndSector(CS.Cyclones,datelist,lags=lags,sector=sector,when=when)
                weights.append(1)
                LAT0.append([cyclone.lats[0]  for cyclone in Cyclones0])
                LON0.append([cyclone.lons[0]  for cyclone in Cyclones0])
                LAT1.append([cyclone.lats[-1] for cyclone in Cyclones0])
                LON1.append([cyclone.lons[-1] for cyclone in Cyclones0])
                lat_,lon_ = [],[]
                for cyclone in Cyclones0:
                        for i in range(len(cyclone.lats)):
                                lat_.append(cyclone.lats[i])
                                lon_.append(cyclone.lons[i])
                LAT2.append(lat_)
                LON2.append(lon_)
                FC.append(Fclim)
                GC.append(Gclim)
                LC.append(Lclim)
        FC,GC,LC = np.array(FC).mean(axis=0),np.array(GC).mean(axis=0),np.array(LC).mean(axis=0)
        weights = 1./np.array(weights)
        # set up lat/lon grid
        res           = 3.
        nlon          = 360/res
        nlat          = 180/res
        lon           = np.arange(nlon)*res - 180.
        lat           = np.arange(nlat+1)*res - 90.
        lon,lat       = np.meshgrid(lon,lat)
        C0,C1,C2      = [],[],[]
        for i in range(len(LAT0)):
                Count0 = np.zeros(lon.shape)
                for j in range(len(LAT0[i])):
                        Count0 += np.where(CS.haversine(LON0[i][j],LAT0[i][j],lon,lat) < 564, 1, 0)
                C0.append(Count0)
        for i in range(len(LAT1)):
                Count1 = np.zeros(lon.shape)
                for j in range(len(LAT1[i])):
                        Count1 += np.where(CS.haversine(LON1[i][j],LAT1[i][j],lon,lat) < 564, 1, 0)
                C1.append(Count1)
        for i in range(len(LAT2)):
                Count2 = np.zeros(lon.shape)
                for j in range(len(LAT2[i])):
                        Count2 += np.where(CS.haversine(LON2[i][j],LAT2[i][j],lon,lat) < 564, 1, 0)
                C2.append(Count2)
        Count0        = (np.array(C0)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count1        = (np.array(C1)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count2        = (np.array(C2)*weights[:,np.newaxis,np.newaxis]).sum(axis=0)/weights.sum()
        Count0,Count1 = CS.proj(Count0,lon[0,:],lat[:,0]),CS.proj(Count1,lon[0,:],lat[:,0])
        Count2        = CS.proj(Count2,lon[0,:],lat[:,0])
        FC,GC,LC      = CS.proj(FC,lon[0,:],lat[:,0]),CS.proj(GC,lon[0,:],lat[:,0]),CS.proj(LC,lon[0,:],lat[:,0])
        #cseqf,cseql,cmap = 14,14,pl.cm.OrRd
        cseqf,cseql,cmap = np.arange(0,8+1,1),np.arange(0,0.8+0.1,0.1),pl.cm.OrRd
        #cseqf,cseql,cmap = np.arange(0,0.35+0.025,0.025),np.arange(0,0.03+0.0025,0.0025),pl.cm.OrRd
        # Plot
        Count0,x,y = CS.interp2d(Count0,CS.x,CS.y,6,kind='linear')
        Count1,x,y = CS.interp2d(Count1,CS.x,CS.y,6,kind='linear')
        Count2,x,y = CS.interp2d(Count2,CS.x,CS.y,6,kind='linear')
        #Count0,Count1   = np.ma.masked_where(Count0<1,Count0),np.ma.masked_where(Count1<1,Count1)
        #Count2          = np.ma.masked_where(Count2<1,Count2)
        pl.figure(1)
        cf   = pl.contourf(x,y,Count0/30.,cseql,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count0.sum()))
        pl.savefig('figs/cyclogensis_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.figure(2) 
        cf   = pl.contourf(x,y,Count1/30.,cseql,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count1.sum()))
        pl.savefig('figs/cyclolysis_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.figure(3)
        cf   = pl.contourf(x,y,Count2/30.,cseqf,cmap=cmap,extend='max')
        cbar = pl.colorbar(cf)
        cbar.set_label(r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('%s' % (Count2.sum()))
        pl.savefig('figs/feature_%s-%sN_%s_%s.pdf' % (sector[1][0],sector[1][1],when,case),format='pdf')
        pl.show()
        """

        from Intrusions import *
        YearRange         = (1979,2000)
        years             = range(YearRange[0],YearRange[1]+1,1)
        IN                = Intrusions(OpenTraj=True)
        Nl                = 4*IN.intrusionLaggedDates(lags=[-10,0],case='warm',rankN=50,plot=False,YearRange=YearRange)
        print Nl.shape
        #methods          = ['M02','M03','M08','M09','M10','M12','M13','M14','M15','M16','M18','M20','M21','M22']
        methods           = ['M13']
        datelist1         = unpick('dates.warm.1000.ERAInt.p')[:50]
        datelist2         = unpick('dates.cold.1000.ERAInt.p')[:50]
        datelist1         = [date for date in datelist1 if date[0] in years]
        datelist2         = [date for date in datelist2 if date[0] in years]
        Months            = [11,12,1,2,3]
        FFC1,FFC2         = [],[]
        GGC1,GGC2         = [],[]
        LLC1,LLC2         = [],[]
        FCLIM,GCLIM,LCLIM = [],[],[]
        for method in methods:
                # CycloneServer and climatologies
                CS                             = CycloneServer(method=method)
                h,lon,lat,FCclim,GCclim,LCclim = CS.makeClim()
                xs                             = [i for i in range(len(CS.datelist)) if (CS.datelist[i][1] in Months) and (CS.datelist[i][0] in years) ]
                Fclim,Gclim,Lclim              = FCclim[xs,:,:].mean(axis=0),GCclim[xs,:,:].mean(axis=0),LCclim[xs,:,:].mean(axis=0)
                # Use only date0s falling in cyclone range
                datelist10 = [date for date in datelist1 if (CS.hourlist[0]<=CS.reds.getHours(*date)<=CS.hourlist[-1])]
                datelist20 = [date for date in datelist2 if (CS.hourlist[0]<=CS.reds.getHours(*date)<=CS.hourlist[-1])]
                print len(datelist10),len(datelist20)
                FC1,FC2    = [],[]
                GC1,GC2    = [],[]
                LC1,LC2    = [],[]
                for date0 in datelist10:
                        x     = CS.datelist.index(date0)
                        fc    = FCclim[x-41+2:x+40+1,:,:]       # -9.75 to 10.0 days; len = 80/4 = 20 days
                        gc    = GCclim[x-41+2:x+40+1,:,:]
                        lc    = LCclim[x-41+2:x+40+1,:,:]
                        FC1.append(fc)
                        GC1.append(gc)
                        LC1.append(lc)
                for date0 in datelist20:
                        x     = CS.datelist.index(date0)
                        fc    = FCclim[x-41+2:x+40+1,:,:]
                        gc    = GCclim[x-41+2:x+40+1,:,:]
                        lc    = LCclim[x-41+2:x+40+1,:,:]
                        FC2.append(fc)
                        GC2.append(gc)
                        LC2.append(lc)
                FC1,FC2 = np.array(FC1).mean(axis=0)-Fclim[np.newaxis,:,:],np.array(FC2).mean(axis=0)-Fclim[np.newaxis,:,:]
                GC1,GC2 = np.array(GC1).mean(axis=0)-Gclim[np.newaxis,:,:],np.array(GC2).mean(axis=0)-Gclim[np.newaxis,:,:]
                LC1,LC2 = np.array(LC1).mean(axis=0)-Lclim[np.newaxis,:,:],np.array(LC2).mean(axis=0)-Lclim[np.newaxis,:,:]
                #FC1,FC2 = np.array(FC1).mean(axis=0),np.array(FC2).mean(axis=0)
                #GC1,GC2 = np.array(GC1).mean(axis=0),np.array(GC2).mean(axis=0)
                #LC1,LC2 = np.array(LC1).mean(axis=0),np.array(LC2).mean(axis=0)
                FFC1.append(FC1)
                FFC2.append(FC2)
                GGC1.append(GC1)
                GGC2.append(GC2)
                LLC1.append(LC1)
                LLC2.append(LC2)
                FCLIM.append(Fclim)
                GCLIM.append(Gclim)
                LCLIM.append(Lclim)
        #Nl               = CS.proj(Nl,lonNl[0,:],latNl[:,0])
        FFC1,FFC2         = np.array(FFC1),np.array(FFC2)
        GGC1,GGC2         = np.array(GGC1),np.array(GGC2)
        LLC1,LLC2         = np.array(LLC1),np.array(LLC2)
        FCLIM,GCLIM,LCLIM = np.array(FCLIM).mean(axis=0),np.array(GCLIM).mean(axis=0),np.array(LCLIM).mean(axis=0)
        FFC1,FFC2         = CS.proj(FFC1,lon[0,:],lat[:,0]),CS.proj(FFC2,lon[0,:],lat[:,0])
        GGC1,GGC2         = CS.proj(GGC1,lon[0,:],lat[:,0]),CS.proj(GGC2,lon[0,:],lat[:,0])
        LLC1,LLC2         = CS.proj(LLC1,lon[0,:],lat[:,0]),CS.proj(LLC2,lon[0,:],lat[:,0])
        FFC1m,FFC2m       = 605*FFC1.mean(axis=0),605*FFC2.mean(axis=0)
        GGC1m,GGC2m       = 605*GGC1.mean(axis=0),605*GGC2.mean(axis=0)
        LLC1m,LLC2m       = 605*LLC1.mean(axis=0),605*LLC2.mean(axis=0)
        #cseqf,cseql      = np.arange(-60,60+10,10),np.arange(-9,9+1.5,1.5)
        cseqf,cseql      = np.arange(-90,90+15,15),np.arange(-12,12+1.5,1.5)
        cmap              = pl.cm.RdBu_r

        segs  = [[-60, -40], [-72, -52], [ -51, -30], [-29, -9]]
        sdays = [[ -5,   0], [ -7,  -3], [-2.5, 2.5], [  3,  7]]
        for ii in range(len(segs)):
                # Mean composite between day -5 to 0
                CS.plotCounts1(FFC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),FFC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='figs/lagged/Feature/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseqf,clabel=r'Cyclone freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
                CS.plotCounts1(GGC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),GGC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='figs/lagged/Genesis/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseql,clabel=r'Genesis freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
                CS.plotCounts1(LLC1m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),LLC2m[segs[ii][0]:segs[ii][1],:,:].mean(axis=0),savename='figs/lagged/Lysis/mean_%s-%s.pdf' % (sdays[ii][0],sdays[ii][1]),\
                                num='day %s to %s' % (sdays[ii][0],sdays[ii][1]),stip1=None,stip2=None,cseq=cseql,clabel=r'Lysis freq. {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)

        # Lagged composites between day -10 to 0 over 2-day windows and mean composite between day -5 to 0
        #for i in range(len(FFC1m)-8):
        #        CS.plotCounts2(FFC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='figs/lagged/Feature/%s.pdf' % (i+1),\
        #                        num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseqf,\
        #                       clabel=r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        #CS.plotCounts2(FFC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='figs/lagged/Feature/mean_int.pdf',\
        #                num='day -5 to 0',stip1=None,stip2=None,cseq=cseqf,clabel=r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        #for i in range(len(GGC1m)-8):
        #        CS.plotCounts2(GGC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='figs/lagged/Genesis/%s.pdf' % (i+1),\
        #                        num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseql,
        #                       clabel=r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        #CS.plotCounts2(GGC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='figs/lagged/Genesis/mean_int.pdf',\
        #                num='day -5 to 0',stip1=None,stip2=None,cseq=cseql,clabel=r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        #for i in range(len(LLC1m)-8):
        #        CS.plotCounts2(LLC1m[i:i+8,:,:].mean(axis=0),Nl[i:i+9,:,:].mean(axis=0),savename='figs/lagged/Lysis/%s.pdf' % (i+1),\
        #                        num='day %s to %s' % (0.25*(i-40),0.25*(i-40)+2),stip1=None,stip2=None,cseq=cseql,\
        #                       clabel=r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)
        #CS.plotCounts2(LLC1m[-20:,:,:].mean(axis=0),Nl[-20:,:,:].mean(axis=0),savename='figs/lagged/Lysis/mean_int.pdf',\
        #                num='day -5 to 0',stip1=None,stip2=None,cseq=cseql,clabel=r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$',cmap=cmap)

        # Plot CLIMs
        cseqf,cseql = np.arange(0,210+15,15),np.arange(0,18+1.5,1.5)
        #cseqf,cseql = np.arange(0,0.35+0.025,0.025),np.arange(0,0.03+0.0025,0.0025)
        FCLIM     = CS.proj(FCLIM,lon[0,:],lat[:,0])
        GCLIM     = CS.proj(GCLIM,lon[0,:],lat[:,0])
        LCLIM     = CS.proj(LCLIM,lon[0,:],lat[:,0])
        FCLIM,x,y = CS.interp2d(FCLIM,CS.x,CS.y,6,kind='linear')
        GCLIM,x,y = CS.interp2d(GCLIM,CS.x,CS.y,6,kind='linear')
        LCLIM,x,y = CS.interp2d(LCLIM,CS.x,CS.y,6,kind='linear')
        cf    = pl.contourf(x,y,605*FCLIM,cseqf,cmap=pl.cm.OrRd,extend='max')
        cbar  = pl.colorbar()
        cbar.set_label(r'Cyclone frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80,85],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('figs/IMILAST.fclim.pdf',format='pdf')
        pl.close()
        cf   = pl.contourf(x,y,605*GCLIM,cseql,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar()
        cbar.set_label(r'Genesis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('figs/IMILAST.gclim.pdf',format='pdf')
        pl.close()
        cf   = pl.contourf(x,y,605*LCLIM,cseql,cmap=pl.cm.OrRd,extend='max')
        cbar = pl.colorbar()
        cbar.set_label(r'Lysis frequency {$2 \times 10^{6}$ km$^{2}$ NDJFM}$^{-1}$')
        drawCoastlinesNoRivers(CS.proj.m)
        CS.proj.m.drawparallels([70,80],latmax=90)
        pl.title('Climatology NDJFM 1989-2010')
        pl.savefig('figs/IMILAST.lclim.pdf',format='pdf')
        pl.close()
