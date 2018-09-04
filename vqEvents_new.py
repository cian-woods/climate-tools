from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from TrackDataServer import DataServer as TRDataServer
from UnPickle import *
from scipy import interpolate,stats
from UnPickle import *
from toPick import *
from netCDF4 import Dataset
from stipling import stipling

import numpy as np
import matplotlib.pyplot as pl
import sys,os

def myround(x, base=5):
	# Rounding function
        return int(base * round(float(x)/base))

def filterInjections(fname,blat,Source,blat0):
	# Filtered filename
	fname_filt = '%s.filtered.%sN.p' % (fname[:-2],blat)
	print fname_filt
	if not os.path.isfile(fname_filt):
		# Open unfiltered injections
		G,Q,D = unpick(fname)
		# TrackDataServer and intrusion trajectories
		td        = TRDataServer(Source=Source,Type='fwrd',blat=blat0,steps=21)
		LON,LAT,P = td.getIntrusionTrajectories(G,D)
        	g,q,d     = [],[],[]
        	for i in range(len(LAT)):
        	        n0,n1 = 0,0
        	        for t in range(len(LAT[i])):
        	                for j in range(len(LAT[i][t])):
        	                        n0 = n0 + 1
        	                        if (LAT[i][t][j] >= blat).any():
        	                                n1 = n1 + 1
        	        if 1.*n1/n0 >= 0.4:
        	                g.append(G[i])
        	                q.append(Q[i])
        	                d.append(D[i])
		toPick([g,q,d],fname_filt)
		print 'Filtered injections from %s to %s' % (len(G),len(g))
	else:
		print 'File %s already exists. Did nothing.' % (fname_filt)
		g,q,d = unpick(fname_filt)
	return g,q,d

def testInj(G,D):
	# Test the injection data for any duplicate points in space-time
	gs = {}
	for i in range(len(D)):
	        for j in range(len(D[i])):
	                try:
	                        gs[D[i][j]].append(G[i][j])
	                except:
	                        gs[D[i][j]] = []
	                        gs[D[i][j]].append(G[i][j])
	Gs = {}
	for i in gs:
	        if len(gs[i]) > 1:
	                inter = list(set.intersection(*map(set,gs[i])))
	                if inter != []:
	                        Gs[i] = [[inte for inte in inter] for kk in range(len(gs[i])-1)]
	if len(Gs) == 0:
		print 'There were no duplicate points: %s injection events' % (len(G))
	else:
		print '%s duplicate points detected' % (len(gs))

def interpolateFluxFile(vq):
	# Interpolates file of shape (ntime,nlon) to (ntime,360)
        n0,n1 = vq.shape
        print 'Interpolating data to 1 degree resolution...'
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
        return vq,np.arange(360)

def getEventIndices(x,flux,min_width,dcrit):
	i0,i1 = False,False
	if x.max() < flux:   return []
	if (x > flux).all(): return [[0,len(x)-1]]
	# find indexes were x crosses flux threshold upward
	up = (np.nonzero( (x[:-1]<flux) & (x[1:]>flux) )[0] + 1).tolist()
	if x[0] > flux: up,i0 = [0] + up, True
	# find indexes were x crosses flux threshold downward
	dn = np.nonzero( (x[:-1]>flux) & (x[1:]<flux) )[0].tolist()
	if x[-1] > flux: dn,i1 = dn + [len(x)-1], True
	# if first crossing is downward, drop it
	if dn[0] < up[0]: dn.pop(0)
	if len(up) > len(dn): up.pop(-1)
	events = zip(up,dn)
	events = [range(event[0],event[1]+1,1) for event in events]
        # If flux crosses 0E then stitch sections together
        if i0 and i1:
                events.append(events[-1]+events[0])
                del events[-2]
                del events[0]
	# Combine sections < dcrit degrees longitude apart
	events = [e for e in events if len(e) >= min_width]
	if len(events) >= 2: events = linkSnap(events,dcrit)
	return events

def linkSections(section0,section1,dcrit):
	sections = [section0,section1]
	ends     = [section0[0],section1[0]]
	x        = np.argmin(ends)
	y        = list(set(range(2))^set([x]))[0]
	if sections[y][0] - sections[x][-1] <= dcrit:
		section  = sections[x] + sections[y]
	else:
		section  = sections[y] + sections[x]
	return section

def indexMatrix(N):
	inds  = [sorted([i,j]) for i in range(N) for j in range(N) if i!=j]
	inds_ = []
	for i in inds:
	        if inds_.count(i) == 0:
	                inds_.append(i)
	return inds_

def linkSnap(events,dcrit):
	# Distance matrix
	inds = indexMatrix(len(events))
	d    = np.array([dist(events[i],events[j]) for i,j in inds])
        if not (d <= dcrit).any():
                return events
	else:
		while d.min() <= dcrit:
			i0,j0  = inds[np.argmin(d)]
			link   = linkSections(events[i0],events[j0],dcrit)
			events = [link] + [events[k] for k in range(len(events)) if k not in [i0,j0]]
			inds   = indexMatrix(len(events))
			d      = np.array([dist(events[i],events[j]) for i,j in inds])
			if len(d) == 0:
				return events
	return events

def linkStep(events0,events1,prop):
	# Pairs the closest sections in each snapshot (unique pairings)
	# Can only pair N = min(len(events0),len(events1)) sections
	ind,d = [],[]
	for i in range(len(events0)):
		for j in range(len(events1)):
			ind.append([i,j])
			d.append(dist(events0[i],events1[j]))
	N       = np.min([len(events0),len(events1)])
	c       = 0
	dE,indE = [],[]
	while c < N:
		i0,j0 = ind[np.argmin(d)]
		indE.append([i0,j0])
		dE.append(np.min(d))
		d   = [  d[i] for i in range(len(ind)) if (i0 != ind[i][0]) and (j0 != ind[i][1])]
		ind = [ind[i] for i in range(len(ind)) if (i0 != ind[i][0]) and (j0 != ind[i][1])]
		c = c + 1
	indE = [indE[i] for i in range(len(indE)) if dE[i]<=prop]
	dE   = [  dE[i] for i in range(len(dE))   if dE[i]<=prop]
	return [[events0[i],events1[j]] for i,j in indE]

def firstFalse(cond):
	for i in range(len(cond)):
		if cond[i] == []: continue
		for j in range(len(cond[i])):
			if not cond[i][j]:
				return i,j
	return -999,-999
	

def linkTime(vq,dates,Source,flux,min_duration,min_width,prop,dt,dcrit):
	# Links event snapshots in time
	if Source == 'ERAInt': ds = MMDataServer(Field='U',Source=Source)
	if Source != 'ERAInt': ds = cmipDataServer(Field='hus' ,LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
	Events = [getEventIndices(vq[i,:],flux=flux,min_width=min_width,dcrit=dcrit) for i in range(len(vq))]
	cond   = [ [False for i in range(len(Events[j]))] for j in range(len(Events)) ]

	# Initialise
	curr         = 0
	i0,j0        = firstFalse(cond)
	G            = [ [ Events[i0][j0]         ] ]
	Q            = [ [ vq[i0][Events[i0][j0]] ] ]
	D            = [ [ dates[i0]              ] ]
	cond[i0][j0] = True

	c = 0
	# Loop through all time and all valid section snapshots (closest sections are linked)
	while (i0,j0) != (-999,-999) and (i0 < len(vq)-1-1):
		for i in range(i0,len(Events)-1):	
			indE = linkStep(Events[i],Events[i+1],prop=prop)
			app  = False
			for j in range(0,len(indE)):
				#dhour = ds.getHours(*dates[i+1]) - ds.getHours(*dates[i])
				dhour = 24*(ds.getDays(*dates[i+1]) - ds.getDays(*dates[i]))
				if (G[curr][-1] in indE[j]) and (dhour <= dt):
					G[curr].append(indE[j][1])
					Q[curr].append(vq[i+1][indE[j][1]])
					D[curr].append(dates[i+1])
					jj0 = Events[i+1].index(indE[j][1])
					cond[i+1][jj0] = True
					app            = True
					break
			if not app:
				# Initialise new injection
				curr = curr + 1
				i0,j0 = firstFalse(cond)
				if (i0,j0) == (-999,-999): break
				G.append([])
				Q.append([])
				D.append([])
				G[curr].append(Events[i0][j0])
				Q[curr].append(vq[i0][Events[i0][j0]])
				D[curr].append(dates[i0])
				cond[i0][j0] = True
				break

	N = [len(i) for i in G]
	G = [G[i] for i in range(len(N)) if N[i] >= min_duration]
	Q = [Q[i] for i in range(len(N)) if N[i] >= min_duration]
	D = [D[i] for i in range(len(N)) if N[i] >= min_duration]
	return G,Q,D

def dist(lons1,lons2):
	# Returns distance between two lists of east-west longitudes
	H = list(set(lons1)&set(lons2))
	if H != []: return -len(H)
	d = []
	d.append( np.abs( np.abs(lons1[0]  -  lons2[0]) - np.floor(  np.abs(lons1[0] -  lons2[0])/180. )*360) )
	d.append( np.abs( np.abs(lons1[0]  - lons2[-1]) - np.floor(  np.abs(lons1[0] - lons2[-1])/180. )*360) )
	d.append( np.abs( np.abs(lons1[-1] -  lons2[0]) - np.floor( np.abs(lons1[-1] -  lons2[0])/180. )*360) )
	d.append( np.abs( np.abs(lons1[-1] - lons2[-1]) - np.floor( np.abs(lons1[-1] - lons2[-1])/180. )*360) )
	return np.min(d)

def datesAtMaxInten(G,Q,Dates):
        dates,lons = [],[]
        for i in range(len(Q)):
                ymax    = [max(Q[i][j]) for j in range(len(Q[i]))]
                yargmax = [np.argmax(Q[i][j]) for j in range(len(Q[i]))]
                k       = np.argmax(ymax)
                l       = yargmax[k]
                dates.append(Dates[i][k])
                lons.append(G[i][k][l])
        return dates,lons

def sectorMaxInten(dates,lons,LonRange=(300,100)):
	dates_,lons_ = [],[]
	for i in range(len(dates)):
		if LonRange[0] > LonRange[1]:
			if (lons[i] < LonRange[1]) or (lons[i] > LonRange[0]):
				dates_.append(dates[i])
				lons_.append(lons[i])
		else:
			if (lons[i] < LonRange[1]) and (lons[i] > LonRange[0]):
                                dates_.append(dates[i])
                                lons_.append(lons[i])
	return dates_,lons_

def makeDataset(G,D,Q,fname):
	dates,lons = datesAtMaxInten(G,Q,D)
	# Open txt file
	File = open(fname,'w')
	File.write('Intrusion|Timestep|Date|LonStart|LonEnd|MeanFlux|MaxFlux|LongitudeOfMaxFlux\n')
	File.write('#|#|(yyyy, mm, dd, hh)|degrees east|degrees east|kg/m/s|kg/m/s|degress east\n')
	File.write('\n')
	# Loop through intrusions
	for i in range(len(G)):
	        for j in range(len(G[i])):
	                date = str(D[i][j])
			File.write('%s|%s|%s|%s|%s|%s|%s|%s\n' % (i+1,j+1,date,G[i][j][0],G[i][j][-1],round(0.304*Q[i][j].mean(),1),round(0.304*Q[i][j].max(),1),G[i][j][np.argmax(Q[i][j])])) 
	File.close()


def lagged_wavenumber(Source,G,Q,D,Freq):
	# DataServer
	if Source == 'ERAInt':
		ds  = reDataServer()
		Dir = '/mnt/climstorage/cian/vq_%s' % (Freq)
	if Source != 'ERAInt':
		ds  = cmipDataServer(Field='hus' ,LevType='plev',DataFreq='day',Source=Source,LevRange=(0,1000),LatRange=(0,90),ExpType='rcp85')
		Dir = '/mnt/climstorage/cian/rcp85/%s/mon/surface/vq_%s' % (Source,Freq)
	# Wavenumber file
	wvd  = Dataset('%s/vq_1979_01.nc' % (Dir),'r')
	lats = wvd.variables['lat'][:]
	wvn  = wvd.variables['wavenumber'][:]
	wvd.close()
	years,months = range(1979,2017+1,1),range(1,12+1,1)
	Time,vq = np.zeros((0)),np.zeros((0,len(lats),len(wvn)))
	for year in years:
		for month in months:
			wvd  = Dataset('%s/vq_%s_%02d.nc' % (Dir,year,month),'r')
			Time = np.append(Time,wvd.variables['time'][:]/24.)
			vq   = np.append(vq,wvd.variables['vq'][:],axis=0)
	#vqsc0 = vq.reshape((-1,365,86,16)).mean(axis=0)
	#vqsc  = np.zeros((0,86,16))
	#for i in range(len(years)):
	#	vqsc = np.append(vqsc,vqsc0,axis=0)
	#vqa   = vq - vqsc
	#prcnt = stats.scoreatpercentile(vqa.sum(axis=-1),95,axis=0)

	Dates = [ds.getDate(day,nd=24) for day in Time]
	# Intrusions
	dates,lons = datesAtMaxInten(G,Q,D)
	dates,lons = sectorMaxInten(dates,lons,LonRange=(0,360))
	if Freq == 'dailymean':
		dt = 1
	if Freq == '6hourly':
                dt = 0.25
	days = [ds.getDays(*date) for date in dates]
	# Loop through intrusions
	lags = np.arange(-20,20+dt,dt)
	VQa  = np.zeros((0,len(lags),len(lats),len(wvn)))
	for day in days:
		x = np.argmin((Time-day)**2)
		if (Time[x] == day) and (len(Time[x-(len(lags)/2 + 1):x+len(lags)/2])==len(lags)):
			try:
				date0  = Dates[x]
				dates_ = [(y,date0[1],date0[2],date0[3]) for y in years]
				xs     = [Dates.index(date) for date in dates_]
				sclim  = np.array([vq[x_-(len(lags)/2 + 1):x_+len(lags)/2,:,:] for x_ in xs if len(Time[x_-(len(lags)/2 + 1):x_+len(lags)/2])==len(lags)]).mean(axis=0)
				VQa    = np.append(VQa,vq[np.newaxis,x-(len(lags)/2 + 1):x+len(lags)/2,:,:] - sclim,axis=0)
			except:
				pass
	# Plot
	cseq = np.arange(-4,4+0.5,0.5)

	pl.figure(1)
	flux  = VQa[:,:,:,0:15+1].sum(axis=-1)
	sx,sy = stipling(flux,xx=None,yy=None,x=lats,y=lags,thresh=0.65)
	pl.contourf(lats,lags,flux.mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
	#pl.plot(sx,sy,'k.',alpha=0.5)
	pl.colorbar()
	pl.title('%s intrusions: k = 0-15' % (len(VQa)))

        pl.figure(2)
        flux  = VQa[:,:,:,0]
        sx,sy = stipling(flux,xx=None,yy=None,x=lats,y=lags,thresh=0.65)
        pl.contourf(lats,lags,flux.mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
	#pl.plot(sx,sy,'k.',alpha=0.5)
        pl.colorbar()
        pl.title('%s intrusions: k = 0' % (len(VQa)))

        pl.figure(3)
        flux  = VQa[:,:,:,1]
        sx,sy = stipling(flux,xx=None,yy=None,x=lats,y=lags,thresh=0.65)
        pl.contourf(lats,lags,flux.mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
	#pl.plot(sx,sy,'k.',alpha=0.5)
        pl.colorbar()
        pl.title('%s intrusions: k = 1' % (len(VQa)))

        pl.figure(4)
        flux  = VQa[:,:,:,2]
        sx,sy = stipling(flux,xx=None,yy=None,x=lats,y=lags,thresh=0.65)
        pl.contourf(lats,lags,flux.mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
	#pl.plot(sx,sy,'k.',alpha=0.5)
        pl.colorbar()
        pl.title('%s intrusions: k = 2' % (len(VQa)))

        pl.figure(5)
        flux  = VQa[:,:,:,6:15+1].sum(axis=-1)
        sx,sy = stipling(flux,xx=None,yy=None,x=lats,y=lags,thresh=0.65)
        pl.contourf(lats,lags,flux.mean(axis=0),cseq,cmap=pl.cm.RdBu_r,extend='both')
        #pl.plot(sx,sy,'k.',alpha=0.5)
        pl.colorbar()
        pl.title('%s intrusions: k = 6-15' % (len(VQa)))

	pl.show()

if __name__ == "__main__": 

	# Input parameters from command line
	Source       = str(sys.argv[1])
	min_duration = int(sys.argv[2])
	min_width    = int(sys.argv[3])
	flux         = int(sys.argv[4])
	dt           = int(sys.argv[5])
	time_freq    = int(sys.argv[6])
	dcrit        = int(sys.argv[7])
	prop         = int(sys.argv[8])
	Season       = str(sys.argv[9])
	LevRange     = (int(sys.argv[10]),int(sys.argv[11]))
	YearRange    = (int(sys.argv[12]),int(sys.argv[13]))
	blat         = int(sys.argv[14])
	fblat        = int(sys.argv[15])

	# Make injection file
	fname = '/mnt/climstorage/cian/newintrusions/%s_intrusions.%s-%s.%s.%sx%shr.%sdeg.%s.%sdt.%s.%s.%sN.p' % (Source,YearRange[0],YearRange[1],Season,min_duration,time_freq,min_width,flux,dt,dcrit,prop,blat)
	if not os.path.isfile(fname):
		vq,dates = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/%s/%s/%s.moist.%s-%s.%s-%shPa.%sN.%s.p' % (blat,Source,Source,YearRange[0],YearRange[1],LevRange[0],LevRange[1],blat,Season))
		if (time_freq == 24) and (Source=='ERAInt'):
			n0,n1 = vq.shape
			vq    = vq.reshape((-1,4,n1)).mean(axis=1)
			dates = dates[0::4]
		vq,lons = interpolateFluxFile(vq)
		G,Q,D   = linkTime(vq,dates,Source=Source,flux=flux,min_duration=min_duration,min_width=min_width,prop=prop,dt=dt,dcrit=dcrit)
		toPick([G,Q,D],fname)
		print 'Made injection file %s' % (fname)
	else:
		print 'Injection file %s already exists. Did nothing.' % (fname)
		G,Q,D = unpick(fname)
	testInj(G,D)
	G_,Q_,D_ = filterInjections(fname,blat=fblat,Source=Source,blat0=blat)
