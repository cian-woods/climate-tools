import sys,os
from numpy import *
from sets import *
from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from scipy import interpolate

def vqEvents(	vq,
		dates,
		Server  =  'CMIP',
		model	=  'MPI-ESM-LR',
		ExpType =  'historical',
		dur     =   2,
		box     =   11.25,
		flux    =   200,
		del_t   =   6,
		dist1   =   10,
		dist2   =   15	):

	# Takes in vertically integrated vq data (time)x(longitude) and search parameters

	# dur = min no. of timesteps in  event
	# box = min longitude span of event
	# flux (Mt/day) = min flux in  in each snapshot of event; del_t = min #snapshots between events
	# del_t = max time allowed between event snapshots in hours
	# one timesptep max in between etc...
		# dist1 = max degrees of longitude "same event" can be seperated in a  snapshot.
		# dist2 = max distance an event can move between snapshots

		# Goes through each longitude (wrap around at 0/359): If vq > flux then begin counting.
		# If 'vq' stays >= 'flux' for  at least 'dur' timesteps and 'box' degrees of longitude
		# then it is an event.

	# Returns timeseries of event snapshots. Secondary script (eg. dur = 3) will 'stitch' events together between snapshots
	# N.B: Flux area can be moving between snapshots. Most likely moving eastward.

	# Start search from argmin of each snapshot. Prevents starting in the middle of an event

	vq = array(vq)
#       tmp = '''
        print 'Interpolating data to 1 degree resolution...'
        # Add wrap value to vq
        xn   = len(vq[0])
        ends = vq[:,-1]
        vq   = append(vq,ends[:,newaxis],1)

        # Interpolate vq to 1 degree resolution
        yold = arange(len(dates))
        ynew = arange(len(dates))
        xold = array(list(arange(0,360,360./xn))+[360])
        xnew = arange(0,360,1.)
        print vq.shape
	print len(xold)
	print len(yold)
        f   = interpolate.interp2d(xold,yold,vq,kind='linear',bounds_error=True)
	vq  = f(xnew,ynew)
	print vq.shape

	if Server == 'ERAInt':
		tt = reDataServer(Field='T2',LevType='surface_analysis',Source='ERAInt')
	elif Server == 'CMIP':
		tt = cmipDataServer(Field='psl',LevType='surface',DataFreq='day',Source=model,ExpType=ExpType)
	elif Server == 'NCEP1':
		tt = reDataServer(Field='shum',Source='NCEP1')
        elif Server[0:3] == 'CAM':
                tt = MMDataServer(Field='Ts',Source=Server)
	
	lonlen = len(vq[0])
	res    = 360./lonlen

	print 'Resolution of data is '+str(res)+' degrees'
#	var   = input('Enter minimum number of gridpoints across width of intrusion: ')
	width = floor(1.*box/res)
	print 'Minimum width of intrusion is ' + str(width*res) + ' degrees'

	H   = []
	G   = []
	V   = []
	ind = []
	c1 = 0

	# Find events exceeding flux over longitude span > box in each snapshot
	for i in range(len(vq)):
		x = argmin(vq[i])
		H.append([])
		V.append([])
		c = 0
		while c < len(vq[i]):
			a = 0
#			f = 0
			f = []
			h = []
			y = x + c
			y = int(y-(lonlen*floor(1.*y/lonlen)))
			while vq[i][y] >= flux:
#				f = f + vq[i][y]
				f.append(vq[i][y])
				a = 1
#				h.append(tt.lon[y])
				h.append(xnew[y])
				c = c + 1
				y = x + c
				y = int(y-(lonlen*floor(1.*y/lonlen)))
#				print i,y,c
#				print vq[i][y]
			if len(h) >= width:
				H[c1].append(h)
				V[c1].append(f)
				ind.append(i)
			c = c + 1 - a
		c1 = c1 + 1


	
	# Create W, which is the same as H but events close enough are joined i.e. <= dist
	if dist1 != None:
		WH    = []
		WV    = []
		ind_w = []
		for i in range(len(H)):
			# D list is for distances between events in same snapshot
			D = []
			if len(H[i]) == 1:
				WH.append([H[i][0]])
				WV.append([V[i][0]])
				ind_w.append(i)
			elif len(H[i]) >= 2:
				# Use i_x to get anti-clockwise distances between all events 
				i_x = range(len(H[i]))+[0]
				for j in range(len(H[i])):
					i1,i2 = i_x[j],i_x[j+1]
					d1,d2 = H[i][i2][0],H[i][i1][-1]
					if d1 < d2:
						d = 360 + (d1-d2)
					else:
						d = d1-d2
					D.append(d)

				# Joing events from here on
				w = [H[i][0][:]]
				v = [V[i][0][:]]
				c = 0
				for j in range(len(D)):
					l = j + 1
					if l > len(H[i])-1:
						l = 0
					if (D[j] <= dist1) and (l != 0):
						w[c] = w[c] + H[i][l][:]
						v[c] = v[c] + V[i][l][:]
					elif (l != 0):
						w.append([])
						v.append([])
						c = c + 1
						w[c] = w[c] + H[i][l][:]
						v[c] = v[c] + V[i][l][:]
					elif (D[j] <= dist1) and (l == 0):
						w[0] = w[0] + w[c]
						v[0] = v[0] + v[c]
				WH.append(w)
				WV.append(v)
				ind_w.append(i)
			else:
				WH.append([])
				WV.append([])

		H   = WH[:]
		V   = WV[:]
		ind = ind_w[:]


	# Stitch snapshots of same events together
	if dur != None:
		G = []
		Q = []
		J = []
		S = []
		D = []
		M = range(len(H[0]))
		c = 0
		while c < len(ind)-1:
			x,y = ind[c],ind[c+1]
			n   = x
			h   = []
			k   = []
			C   = []
			for j in M:
#			for j in range(len(H[n])):
				c1  = c
				ic  = j
				P   = [H[n][j]]
				PV  = [V[n][j]]
				x,y = ind[c],ind[c+1]
				s   = [x]
				d   = [dates[x]]
				while (tt.getHours(*dates[y]) - tt.getHours(*dates[x]) <= del_t) and (c1 < len(ind)-1):
					l  = dist([H[x][ic]],H[y],dist2)
					ic = l[0]
					if ic != None:
						P.append(H[y][ic])
						PV.append(V[y][ic])
						s.append(y)
						d.append(dates[y])
					else:
						break
					if len(H[y])>len(H[x]):
						s1 = set(dist(H[x],H[y],dist2))
						s2 = set(range(len(H[y])))
						h.append(c1+1)
						k.append(list(s2-s1))
					C.append(c1+1)
					c1 = c1 + 1
					try:
						x,y = ind[c1],ind[c1+1]
					except:
						break
				if len(P) >= dur:
					G.append(P)
					Q.append(PV)
					J.append(ind[c])
					S.append(s)
					D.append(d)
					
			if C != []:
				c = max(C)
				try:
					M = range(len(H[ind[c]]))
				except:
					break
			if h != []:
				o = argmin(h)
				c = h[o]
				M = k[o]
			else:
				c = c + 1
				try:
					M = range(len(H[ind[c]]))
				except:
					break

		# Get rid of any potential double events in data
		Y = []
		X = []
		R = []
		L = []
		E = []
		for i in range(len(G)-1):
			t1 = set(range(J[i],J[i]+len(G[i]),1))
			t2 = set(range(J[i+1],J[i+1]+len(G[i+1]),1))
			t  = sort(list(t1&t2))
			if t != []:
				t1,t2 = sort(list(t1)),sort(list(t2))
				x1 = argmin((t1-t[0])**2)
				y1 = argmin((t1-t[-1])**2)
				x2 = argmin((t2-t[0])**2)
				y2 = argmin((t2-t[-1])**2)
				r1 = G[i][x1:y1+1]
				r2 = G[i+1][x2:y2+1]
				if r1 == r2:
					R.append(i+1)
		R = set(R)
		for i in range(len(G)):
			if list(set([i])&R) == []:
				Y.append(G[i])
				X.append(Q[i])
				L.append(S[i])
				E.append(D[i])
		G = Y[:]
		Q = X[:]
		S = L[:]
		D = E[:]
				

	return G,Q,D,vq
#	return H,V,ind,G,J,S,Q,D


def dist(H1,H2,dist2):
	# H1 and H2 are lists of events
	l = []
	for i in range(len(H1)):
		D = []
		J = []
		for j in range(len(H2)):
			d1 = abs(H1[i][0]-H2[j][0])
			d2 = abs(H1[i][-1]-H2[j][-1])
			d1 = abs(d1-floor(d1/180.)*360)
			d2 = abs(d2-floor(d2/180.)*360)
			d  = min([d1,d2])
			D.append(d)
			J.append(j)
		x = argmin(D)
		if D[x] <= dist2:
			l.append(J[x])
		else:
			l.append(None)
	return l
			
	

if __name__ == "__main__":
	from UnPickle import *
	from toPick import *
	from scipy.stats import stats as stats

	model     = str(sys.argv[1])
	exptype   = str(sys.argv[2])
	dur       = int(sys.argv[3])
	box       = int(sys.argv[4])
	flux      = int(sys.argv[5])
	del_t     = int(sys.argv[6])
	dist1     = int(sys.argv[7])
	dist2     = int(sys.argv[8])	
	Season    = str(sys.argv[9])
	YearRange = (int(sys.argv[10]),int(sys.argv[11]))

	# Flux file
	if model == 'ERAInt':
		p      = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/70/ERAInt/ERAInt.moist.%s-%s.30-1000hPa.70N.%s.p' % (YearRange[0],YearRange[1],Season))
		Server = 'ERAInt'
        elif model == 'NCEP1':
                p      = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/70/NCEP1/NCEP1.moist.%s-%s.300-1000hPa.70N.%s.p' % (1950,2016,Season))
                Server = 'NCEP1'
	elif model[0:3] == 'CAM':
		p      = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/70/%s/%s.moist.%s-%s.0-1000hPa.70N.%s.p' % (model,model,YearRange[0],YearRange[1],Season))
		Server = model
	else:
		p      = unpick('/mnt/climstorage/cian/scripts/newfluxfiles/70/%s/%s.moist.2075-2100.30-1000hPa.70N.%s.p' % (model,model,Season))
		Server = 'CMIP'

	# Filename
#	filename = '/mnt/climstorage/cian/intrusions/%s_intrusions.%s.%sx%shr.%sdeg.%s.%sdt.%s.%s.p' % (model,Season,dur,del_t,box,flux,del_t,dist1,dist2)	
	filename = '/mnt/climstorage/cian/intrusions/%s_intrusions.%s-%s.%s.%sx6hr.%sdeg.%s.%sdt.%s.%s.p' % (model,YearRange[0],YearRange[1],Season,dur,box,flux,del_t,dist1,dist2)
	print filename

	# Intrusions
	G,Q,D,vq = vqEvents(p[0],p[1],Server=Server,model=model,ExpType=exptype,dur=dur,box=box,flux=flux,del_t=del_t,dist1=dist1,dist2=dist2)
	print len(G)

	# Stats
        l = [j for i in vq for j in i if j > 0]
        print stats.scoreatpercentile(l,90)
	q = 0
	for i in Q:
		for j in i:
			q = q + sum(j)
	print 100.*q/sum(l)


	# Save file
	if os.path.isfile(filename)==False:
		print 'Saving intrusions...'
		toPick([G,Q,D],filename)
	else:
		print 'intrusion file exists already'
	print '\n\n'
