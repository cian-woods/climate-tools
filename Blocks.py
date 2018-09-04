import glob,sys
import numpy as np
import matplotlib.pyplot as pl

from drawCoastlinesNoRivers import drawCoastlinesNoRivers

sys.path.insert(0, '../scripts')
from Cyclones import CycloneServer
from LambertProjector import *
from ReanalysisDataServer import DataServer as reDataServer
from ReanalysisMonthlyMeanDataServer import DataServer as MMDataServer
from UnPickle import *
from toPick import *
from scipy import stats

class Blocks:

	def __init__(	self	):
		# Projector and Dataset
		self.proj = LambertProjector(boundinglat=45,resolution=600.)
		self.x    = self.proj.x[0,:]
		self.y    = self.proj.y[:,0]
		self.ds   = reDataServer(Field='T2',LevType='surface_analysis')
                # Cyclones
                self.CS = CycloneServer(method='M13',hemisphere='NH')
		# Files
		self.files = glob.glob('/qnap/cian/IMILAST/blockings/*.asc')
		self.lats  = [float(f[29:33]) for f in self.files]
		self.nlats = len(self.lats)
		self.Files = sorted(zip(self.lats,self.files))
		self.lats  = np.array([File[0] for File in self.Files])
		self.Files = [File[1] for File in self.Files]
		self.lons  = np.arange(0,360,360/320.)
		self.nlons = len(self.lons)
		self.lonsg,self.latsg = np.meshgrid(self.lons,self.lats)
		self.lonsz,self.latsz = self.lonsg.reshape(-1),self.latsg.reshape(-1)
		# Dates
		self.dates = []
		for i in range(1979,2015+1,1):
			self.dates = self.dates + self.ds.getDateList(Year=i,Season='Annual')
		self.dates = self.dates + self.ds.getDateList(Year=2016,Season='JFMA')
		self.dates = self.dates[::4]
		# Make blocking array
		self.T = []
		for File in self.Files:
			print File
			l,t = list(open(File)),[]
			for i in l:
				t.append([int(j) for j in i[1:-1].split(' ')])
			self.T.append(t)
		self.T     = np.array(self.T)
		self.T     = np.rollaxis(self.T,1,0)
		# Extract Season Clim
		self.F = np.empty((0,self.nlats,self.nlons))
		for y in range(1980,2016+1,1):
			datelist = self.ds.getDateList(Year=y,Season='NDJFM')
			year1,month1,day1,hour1 = datelist[0]
			year2,month2,day2,hour2 = datelist[-1]
			x,y    = self.dates.index((year1,month1,day1,0)),self.dates.index((year2,month2,day2,0))
			self.F = np.append(self.F,self.T[x:y+1,:,:],axis=0)
		self.Tclim    = self.F.mean(axis=0)
		self.Cclim    = np.array([self.blockDensity(self.F[i]) for i in range(len(self.F))]).mean(axis=0)	
		# Plot clim
		cseq   = np.arange(0,3.5+0.25,0.25)
		cc,x,y = self.CS.interp2d(self.Cclim,self.x,self.y,6,kind='linear')
		cf     = pl.contourf(x,y,cc,cseq,cmap=pl.cm.OrRd,extend='max')
		cbar   = pl.colorbar(cf)
		cbar.set_label('Block frequency')
		self.proj.m.drawparallels([45,50,55,60,65,70],latmax=90)
		drawCoastlinesNoRivers(self.proj.m)
		pl.show()

	def blockDensity(self,T):
		index     = T.reshape(-1)
		lats,lons = self.latsz[np.where(index==1)],self.lonsz[np.where(index==1)]
		xs,ys     = self.proj.m(lons,lats)
                Count0    = np.zeros((self.proj.nx,self.proj.ny))
                for i in range(len(xs)):
                        xi,yi    = np.argmin(abs(self.x-xs[i])),np.argmin(abs(self.y-ys[i]))
                        Count0[yi,xi] = Count0[yi,xi] + 1
		return Count0

	def getDataDateList(self,datelist):
		inds = [self.dates.index(date) for date in datelist]
		T    = self.T[inds,:,:]
		C    = [self.blockDensity(T[i]) for i in range(len(T))]
		return C

	def laggedBlocks(self):
		dates = unpick('dates.warm.50.ERAInt.p')
		days  = np.arange(-10,0+1,1)
		C     = []
		for date0 in dates:
			if (date0[0],date0[1],date0[2],0) in self.dates:
				hour0    = self.ds.getHours(*date0) - date0[3]	# Make hour0 of date0 = 0 hrs, same as blocking dates
				datelist = [self.ds.getDate(hour0 + ii*24) for ii in days]
				C.append(self.getDataDateList(datelist))
		C = np.array(C).mean(axis=0) - self.Cclim[np.newaxis,:,:]
		cseq = np.arange(-2,2+0.25,0.25)
		for ii in range(len(C)):
			cc,x,y = self.CS.interp2d(C[ii],self.x,self.y,6,kind='linear')
			cf     =  pl.contourf(x,y,cc,cseq,cmap=pl.cm.RdBu_r,extend='both')
                	cbar   =  pl.colorbar(cf)
                	cbar.set_label('Block frequency anomaly')
                	self.proj.m.drawparallels([45,50,55,60,65,70],latmax=90)
                	drawCoastlinesNoRivers(self.proj.m)
			pl.title('day %s' % (days[ii]))
                	pl.show()
	
	def lonsMedian(self,lons,LonRange=(330,105)):
		if LonRange[0] < LonRange[1]: lon0 = stats.scoreatpercentile(lons,50)
		if LonRange[0] > LonRange[1]:
			LON = []
			for lon in lons:
				if lon > 180:
					LON.append(lon-360)
				else:
					LON.append(lon)
			lon0 = stats.scoreatpercentile(LON,50)
			if lon0 < 0: lon0 = lon0 + 360
		return lon0

	def slpSnaps(self,date0,N=10):
		days  = np.arange(-N,N+1,1)
		dates = [self.ds1.getDate(self.ds1.getHours(*date0) + ii*24) for ii in days]
		slp   = np.array([self.ds1.snapshot(*date).squeeze() for date in dates])
#		slp   = slp - slp.mean(axis=2)[:,:,np.newaxis]
		return slp

        def datesAtMaxInten(self,GQD=None):
		if GQD != None:
			G,Q,D = GQD
		else:
			G,Q,D = self.G,self.Q,self.D
		dates,lons = [],[]
		for i in range(len(Q)):
			ymax    = [max(Q[i][j]) for j in range(len(Q[i]))]
			yargmax = [np.argmax(Q[i][j]) for j in range(len(Q[i]))]
			k       = np.argmax(ymax)
			l       = yargmax[k]
			dates.append(D[i][k])
			lons.append(G[i][k][l])
		return lons,dates

	def intrusionsWithBlocks(self):
		lons,dates  = self.datesAtMaxInten(GQD=None)
		gg1,qq1,dd1 = [],[],[]
		gg2,qq2,dd2 = [],[],[]
		for i in range(len(lons)):
			year,month,day,hour = dates[i]
			x     = np.argmin((self.lons - lons[i])**2)
			y     = self.dates.index((year,month,day,0))
			if self.nlons - x >= 20: tsnap = self.T[y,:,x:x+20]
			if self.nlons - x <  20: tsnap = np.append(self.T[y,2:,x:],self.T[y,2:,0:20-self.nlons+x],axis=1)
			if tsnap.mean() > 0:
				gg1.append(self.G[i])
				qq1.append(self.Q[i])
				dd1.append(self.D[i])
			else:
                                gg2.append(self.G[i])
                                qq2.append(self.Q[i])
                                dd2.append(self.D[i])
		return gg1,qq1,dd1,gg2,qq2,dd2

	def plot(self,index,T,title='jim',glons=None,slp=None,lonm=None,uv=None):
#		cseq1 = np.arange(-0.16,0.16+0.02,0.02)
		cseq1 = np.arange(0,0.5+0.05,0.05)
#		cseq1 = np.arange(0,0.2+0.01,0.01)
#		cseq2 = np.arange(1000,1030+2,2)
#		cseq2 = np.arange(-1200,1200+200,200)
		cseq2 = 11
		j,k   = self.proj.m(*np.meshgrid(self.lons,self.lats))
		self.proj.m.contourf(j,k,T,cseq1,cmap=pl.cm.OrRd,extend='both')
		pl.colorbar()
		if glons != None:
			j,k = self.proj.m(glons,[70 for i in glons])
			pl.plot(j,k,'r',linewidth=3)
		if slp != None:
			j,k = self.proj.m(*np.meshgrid(self.ds1.lon,self.ds1.lat))
			cl  = self.proj.m.contour(j,k,slp,cseq2,colors='0.5',linewidths=0.5)
#			cl  = pl.contour(self.proj.x,self.proj.y,slp,cseq2,colors='0.5',linewidths=0.5)
			pl.clabel(cl,fmt='%4.0f',colors='k',fontsize=7)
		if lonm != None:
			j,k = self.proj.m([lonm for i in range(10)],np.linspace(60,80,10))
			self.proj.m.plot(j,k,'k',linewidth=2)
		if uv != None:
			Q  = pl.quiver(self.proj.x[::3,::3],self.proj.y[::3,::3],uv[0][::3,::3],uv[1][::3,::3],units='inches',scale=30,\
                                scale_units='inches',headwidth=3,headlength=5,headaxislength=4.5,pivot='tail')
                        qk = pl.quiverkey(Q, 0.2, 1.02, 15, '%s%s' % (15,'m s$^{-1}$'), labelpos='W',fontproperties={'weight': 'bold'})
		self.proj.m.drawcoastlines()
		self.proj.m.drawparallels([60,70,80],latmax=90)
		pl.title(title)
		pl.savefig('blockcomps/yesblock/%s.pdf' % (index))
		pl.show()
#		pl.close()

	def compInts(self):
#		lon0,dates     = [i[0][-1] for i in self.G],[i[0] for i in self.D]
		G,Q,D          = self.G,self.Q,self.D
		G,Q,D,G1,Q1,D1 = self.intrusionsWithBlocks()
		q  = [sum([sum(i) for i in Q[j]])  for j in range(len(Q))]
		q1 = [sum([sum(i) for i in Q1[j]]) for j in range(len(Q1))]
		print np.mean(q),np.mean(q1)
		lon0,dates     = self.datesAtMaxInten(GQD=[G1,Q1,D1])
#		lon0,dates     = self.filterInjections(lon0,dates,Sector=(105,280))
#		lonm           = self.lonsMedian(lon0,LonRange=(0,360))
#		lon0,dates     = self.datesAtMaxInten(GQD=None)
		data,SLP       = [],[]
		days = np.arange(-10,10+1,1)
		for i in range(len(lon0)):
			x1 = np.argmin((self.lons - lon0[i])**2)
			x2 = np.argmin((self.ds1.lon - lon0[i])**2)
			print '%s of %s ... (%s %s)' % (i+1,len(lon0),self.lons[x1],dates[i])
			year,month,day,hour = dates[i]
			y   = self.dates.index((year,month,day,0))
#			slp = self.ds1.snapshot(*dates[i])
#			slp = slp - slp.mean(axis=1)[:,np.newaxis]
#			t   = self.T[y,:,:] # - self.Tclim[:,:]
#			self.plot(1,t,'jim',glons=[lon0[i]+ij for ij in range(10)],slp=slp)
#			d,l = shiftgrid(self.lons[x1],t,self.lons,start=True,cyclic=360.)
#			self.plot(d,glons=self.G[i][0])
			slp = self.slpSnaps(date0=(year,month,day,0),N=10)
#			slp = self.wave1(slp,nlon=self.ds1.nlon)
			t   = self.T[y-10:y+10+1,:,:] - self.Tclim[np.newaxis,:,:]
			slp = np.array([shiftgrid(self.ds1.lon[x2],sp,self.ds1.lon,start=True,cyclic=360.)[0] for sp in slp])
			t   = np.array([shiftgrid(self.lons[x1],ff,self.lons,start=True,cyclic=360.)[0] for ff in t])
#			slp = np.array([shiftgrid(360-lonm,sp,self.ds1.lon,start=True,cyclic=360.)[0] for sp in slp])
#			t   = np.array([shiftgrid(360-lonm,ff,self.lons,start=True,cyclic=360.)[0] for ff in t])
			data.append(t)
			SLP.append(slp)
		data = np.array(data).mean(axis=0)
		SLP  = np.array(SLP).mean(axis=0)
		for l in range(len(data)):
			self.plot(index=l+1,T=data[l],title=days[l],slp=SLP[l],lonm=None)

	def wave1(self,z,nlon):
	        # Takes in (time)x(lat)x(lon) v'T' and returns the wave1 pattern

	        # Divide by n becasue discrete transform and extract wave1 (241 index after shift)
	        Fk = np.fft.fft(z,axis=2)/nlon
	        Fk = np.fft.fftshift(Fk,axes=(2,))[:,:,nlon/2+1]
	        # Calcutae amplitide and phase of wave1
	        amp = 2*np.abs(Fk)
	        ang = np.arctan2(Fk.imag,Fk.real)
	        # Create wave pattern
	        S = np.array([-1*amp*np.cos(i + ang) for i in np.linspace(-np.pi,np.pi,nlon)])
	        S = np.rollaxis(S,0,3)
	        return S

	def filterInjections(self,lon0,dates,Sector=(330,105)):
		dd,ll = [],[]
		for i in range(len(lon0)):
			if Sector[0] > Sector[1]:
				if (Sector[0] <= lon0[i] <= 360) or (0 <= lon0[i] < Sector[1]):
					dd.append(dates[i])
					ll.append(lon0[i])
			elif Sector[0] < Sector[1]:
				if (Sector[0] <= lon0[i] < Sector[1]):
					dd.append(dates[i])
					ll.append(lon0[i])
                return ll,dd

	def compBlocks(self):
		dates = unpick('datesATL.p')
		data  = []
		days  = np.arange(-10,10+1,1)
		for date in dates:
			datelist = [self.dsu.getDate(self.dsu.getHours(*date) + ii*24) for ii in np.arange(0,2+0.25,0.25)]
			year,month,day,hour = date
			y   = self.dates.index((year,month,day,0))
			U   = np.array([self.proj(self.dsu.snapshot(*datei).squeeze(),self.dsu.lon,self.dsu.lat) for datei in datelist]).mean(axis=0)
			V   = np.array([self.proj(self.dsv.snapshot(*datei).squeeze(),self.dsv.lon,self.dsv.lat) for datei in datelist]).mean(axis=0)
			urot,vrot = self.proj.m.rotate_vector(U,V,self.proj.lon,self.proj.lat,returnxy=False)
			slp = self.proj(self.ds1.snapshot(*date).squeeze(),self.ds1.lon,self.ds1.lat)
			slp = slp - slp.mean(axis=1)[:,np.newaxis]
			self.plot(1,self.T[y-4:y,:,:].mean(axis=0),'jim',slp=slp,uv=[urot,vrot])
			t   = self.T[y-10:y+10+1,:,:] - self.Tclim[np.newaxis,:,:]
			data.append(t)
		data = np.array(data).mean(axis=0)
		for l in range(len(data)):
			self.plot(index=l+1,T=data[l],title=days[l],slp=None,lonm=None)

if __name__ == "__main__":

	b = Blocks()
	b.laggedBlocks()
#	b.compBlocks()
#	b.compInts()
