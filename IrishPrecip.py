import numpy as np
import matplotlib.pyplot as pl

SMs = {}
SMs['JJA']    = [6,7,8]
SMs['DJF']    = [12,1,2]
SMs['Annual'] = range(1,12+1,1)

class IrishPrecip:

	def __init__(	self	):

		self.File   = open('IOI_1711_SERIES.CSV','r')
		self.List   = list(self.File)
		self.prec   = np.array([float(self.List[i].split(',')[-1][0:-2]) for i in range(1,len(self.List))])
		self.years  = range(1711,2016+1,1)
		self.annual = range(1,12+1,1)
		self.months = [self.annual[i%12] for i in range(len(self.prec))]

	def getSeason(self,Year=1712,Season='JJA'):
		if SMs[Season][0] < SMs[Season][-1]:
			ix0 = 0 + self.years.index(Year)*12
			f   = self.prec[ix0+SMs[Season][0]-1:ix0+SMs[Season][0]-1+len(SMs[Season])].sum()
		else:
			ix0    = 0 + self.years.index(Year-1)*12
			prec   = self.prec[ix0:ix0+24]
			f      = prec[SMs[Season][0]-1:SMs[Season][0]-1+len(SMs[Season])].sum()
		return f

if __name__ == "__main__":

	ip = IrishPrecip()

	years = range(1980,2016+1,1)	
	p     = [ip.getSeason(Year=year,Season='JJA') for year in years]

	pl.plot(years,p)
	pl.show()
