from numpy import zeros
from scipy import stats

import sys

def datesToN(dates,YearRange=(1990,2012)):

	years = range(YearRange[0],YearRange[1]+1)
	n     = zeros(len(years))

	for i in dates:
		for j in range(len(years)):
			if i[0] == years[j]:
				n[j] = n[j] + 1

	return years,n

def datesToSeason(dates,YearRange,Season='DJF'):

	years = range(YearRange[0],YearRange[1]+1)
        n     = zeros(len(years))

	ind = [[] for i in years]
	d   = [[] for i in years]
	if (Season == 'DJF') or (Season == 'DJ'):
		for ii in range(len(dates)):
			date = dates[ii]
			for j in range(len(years)):
				if (date[0] == years[j]) and ((date[1] == 1) or (date[1] == 2)):
					n[j] = n[j] + 1
					d[j].append(date)
					ind[j].append(ii)
				elif (date[0] == years[j]-1) and (date[1] == 12):
					n[j] = n[j] + 1
					d[j].append(date)
					ind[j].append(ii)
	elif Season == 'ON':
                for ii in range(len(dates)):
			date = dates[ii]
                        for j in range(len(years)):
                                if (date[0] == years[j]-1) and ((date[1]) == 10 or (date[1] == 11)):
                                        n[j] = n[j] + 1
					d[j].append(date)
					ind[j].append(ii)

	elif Season == 'ONDJ':
                for ii in range(len(dates)):
			date = dates[ii]
                        for j in range(len(years)):
                                if (date[0] == years[j]) and (date[1] == 1):
                                        n[j] = n[j] + 1
					d[j].append(date)
					ind[j].append(ii)
                                elif (date[0] == years[j]-1) and (date[1] >= 10):
                                        n[j] = n[j] + 1
					d[j].append(date)
					ind[j].append(ii)

        elif Season == 'NDJFM':
                for ii in range(len(dates)):
                        date = dates[ii]
                        for j in range(len(years)):
                                if (date[0] == years[j]) and (date[1] <= 3):
                                        n[j] = n[j] + 1
                                        d[j].append(date)
                                        ind[j].append(ii)
                                elif (date[0] == years[j]-1) and (date[1] >= 10):
                                        n[j] = n[j] + 1
                                        d[j].append(date)
                                        ind[j].append(ii)

	return years,n,d,ind

if __name__ == "__main__":

	from UnPickle import *
	import matplotlib.pyplot as pl

	flux    = str(sys.argv[1])
	year1   = int(sys.argv[2])
	year2   = int(sys.argv[3])

	G,Q,D   = unpick('intrusionfiles/ERAInt_intrusions.%s.6x6hr.9deg.%s.6dt.20.5.filtered.23steps.80N.0.4.-100.0.full.p' % ('ONDJ',flux))
	dates   = [d[0] for d in D]

#	dates   = unpick('intdates1.p') + unpick('intdates2.p')

	years,n1,dys1 = datesToSeason(dates,(year1,year2),'ON')
	years,n2,dys2 = datesToSeason(dates,(year1,year2),'DJ')

	print len(dates),n1.sum(),n2.sum()

	slope, intercept, r_value, p_value, std_err = stats.linregress(years,n1)
	line1 = [slope*year + intercept for year in years]
	print slope*10,100*(slope*10./line1[0]),p_value/2.
	slope, intercept, r_value, p_value, std_err = stats.linregress(years,n2)
	line2 = [slope*year + intercept for year in years]
	print slope*10,100*(slope*10./line2[0]),p_value/2.
	slope, intercept, r_value, p_value, std_err = stats.linregress(years,n1+n2)
	line3 = line1 + line2
	print slope*10,100*(slope*10./line3[0]),p_value/2.

	PP,MM = [],[]
	pl.figure(1)
	pp, = pl.plot(years,n1,'b',linewidth=2)
	pl.plot(years,line1,'b--',linewidth=1)
	PP.append(pp)
	MM.append('ON')
	pp, = pl.plot(years,n2,'k',linewidth=2)
	pl.plot(years,line2,'k--',linewidth=1)
	PP.append(pp)
	MM.append('DJ')

	lg = pl.legend(PP,MM,loc=0,frameon=False,prop={'size':10},ncol=1,columnspacing=1,handletextpad=0.2,title=None)
	pl.setp(lg.get_title(),fontsize=10)

	pl.ylim(2,18)
	pl.xlim(year1,year2)
	pl.xlabel('Year')
	pl.ylabel('No. of Intrusions')

	pl.savefig('N.series.pdf',format='pdf')

	pl.show()
#	pl.close()




