from ReanalysisDataServer import *

import matplotlib.pyplot as pl
import numpy as np

dnl = DataServer(Field='Ts',Source='CAM4xCO2')
dwl = DataServer(Field='U',Source='ERAInt')

years = range(1900,1919+1,1)
for year in years:
	dates = dnl.getDateList(Year=year,Season='Annual')
	print dates[230:240]
	hours = np.array([dnl.getHours(*i) for i in dates])
	time  = Dataset('/mnt/climstorage/cian/EXPS/Expt_4xCO2_year/h1/TS/TS_%s.nc' % (year),'r').variables['time'][:]

	print dnl.getHours(year,3,1,0),time[236],dwl.getDate(dnl.getHours(year,3,1,0)),dnl.getDate(dnl.getHours(year,3,1,0))
	print (hours==time).all(),len(time),len(hours)

	pl.plot(time,'k')
	pl.plot(hours,'r')
	pl.title(year)
	pl.show()
