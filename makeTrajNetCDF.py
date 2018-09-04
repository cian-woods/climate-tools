from Trajectories import *

tr = Trajectories()

for year in years:
	for month in months:
		tracks = tr.getMonth(Year=year,Month=month)
