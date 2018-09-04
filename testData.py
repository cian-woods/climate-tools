from cmipDataServer import *

#Models = [g[9:] for g in glob.glob('../rcp85/*')]
Models = ['MRI-CGCM3']
years  = range(1980,2100+1,1)

for Model in Models:

	d = DataServer(Field='va',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='day',LatRange=(70,70),LonRange=(0,0),LevRange=(500,500))
	for year in years:
		print 'va',year
		s = d.getDataSnaps(Year=year,Season='Annual')
	d.closeFiles(d.Handles)

        d = DataServer(Field='hus',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='day',LatRange=(70,70),LonRange=(0,0),LevRange=(500,500))
        for year in years:
                print 'hus',year
                s = d.getDataSnaps(Year=year,Season='Annual')
        d.closeFiles(d.Handles)

        #d = DataServer(Field='ta',LevType='plev',Source=Model,ExpType='rcp85',DataFreq='day',LatRange=(70,70),LonRange=(0,0),LevRange=(500,500))
        #for year in years:
        #        print 'ta',year
        #        s = d.getDataSnaps(Year=year,Season='Annual')
        #d.closeFiles(d.Handles)

