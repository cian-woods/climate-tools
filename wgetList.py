import sys,os

def wgetList(wgetfile,Year0=2050,Year1=2100,exptype='rcp85'):
# Opens wget file from cmip5 site and creates new file with location URLs of all files for the requested fields for year range

	# List of variables and lower time bound for data
	fields = ['za']

	# Open wget file
	f = open(wgetfile,'r')
	l = list(f)

	# Model name
	model = l[27].split('_')[2]

	# Open new files for URLs and model names (overwritten each time)
	#base1    = '../'+exptype+'/'+model+'/'+'wgetURLs'
	#base2    = '../'+exptype+'/'+model+'/'+'models'
        base1    = '/mnt/climstorage/cian/'+exptype+'/'+model+'/'+'wgetURLs'
        base2    = '/mnt/climstorage/cian/'+exptype+'/'+model+'/'+'models'
	cc       = 0
	name1    = base1 + '_' + str(cc) + '.txt'
	name2    = base2 + '_' + str(cc) + '.txt'
	while os.path.isfile(name1) == True:
		cc    = cc + 1
		name1 = base1 + '_' + str(cc) + '.txt'
		name2 = base2 + '_' + str(cc) + '.txt'
	wgetURLs = open(name1,'w')
	models   = open(name2,'w')

	# Variables start at line 27
	c = 27
	while l[c][0] == "'":
		for ii in fields:
			# year = end year of CMIP dataset. Only want data for Year0 and after.
			year = int(l[c].split(' ')[0][1:-1].split('_')[-1].split('-')[1][0:4])
			if (l[c][1:len(ii)+2]==ii+'_') and (year>=Year0) and (year<=Year1):
				wgetURLs.write(l[c].split(' ')[1][1:-1]+'\n')
				models.write(l[c].split(' ')[0].split('_')[2]+'\n')
		c = c + 1

	wgetURLs.close()
	models.close()

	return name1,name2,cc

if __name__ == "__main__":
	
	wgetfile = str(sys.argv[1])
	try:
		Year0	 = int(sys.argv[2])
		Year1	 = int(sys.argv[3])
		exptype  = str(sys.argv[4])
	except:
		Year0   = 1979
		Year1   = 2005
		exptype = 'historical'

	n1,n2,cc = wgetList(wgetfile,Year0,Year1,exptype=exptype)

	print os.path.basename(n1)
	print os.path.basename(n2)
	print cc
