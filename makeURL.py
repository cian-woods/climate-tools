base = {}
base['MPI-ESM-LR'] = 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/atmos/day/r1i1p1/v20111014'
base['FGOALS-g2']  = 'http://esgf2.dkrz.de/thredds/fileServer/cmip5/output1/LASG-CESS/FGOALS-g2/rcp85/day/atmos/day/r1i1p1/v2'


Model,Field = 'FGOALS-g2','hus'
f = open('wgetURLs_1.txt','w')
for year in range(2006,2100+1,1):
	File   = '%s_day_%s_rcp85_r1i1p1_%s0101-%s1231.nc\n' % (Field,Model,year,year)
	string = '%s/%s/%s' % (base[Model],Field,File)
	f.write(string)
f.close()
