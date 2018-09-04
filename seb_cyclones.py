import pickle

fname = '/mnt/climstorage/sebastian/tracks/einterim/einterim_tracks_19790101_20161231[-180, 180, 80, 20]_formatted.pkl'
sname = '/mnt/climstorage/cian/cyclones/ERAinterim_1.5_NH_MSS_19790101_20161231_ST_MTEX.txt'

f = open(sname,'w')
f.write('99 Code,CycloneNo,StepNo,DateI10,Year,Month,Day,Time,LongE,LatN,SLP_Central_Pressure,Intensity2,Intensity3\n')

print('Opening track_data file %s ...' % (fname))
track_data = pickle.load(open(fname,'rb'))
datestamps = list(track_data.keys())[:]

n = 1
for datestamp in datestamps:
	keys = list(track_data[datestamp].keys())
	for key in keys:
		lon  = track_data[datestamp][key].lon.values
		lat  = track_data[datestamp][key].lat.values
		msl  = track_data[datestamp][key].msl.values
		time = track_data[datestamp][key].time.values
		line = '90 %06d %03d\n' % (n,len(lon))
		f.write(line)
		for i in range(len(lon)):
			date                = str(time[i])[0:16]
			year,month,day,hour = date[0:4],date[5:7],date[8:10],date[11:13]
			yearmonthdayhour    = date[0:4]+date[5:7]+date[8:10]+date[11:13]
			line = 'SB %06d %03d %s %s %s %s %s %s %s %s -999.999 -999.999\n' % (n,i+1,yearmonthdayhour,year,month,day,hour,round(lon[i],2),round(lat[i],2),round(msl[i],1))
			f.write(line)
		n = n + 1
f.close()
