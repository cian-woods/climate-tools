f=open('AssociatedIntrusions.txt','r')
g=list(f)
Rank,Intrusion,Timestep,Date,Lon0,Lon1 = {},{},{},{},{},{}
					 
for i in range(1,len(g)):
#	Rank,Intrusion,Timestep,Date,Lon0,Lon1 = [],[],[],[],[],[]
	rank,intrusion,timestep,date,lon0,lon1 = g[i].split('|')
	rank,intrusion,timestep,date,lon0,lon1 = int(rank),int(intrusion),int(timestep),date,float(lon0),float(lon1[:-1])
	Rank[r].append(rank)
	Intrusion[r].append(intrusion)
	Timestep[r].append(timestep)
	Date[r].append(date)
	Lon0[r].append(lon0)
	Lon1[r].append(lon1)
print Lon1
