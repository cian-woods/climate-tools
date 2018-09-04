import cPickle as pk
def toPick1(data,filename):
	f = open(filename,'a')
	for i in data:
		pk.dump(i,f)
	pk.dump('End',f)
