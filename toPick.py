import cPickle as pk
def toPick(data,filename):
	for i in data:
		pk.dump(i,open(filename,'a'))
	pk.dump('End',open(filename,'a'))
