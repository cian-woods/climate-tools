import cPickle as pickle

def unpick(file):
    pkl_file=open(file,'r')
    #pkl_file=open(file,'rb')
    cycloneList=[]
    count=0
    while True:
        data=pickle.load(pkl_file)
        if (type(data)==str) and (data=='End'):
            break
        else:
            cycloneList.append(data)
        count=count+1
    return cycloneList
