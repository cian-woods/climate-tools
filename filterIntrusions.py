from UnPickle import *
from toPick import *

import numpy as np
import os,sys

YearRange    = (int(sys.argv[1]),int(sys.argv[2]))
Season       = str(sys.argv[3])
Model        = str(sys.argv[4])
G,Q,D        = unpick('../intrusions/%s_intrusions.%s-%s.%s.6x6hr.9deg.200.6dt.20.5.p' % (Model,YearRange[0],YearRange[1],Season))
LONf,LATf,Pf = unpick('intrusiontracks/%s_fwrd_tracks.%s-%s.%s.6hour.p' % (Model,YearRange[0],YearRange[1],Season))
LONb,LATb,Pb = unpick('intrusiontracks/%s_back_tracks.%s-%s.%s.6hour.p' % (Model,YearRange[0],YearRange[1],Season))

gg,qq,dd     = [],[],[]
lonf,latf,pf = [],[],[]
lonb,latb,pb = [],[],[]
for ii in range(len(G)):
	Ntot,Nlat = 0,0
	for t in range(len(G[ii])):
		Ntot = Ntot + len(G[ii][t])
		Nlat = Nlat + len(np.where((LATf[ii][t]>=80).mean(axis=1))[0])
	if 1.*Nlat/Ntot > 0.4:
		gg.append(G[ii])
		qq.append(Q[ii])
		dd.append(D[ii])
		lonf.append(LONf[ii])
		latf.append(LATf[ii])
		pf.append(Pf[ii])
                lonb.append(LONb[ii])
                latb.append(LATb[ii])
                pb.append(Pb[ii])

print '%s of %s intrusions passed filtering' % (len(gg),len(G))
fname1 = '../intrusions/%s_intrusions.%s-%s.%s.6x6hr.9deg.200.6dt.20.5.filtered.p' % (Model,YearRange[0],YearRange[1],Season)
fname2 = 'intrusiontracks/%s_fwrd_tracks.%s-%s.%s.6hour.filtered.p' % (Model,YearRange[0],YearRange[1],Season)
fname3 = 'intrusiontracks/%s_back_tracks.%s-%s.%s.6hour.filtered.p' % (Model,YearRange[0],YearRange[1],Season)
if not os.path.isfile(fname1):
	toPick([gg,qq,dd],fname1)
else:
	print '%s already exists' % (fname1)
if not os.path.isfile(fname2):
	toPick([lonf,latf,pf],fname2)
else:
        print '%s already exists' % (fname2)
if not os.path.isfile(fname3):
        toPick([lonb,latb,pb],fname3)
else:
        print '%s already exists' % (fname3)
