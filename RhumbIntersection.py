from mpl_toolkits.basemap import Basemap as Basemap
import matplotlib.pyplot as pl

def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False

m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')

lat11,lon11 = 70,-180
lat12,lon12 = 70,180
lat21,lon21 = 75,50
lat22,lon22 = 65,-50

j11,k11 = m(lon11,lat11)
j12,k12 = m(lon12,lat12)
j21,k21 = m(lon21,lat21)
j22,k22 = m(lon22,lat22)

L1,L2 = line([j11,k11],[j12,k12]),line([j21,k21],[j22,k22])
R     = intersection(L1, L2)
if R:
	x,y     = R
	lon,lat = m(*R,inverse=True)
	print "Intersection detected: (%s, %s)" % (lon,lat)
	m.plot([j11,j12],[k11,k12],'k')
	m.plot([j21,j22],[k21,k22],'r')
	m.plot(x,y,'b.')
	m.drawcoastlines()
	pl.show()
else:
	print "No single intersection point detected"
