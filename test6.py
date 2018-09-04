from LambertProjector import *
from ReanalysisDataServer import *

d    = DataServer()
proj = LambertProjector(boundinglat=50,resolution=100.)

s = d.snapshot(1998,1,1,0)
print s.shape
s = proj(s,d.lon,d.lat)
print s.shape

print proj.lat
print proj.lon

print proj.x
print proj.y
