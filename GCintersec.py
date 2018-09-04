import numpy as np
from math import atan2,asin

# radians per degree
rpd   = np.pi/180.
# Points on great circle 1.
glat1 = rpd*80
glon1 = rpd*0
glat2 = rpd*60
glon2 = rpd*15
# Points on great circle 2.
cglat1 = rpd*70
cglon1 = rpd*0
cglat2 = rpd*70
cglon2 = rpd*360
# 1. Put in polar coords.
x1 = np.cos(glat1) * np.cos(glon1)
y1 = np.cos(glat1) * np.sin(glon1)
z1 = np.sin(glat1)
x2 = np.cos(glat2) * np.cos(glon2)
y2 = np.cos(glat2) * np.sin(glon2)
z2 = np.sin(glat2)
cx1 = np.cos(cglat1) * np.cos(cglon1)
cy1 = np.cos(cglat1) * np.sin(cglon1)
cz1 = np.sin(cglat1)
cx2 = np.cos(cglat2) * np.cos(cglon2)
cy2 = np.cos(cglat2) * np.sin(cglon2)
cz2 = np.sin(cglat2)
# 2. Get normal to planes containing great circles.
#    It's the np.cross product of vector to each point from the origin.
N1 = np.cross([x1, y1, z1], [x2, y2, z2])
N2 = np.cross([cx1, cy1, cz1], [cx2, cy2, cz2])
# 3. Find line of intersection between two planes.
#    It is normal to the poles of each plane.
L = np.cross(N1, N2)
# 4. Find intersection points.
X1   = L / np.sqrt(L[0]**2 + L[1]**2 + L[2]**2)
X2   = -X1
ilat1 = asin(X1[2]) * 180./np.pi
ilon1 = atan2(X1[1], X1[0]) * 180./np.pi
ilat2 = asin(X2[2]) * 180./np.pi
ilon2 = atan2(X2[1], X2[0]) * 180./np.pi

print ilat1,ilon1
print ilat2,ilon2
