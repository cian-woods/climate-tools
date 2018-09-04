import numpy as np
import sys

dT,Phi0,dPhi = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])

omega = 7.29e-5
f     = round(2*omega*np.sin( ( (Phi0+(dPhi/2) )*np.pi)/180. ),6)
R     = 287
g     = 9.8
T     = 250
H     = R*T/g

dTdy  = round((dT/dPhi)*(1/(111.*1000)),10)
dUdz  = round(-1000.*(R/(f*H))*dTdy,3)	# m s**-1 km**-1

print 'Between %sN and %sN, with dT = %s K; f(phi=%s) = %s s**-1\ndTdy = %s K Mm**-1 and dU/dz = %s m s**-1 km**-1' % (Phi0,Phi0+dPhi,dT,Phi0+dPhi/2,f,1e06*dTdy,dUdz)
