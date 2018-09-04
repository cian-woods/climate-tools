from netCDF4 import Dataset
from LambertProjector import *
from drawCoastlinesNoRivers import *
import matplotlib.pyplot as pl

Projection = 'npstere'
proj = LambertProjector(boundinglat=20,resolution=300.,Projection=Projection)
print np.diff(proj.x[0,:])/1000.
pl.plot(proj.x,proj.y,'k.',markersize=2.5)
drawCoastlinesNoRivers(proj.m)
pl.title('resolution = %s km' % (proj.res))
pl.savefig('../WarmArctic/figs/%s_350x350km.pdf' % (Projection),format='pdf')
pl.show()

"""
fname = 'nplaea_axes.nc'
if not os.path.isfile(fname):
	# Create file
	File = Dataset(fname,'w',format='NETCDF3_CLASSIC')
	# Define some global attribs
	File.Conventions='COARDS'
	# Horizontal axes
	File.createDimension('X',proj.nx)
	var           = File.createVariable('X','f',('X',))
	var.long_name = 'Distance along x-axis'
	var.units     = 'Meters'
	var[:]        = proj.x[0,:].astype('f')
	File.createDimension('Y',proj.ny)
	var           = File.createVariable('Y','f',('Y',))
	var.long_name = 'Distance along y-axis'
	var.units     = 'Meters'
	var[:]        = proj.y[:,0].astype('f')
	# Create Variables
	var           = File.createVariable('lat','f',('X','Y'))
	var.long_name = 'Latitude'
	var.units     = 'Degrees North'
	var[:]        = proj.lat.astype('f')
	var           = File.createVariable('lon','f',('X','Y'))
	var.long_name = 'Longitude'
	var.units     = 'Degrees East'
	var[:]        = proj.lon.astype('f')
	# Close file
	File.close()
"""
