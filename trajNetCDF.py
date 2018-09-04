from sys import *
from netCDF4 import *
import numpy as np

def trajNetCDF(time,lons,lats,Ps,PVs,steps,lon,FileName):

	ntime,nlon,nsteps = len(time),len(lon),len(steps)

        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'Time initialisation'
        var.units     = 'hours since 1900-01-01 00:00:0.0'
        var[:]        = time
        # Create axes
        File.createDimension('lon',nlon)
        var           = File.createVariable('lon','f',('lon',))
        var.long_name = 'Longitude initialisation'
        var.units     = 'degrees_east'
        var[:]        = lon 
        File.createDimension('steps',nsteps)
        var           = File.createVariable('steps','d',('steps',))
        var.long_name = 'Time along trajectory'
        var.units     = 'Hours'
        var[:]        = steps
        # Create Variables
        var           = File.createVariable('Latitude','f',('time','lon','steps'))
        var.long_name = 'Latitude of particle position'
        var.units     = 'Degrees North'
        var[:]        = lats
        var           = File.createVariable('Longitude','f',('time','lon','steps'))
        var.long_name = 'Longitude of particle position'
        var.units     = 'Degress East'
        var[:]        = lons
        var           = File.createVariable('Pressure','f',('time','lon','steps'))
        var.long_name = 'Pressure of particle position'
        var.units     = 'hPa'
        var[:]        = Ps
        var           = File.createVariable('pv','f',('time','lon','steps'))
        var.long_name = 'Potential vorticity'
        var.units     = 'PVU'
        var[:]        = PVs
        # Close file
        File.close()
