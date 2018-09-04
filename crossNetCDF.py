from sys import *
from netCDF4 import *
import numpy as np

def crossNetCDF(time,lon,latstart,field,loncross,flux,step,FileName):

	ntime,nlon = len(time),len(lon)
        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'Time at initialisation'
        var.units     = 'hours since 1900-01-01 00:00:0.0'
        var[:]        = time
        # Create axes
        File.createDimension('lon',nlon)
        var           = File.createVariable('lon','f',('lon',))
        var.long_name = 'Longitude initialisation'
        var.units     = 'degrees_east'
        var[:]        = lon 
        # Create Variables
        var           = File.createVariable('Step','f',('time','lon'))
        var.long_name = 'Timesteps since initialisation'
        var.units     = '0-1'
        var[:]        = step
        var           = File.createVariable('Flux','f',('time','lon'))
        var.long_name = 'Meridional moisture flux at 70N at crossing'
        var.units     = 'Tg day**-1 deg**-1'
        var[:]        = flux
        var           = File.createVariable('LatStart','f',('time','lon'))
        var.long_name = 'Initialisation latitude'
        var.units     = 'Degrees North'
        var[:]        = latstart
        var           = File.createVariable('PW','f',('time','lon'))
        var.long_name = 'Initialisation precipitable water'
        var.units     = 'kg m**-2'
        var[:]        = field
        var           = File.createVariable('LonCross','f',('time','lon'))
        var.long_name = 'Crossing longitude at 70N'
        var.units     = 'Degrees East'
        var[:]        = loncross
        # Close file
        File.close()
