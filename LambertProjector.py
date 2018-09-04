import os,sys
from mpl_toolkits import basemap
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from progress import ProgressMeter
    
class LambertProjector:
    def __init__(self,
                 boundinglat = 60.,
                 resolution = 80., # nominal resolution in km
                 FileName   = None,
		 Projection = 'nplaea'	):

        # set up map
        self.m = Basemap(projection=Projection,
                         boundinglat=boundinglat, lon_0=0, resolution='i')
#       self.m = Basemap(projection='nplaea',
#                        boundinglat=boundinglat, lon_0=0, resolution='c',round=False)
        self.nx = int( (90.-boundinglat)*2.*self.m.rmajor*np.pi/180./resolution/1000.)
        # make sure number is even, to avoid point on pole which has undefined longitude
        self.nx = self.nx/2*2 
        self.ny = self.nx
        # dummy lon/lat just to get grid on map projectio
        lon = np.arange(-180,180)
        lat = np.arange(-90,91)
        field = np.ones((len(lat),len(lon)),dtype='float')
        field,self.x,self.y = \
          self.m.transform_scalar(field,lon,lat,self.nx,self.ny,returnxy=True)
        self.lon,self.lat = self.m(self.x,self.y,inverse=True)                        
	self.res = (self.x[0,1] - self.x[0,0])/1000

    def __call__(self,field,lon,lat):
        """
        input field on regular lat-lon grid
        output field on regular projection grid
        """
        if len(field.shape) == 2:
            field,lon = basemap.addcyclic(field,lon)
            field,lon = basemap.shiftgrid(180,field,lon,start=False)
            self.field = self.m.transform_scalar(field, lon, lat, self.nx, self.ny)

        elif len(field.shape) == 3:
            n = field.shape[0]
            self.field = np.zeros((n,self.ny,self.nx),dtype='f')
            for l in range(n):
                field1,lon1 = basemap.addcyclic(field[l],lon)
                field1,lon1 = basemap.shiftgrid(180,field1,lon1,start=False)
                self.field[l] = self.m.transform_scalar(field1, lon1, lat, self.nx, self.ny)

        elif len(field.shape) == 4:
            n0 = field.shape[0]
            n1 = field.shape[1]
            if hasattr(field,'mask'):
                self.field = np.ma.zeros((n0,n1,self.ny,self.nx),dtype=float)
            else:
                self.field = np.zeros((n0,n1,self.ny,self.nx),dtype=float)
            print 'LambertProjector: Projecting 4D field'
            m = ProgressMeter(total=n0*n1)
            for l0 in range(n0):
                for l1 in range(n1):
                    field1,lon1 = basemap.addcyclic(field[l0,l1],lon)
                    field1,lon1 = basemap.shiftgrid(180,field1,lon1,start=False)
                    self.field[l0,l1] = self.m.transform_scalar(field1, lon1, lat, self.nx, self.ny)
                    m.update(1)
                    
        return self.field

    def plot(self,contour_levels=30,title=''):
        self.m.contourf(self.x,self.y,self.field,contour_levels)
        #self.m.imshow(self.field,interpolation='none')
        self.m.drawcoastlines(color='k',linewidth=.5)
        self.m.drawparallels(np.arange(60.,80.,10.))
        self.m.drawmeridians(np.arange(-180.,181.,30.))
        self.m.colorbar()
        plt.title(title)

    def createOutFile(self,
                      FileName,
                      FieldName,
                      FieldLong_name = '-',
                      FieldUnits = '-',
                      TimeUnits = '-'):
        from netCDF4 import Dataset as NetCDFFile
        File = NetCDFFile(FileName,'w',format='NETCDF3_CLASSIC')
        # Create dims
        # time is record dimension
        File.createDimension('time',None)
        var = File.createVariable('time','d',('time',))
        var.long_name = 'time'
        var.units = TimeUnits
        # x
        File.createDimension('x',p.nx)
        var = File.createVariable('x','d',('x',))
        var.long_name = 'distance east'
        var.units = 'km'
        var[:] = (p.x[0]/1000.).astype('f')
        # y
        File.createDimension('y',p.nx)
        var = File.createVariable('y','f',('y',))
        var.long_name = 'distance north'
        var.units = 'km'
        var[:] = (p.y[:,0]/1000.).astype('f')
        # create variables
        var = File.createVariable('lon','f',('y','x'))
        var.long_name = 'longitude'
        var.units = 'degrees_east'
        var[:,:] = p.lon[:,:]
        var = File.createVariable('lat','f',('y','x'))
        var.long_name = 'latitude'
        var.units = 'degrees_north'
        var[:,:] = p.lat[:,:]
        var = File.createVariable(FieldName,'f',('time','y','x'))
        var.long_name = FieldLong_name
        var.units = FieldUnits
        var.missing_value=1.e20
        return File
        
if __name__ == '__main__':
    p = LambertProjector(resolution=80.)
    from ReanalysisMonthlyMeanDataServer import DataServer as MonthlyDataServer
    ci = MonthlyDataServer(Source='ERAInt',
                           LevType='surface_analysis',
                           Field='ci').getSeasonalClimatology('DJF')
    climo = MonthlyDataServer(Source='ERAInt',
                              LevType='surface_analysis',
                              Field='T2').getSeasonalClimatology('DJF')
    plt.ion()
    from ReanalysisDataServer import DataServer
    d = DataServer(Source='ERAInt',LevType='surface_analysis',Field='T2')
    h = d.getHours(1996,1,1,0)
    for i in range(1000):
        print i
        field = np.ma.array(d.snapshot(*d.getDate(h)) - climo)
        field.mask = ci.mask
        j = np.argmin(abs(d.lat-60.))
        field = p(field[j:,:],d.lon,d.lat[j:])
        p.plot(contour_levels=np.arange(-20,22,2),title=('%s '*4) % d.getDate(h))
        h += 6
        plt.show()
    
    ## File = createOutFile('jim.nc',p,d)
    ## File.variables['ci'][0] = field.filled().astype('float')
    ## File.variables['time'][0] = 1
    ## File.close()
