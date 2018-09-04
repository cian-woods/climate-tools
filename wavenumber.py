from ReanalysisDataServer import DataServer as reDataServer
from cmipDataServer import DataServer as cmipDataServer
from scipy import interpolate
from scipy import ndimage as nd
from netCDF4 import Dataset
from LambertProjector import *
from UnPickle import *

import os,sys,glob
import numpy as np
import matplotlib.pyplot as pl

def interpolateND(field,xold,xnew,axis):
        # field in (time)x(lev)x(lon)
        f     = interpolate.interp1d(xold,field,axis=axis,bounds_error=False,fill_value='extrapolate')
        field = f(xnew)
        return field

def fill(data,invalid=None):
         if invalid is None: invalid = np.isnan(data)
         ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
         return data[tuple(ind)]

def waveN(v,N):
        # Takes in v'F' (longitude as last dimension) and returns the wave(k) patterns
        # Divide by n becasue discrete transform and extract wave1 (nlon/2 + 1 index after shift)
        size = len(v.shape)
        nlon = v.shape[-1]/2 + 1
        Fk   = np.fft.fft(v,axis=size-1)/v.shape[-1]
        Fk   = np.fft.fftshift(Fk,axes=(size-1,))[...,nlon:nlon + N]
        # Calculate amplitide and phase of wavenumbers
        amp = 2*np.abs(Fk)
        ang = np.arctan2(Fk.imag,Fk.real)
        return amp,ang

"""
def waveN(v,N):
        # Takes in v'F' (longitude as last dimension) and returns the wave(k) patterns
        # Divide by n becasue discrete transform and extract wave1 (nlon/2 + 1 index after shift)
	size = len(v.shape)
	nlon = v.shape[-1]/2 + 1
        Fk   = np.fft.fft(v,axis=size-1)/v.shape[-1]
        Fk   = np.fft.fftshift(Fk,axes=(size-1,))[:,:,:,nlon:nlon + N]
        # Calculate amplitide and phase of wavenumbers
        amp = 2*np.abs(Fk)
        ang = np.arctan2(Fk.imag,Fk.real)
        return amp,ang
"""

def makeNetCDF1(field,times,wavenumbers,lats,time_units,FileName,case):
        print 'Creating %s ...' % (FileName)

        if case == 'vq': orr,units = 'moisture','kg s**-1 m**-1'
        if case == 'vT': orr,units = 'temperature','K kg s**-1 m**-1'
	if case == 'vv': orr,units = 'kinetic','kg s**-2'
	if case == 'vu': orr,units = 'kinetic','kg s**-2'
        ntime,nwavenumbers,nlat = len(times),len(wavenumbers),len(lats)

        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'Time'
        var.units     = time_units
        var[:]        = times
        # Horizontal axes
        File.createDimension('lat',nlat)
        var           = File.createVariable('lat','f',('lat',))
        var.long_name = 'Latitude'
        var.units     = 'degrees_north'
        var[:]        = lats.astype('f')
        File.createDimension('wavenumber',nwavenumbers)
        var           = File.createVariable('wavenumber','f',('wavenumber',))
        var.long_name = 'Zonal wavenumber'
        var.units     = '%s-%s' % (wavenumbers[0],wavenumbers[-1])
        var[:]        = wavenumbers.astype('f')
        # Create Variables
        var           = File.createVariable(case,'f',('time','lat','wavenumber'))
        var.long_name = 'Vertically integrated zonal-mean meridional %s flux' % (orr)
        var.units     = units
        var[:]        = field
        # Close file
        File.close()

def makeNetCDF2(field1,field2,times,wavenumbers,lats,levs,wvn_units,time_units,FileName,case):
	print 'Creating %s ...' % (FileName)

	if case == 'v':  fieldname,fieldunits = 'meridional wind','m s**-1'
	if case == 'q':  fieldname,fieldunits = 'specific humidity','kg/kg'
	if case == 'vq': fieldname,fieldunits = 'Meridional moisture flux','kg s**-1 m**-1'
	if case ==  'T': fieldname,fieldunits = 'temperature','K'
	ntime,nwavenumbers,nlat,nlev = len(times),len(wavenumbers),len(lats),len(levs)

	# Create file
	File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
	# Define some global attribs
	File.Conventions='COARDS'
	# Time is record dimension
	File.createDimension('time',ntime)
	var           = File.createVariable('time','d',('time',))
	var.long_name = 'Time'
	var.units     = time_units
	var[:]        = times
	# Horizontal axes
        File.createDimension('lev',nlev)
        var           = File.createVariable('lev','f',('lev',))
        var.long_name = 'Pressure level'
        var.units     = 'hPa'
        var[:]        = levs.astype('f')
	File.createDimension('lat',nlat)
	var           = File.createVariable('lat','f',('lat',))
	var.long_name = 'Latitude'
	var.units     = 'degrees_north'
	var[:]        = lats.astype('f')
	File.createDimension('wavenumber',nwavenumbers)
	var           = File.createVariable('wavenumber','f',('wavenumber',))
	var.long_name = 'Zonal wavenumber'
	var.units     = wvn_units
	var[:]        = wavenumbers.astype('f')
	# Create Variables
	var           = File.createVariable('Amp','f',('time','lev','lat','wavenumber'))
	var.long_name = 'Amplitude of %s' % (fieldname)
	var.units     = fieldunits
	var[:]        = field1
        var           = File.createVariable('Phi','f',('time','lev','lat','wavenumber'))
        var.long_name = 'Phase of %s' % (fieldname)
        var.units     = 'Radians'
        var[:]        = field2
	# Close file
	File.close()

def makeNetCDF3(field,times,levs,lats,lons,wavenumbers,time_units,FileName,case):
        print 'Creating %s ...' % (FileName)

        if case == 'v'   : orr,units = 'Meridional velocity','m s**-1'
        if case == 'q'   : orr,units = 'Specific humidity','kg kg**-1'
	if case == 'vq'  : orr,units = 'Meridional specific humidity flux','m s**-1 kg kg**-1'
        ntime,nlev,nwavenumbers,nlat,nlon = len(times),len(levs),len(wavenumbers),len(lats),len(lons)

        # Create file
        File = Dataset(FileName,'w',format='NETCDF3_CLASSIC')
        # Define some global attribs
        File.Conventions='COARDS'
        # Time is record dimension
        File.createDimension('time',ntime)
        var           = File.createVariable('time','d',('time',))
        var.long_name = 'Time'
        var.units     = time_units
        var[:]        = times
        # Spatial axes
        File.createDimension('lev',nlev)
        var           = File.createVariable('lev','f',('lev',))
        var.long_name = 'Pressure level'
        var.units     = 'hPa'
        var[:]        = levs.astype('f')
        File.createDimension('lat',nlat)
        var           = File.createVariable('lat','f',('lat',))
        var.long_name = 'Latitude'
        var.units     = 'degrees_north'
        var[:]        = lats.astype('f')
        File.createDimension('lon',nlon)
        var           = File.createVariable('lon','f',('lon',))
        var.long_name = 'Longitude'
        var.units     = 'degrees_east'
        var[:]        = lons.astype('f')
        File.createDimension('wavenumber',nwavenumbers)
        var           = File.createVariable('wavenumber','f',('wavenumber',))
        var.long_name = 'Zonal wavenumber'
        var.units     = '1-15'
        var[:]        = wavenumbers.astype('f')
        # Create Variables
        var           = File.createVariable(case,'f',('time','lev','lat','lon','wavenumber'))
        var.long_name = '%s' % (orr)
        var.units     = units
        var[:]        = field
        # Close file
        File.close()

def extendField(field):
        fieldend = field[...,-1][...,np.newaxis]
        field    = np.append(field,fieldend,axis=-1)
        field    = np.append(fieldend,field,axis=-1)
        return field

"""
def extendField(field):
	fieldend = field[:,:,:,-1][:,:,:,np.newaxis]
	field    = np.append(field,fieldend,axis=3)
	field    = np.append(fieldend,field,axis=3)
	return field
"""

def test(Source,Year,Month):
	Nwv    = np.arange(1,15+1,1) - 1
	proj   = LambertProjector(boundinglat=20,resolution=200.)
	if Source != 'ERAInt':
		fvq      = Dataset('../rcp85/%s/mon/surface/vq/vq.mon.mean.nc' % (Source),'r')
		fvkqk    = Dataset('../rcp85/%s/mon/plev/vkqk_resynth/vkqk_resynth_%s.nc' % (Source,Year),'r')
		fvqk     = Dataset('../rcp85/%s/mon/plev/vqk_resynth/vqk_resynth_%s.nc' % (Source,Year),'r')
		fvqk_pos = Dataset('../rcp85/%s/mon/plev/vqk+_resynth/vqk+_resynth_%s.nc' % (Source,Year),'r')
		fvqk_neg = Dataset('../rcp85/%s/mon/plev/vqk-_resynth/vqk-_resynth_%s.nc' % (Source,Year),'r')
		dv       = cmipDataServer(Field='va',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		dq       = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
	if Source == 'ERAInt':
		fvq      = Dataset('../vq.mon.mean.nc','r')
		fvkqk    = Dataset('../vkqk_resynth/vkqk_resynth_%s.nc' % (Year),'r')
		fvqk     = Dataset('../vqk_resynth/vqk_resynth_%s.nc' % (Year),'r')
		fvqk_pos = Dataset('../vqk+_resynth/vqk+_resynth_%s.nc' % (Year),'r')
		fvqk_neg = Dataset('../vqk-_resynth/vqk-_resynth_%s.nc' % (Year),'r')
		dv       = reDataServer( Field='V',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		dq       = reDataServer( Field='q',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
	long_name = fvkqk.variables['vq'].long_name
	units     = fvkqk.variables['vq'].units
	levs1     = fvkqk.variables['lev'][:]
	lats1     = fvkqk.variables['lat'][:]
	lons1     = fvkqk.variables['lon'][:]
	wvns1     = fvkqk.variables['wavenumber'][Nwv]

	# Mass of girdbox on each level
	dP_ = []
	dP  = np.diff(levs1)/2.
	dP_.append(dP[0])
	for i in range(len(dP)-1):
	        dP_.append(dP[i]+dP[i+1])
	dP_.append(dP[-1])
	dP_ = np.array(dP_)
	dM  = 100.*dP_/9.80665
        # Mass of girdbox on each level
        dP_ = []
        dP  = np.diff(dv.lev)/2.
        dP_.append(dP[0])
        for i in range(len(dP)-1):
                dP_.append(dP[i]+dP[i+1])
        dP_.append(dP[-1])
        dP_  = np.array(dP_)
        dMe  = 100.*dP_/9.80665

	# Montly mean zonal-mean field
	years    = range(1979,2016+1,1)
	vqzm     = fvq.variables['vq'][:,:,:].reshape(-1,12,len(lats1),16)[years.index(Year),Month-1,:,:].sum(axis=-1)
	# Resynthesised fields
	vq1      = (   fvkqk.variables['vq'][Month-1,:,:,:,Nwv].sum(axis=-1)*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
	vq1k     = (    fvqk.variables['vq'][Month-1,:,:,:,Nwv].sum(axis=-1)*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
	vq1k_pos = (fvqk_pos.variables['vq'][Month-1,:,:,:,Nwv].sum(axis=-1)*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
	vq1k_neg = (fvqk_neg.variables['vq'][Month-1,:,:,:,Nwv].sum(axis=-1)*dM[:,np.newaxis,np.newaxis]).sum(axis=0)
	# Full fields
	v2  = dv.getDataSnaps(Year=Year,Month=Month)[:,:,:,:]
	q2  = dq.getDataSnaps(Year=Year,Month=Month)[:,:,:,:]
	#v2  = v2 - v2.mean(axis=-1)[:,:,np.newaxis]
	#q2  = q2 - q2.mean(axis=-1)[:,:,np.newaxis]
	if (type(v2)==np.ma.core.MaskedArray) and (v2.mask.shape == v2.data.shape):
		v2 = fill(v2.data,invalid=v2.mask)
        if (type(q2)==np.ma.core.MaskedArray) and (q2.mask.shape == q2.data.shape):
                q2 = fill(q2.data,invalid=q2.mask)
	vq2  = ((v2*q2).mean(axis=0)*dMe[:,np.newaxis,np.newaxis]).sum(axis=0)
	#vq2 = vq2 - vq2.mean(axis=-1)[:,np.newaxis]
	#v2  = v2.mean(axis=0)
	#q2  = q2.mean(axis=0)

	# Project to Lambert
	vq1       = proj(vq1,lons1,lats1)
	vq1k_full = proj(vq1k+vqzm[:,np.newaxis],lons1,lats1)
	vq1k      = proj(vq1k,lons1,lats1)
	vq1k_pos  = proj(vq1k_pos,lons1,lats1)
	vq1k_neg  = proj(vq1k_neg,lons1,lats1)
	vq2       = proj(vq2,dv.lon,dv.lat)

	# Plot
	#cseq = np.arange(-160,160+20,20)
	cseq = 14
	pl.figure(1)
	cf = pl.contourf(proj.x,proj.y,vq1,cseq,cmap=pl.cm.RdBu_r,extend='both')
	cbar = pl.colorbar(cf)
	proj.m.drawcoastlines()
	proj.m.drawparallels([70,80],latmax=90)
	pl.title('Resynthesised v(k)q(k): (%s, %s)' % (Year,Month))

	pl.figure(2)
        cf = pl.contourf(proj.x,proj.y,vq2,cseq,cmap=pl.cm.RdBu_r,extend='both')
        cbar = pl.colorbar(cf) 
        proj.m.drawcoastlines()
        proj.m.drawparallels([70,80],latmax=90)
	pl.title('Full: (%s, %s)' % (Year,Month))

        pl.figure(3)
        cf = pl.contourf(proj.x,proj.y,vq1k,cseq,cmap=pl.cm.RdBu_r,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawcoastlines()
        proj.m.drawparallels([70,80],latmax=90)
        pl.title('Resynthesised zonal anomaly vq(k): (%s, %s)' % (Year,Month))

        pl.figure(4)
        cf = pl.contourf(proj.x,proj.y,vq1k_full,cseq,cmap=pl.cm.RdBu_r,extend='both')
        cbar = pl.colorbar(cf)
        proj.m.drawcoastlines()
        proj.m.drawparallels([70,80],latmax=90)
        pl.title('Resynthesised full vq(k): (%s, %s)' % (Year,Month))

	pl.show()

def resynthMonthlyFiles_vq(Source):
        # DataServers and axes from wavenumber file
        proj = LambertProjector(boundinglat=30,resolution=200.)
        if Source != 'ERAInt': dv  = cmipDataServer(Field='va',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
        if Source == 'ERAInt': dv  = reDataServer(  Field='V',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
        f    = Dataset('../%s_decomp/%s_decomp_%s.nc' % ('V','V',1979),'r')
        wvns = f.variables['wavenumber'][0:15]
        lats = f.variables['lat'][:]
        levs = f.variables['lev'][:]
        lons = np.arange(0,360,1)
	f.close()

        years = range(1980,1980+1,1)
        for year in years:
                if Source == 'ERAInt': Dir,fname = '../vkqk_resynth','vkqk_resynth_%s.nc' % (year)
                if Source != 'ERAInt': Dir,fname = '../rcp85/%s/mon/plev/vkqk_resynth' % (Source),'vkqk_resynth_%s.nc' % (year)
                if not os.path.isfile('%s/%s' % (Dir,fname)):
                        SMM   = np.zeros((15,360,0,11,86))
                        times = []
                        print year
                        # Decomposition
                        if Source != 'ERAInt':
				fv = Dataset('../rcp85/%s/day/plev/v_decomp/v_decomp_%s.nc' % (Source,year),'r')
				fq = Dataset('../rcp85/%s/day/plev/q_decomp/q_decomp_%s.nc' % (Source,year),'r')
                        if Source == 'ERAInt':
				fv = Dataset('../v_decomp/v_decomp_%s.nc' % (year),'r')
				fq = Dataset('../q_decomp/q_decomp_%s.nc' % (year),'r')
                        Ampv,Phiv = fv.variables['Amp'][:,:,:,0:15],fv.variables['Phi'][:,:,:,0:15]
			Ampq,Phiq = fq.variables['Amp'][:,:,:,0:15],fq.variables['Phi'][:,:,:,0:15]
			fv.close()
			fq.close()
                        Sv,Sq = [],[]
                        # Resynth
                        for i in wvns:
                                sf = -1
                                if i%2 == 0: sf = 1
                                Sv.append([sf*Ampv[:,:,:,i-1]*np.cos(j + Phiv[:,:,:,i-1]) for j in np.linspace(-i*np.pi,i*np.pi,len(lons)+1)[0:-1]])
				Sq.append([sf*Ampq[:,:,:,i-1]*np.cos(j + Phiq[:,:,:,i-1]) for j in np.linspace(-i*np.pi,i*np.pi,len(lons)+1)[0:-1]])
                        Sv,Sq = np.array(Sv),np.array(Sq)
			S     = Sv*Sq
        
                        datelist = dv.getDateList(Year=year,Season='Annual')
                        if Source == 'ERAInt':
                                datelist = datelist[::4]
                        months   = np.array([date[1] for date in datelist])     
                        for month in range(1,12+1,1):
                                xs  = np.where(months==month)[0]
                                Smm = S[:,:,xs,:,:].mean(axis=2)[:,:,np.newaxis,:,:]
                                SMM = np.append(SMM,Smm,axis=2)
                                times.append([datelist[ix] for ix in xs][14])

                        if Source == 'ERAInt': times = [dv.getHours(*time) for time in times]
                        if Source != 'ERAInt': times = [ dv.getDays(*time) for time in times]
                        SMM = np.rollaxis(SMM,2,0)
                        SMM = np.rollaxis(SMM,3,1)
                        SMM = np.rollaxis(SMM,4,2)
                        SMM = np.rollaxis(SMM,4,3)
                        # Make netcdf file   
                        os.system('mkdir -p %s' % (Dir))
                        FileName = '%s/%s' % (Dir,fname)
                        makeNetCDF3(SMM,times,levs,lats,lons,wvns,dv.time_units,FileName,'vq')
                else:
                        print '\nFile %s/%s already exists. Did nothing...\n' % (Dir,fname)

def resynthMonthlyFiles(Source,case):
        # Use daily amplitude and phase files to resynthesise waves into montly mean files
        # DataServers and axes from wavenumber file
        if case not in ['v','q','vqk']:
                print '\ncase is not suitable. Exiting...\n'
                sys.exit()
        if case == 'v'  : head = 'v'
        if case == 'q'  : head = 'q'
        if case == 'vqk': head = 'vq'
        proj = LambertProjector(boundinglat=30,resolution=200.)
        if Source != 'ERAInt': dv  = cmipDataServer(Field='va',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
        if Source == 'ERAInt': dv  = reDataServer(  Field='V',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
        f    = Dataset('../rcp85/BNU-ESM/day/plev/%s_decomp/%s_decomp_%s.nc' % (head,head,2002),'r')
        wvns = f.variables['wavenumber'][0:15]
        lats = f.variables['lat'][:]
        levs = f.variables['lev'][:]
        lons = np.arange(0,360,1)
        nlev,nlat,nlon,nwvn = len(levs),len(lats),len(lons),len(wvns)

        years = range(1980,2005+1,1)
        for year in years:
                if Source == 'ERAInt': Dir,fname = '../%s_resynth' % (case),'%s_resynth_%s.nc' % (case,year)
                if Source != 'ERAInt': Dir,fname = '../rcp85/%s/mon/plev/%s_resynth' % (Source,case),'%s_resynth_%s.nc' % (case,year)
                if not os.path.isfile('%s/%s' % (Dir,fname)):
                        SMM   = np.zeros((nwvn,nlon,0,nlev,nlat))
                        times = []
                        print year
                        # Decomposition
                        if Source != 'ERAInt': f = Dataset('../rcp85/%s/day/plev/%s_decomp/%s_decomp_%s.nc' % (Source,head,head,year),'r')
                        if Source == 'ERAInt': f = Dataset('../%s_decomp/%s_decomp_%s.nc' % (head,head,year),'r')
                        Amp,Phi = f.variables['Amp'][:,:,:,0:15],f.variables['Phi'][:,:,:,0:15]
                        S = []
                        # Resynth
                        for i in wvns:
                                sf = -1
                                if i%2 == 0: sf = 1
                                S.append([sf*Amp[:,:,:,i-1]*np.cos(j + Phi[:,:,:,i-1]) for j in np.linspace(-i*np.pi,i*np.pi,len(lons)+1)[0:-1]])
                        S = np.array(S)

                        datelist = dv.getDateList(Year=year,Season='Annual')
                        if Source == 'ERAInt':
                                datelist = datelist[::4]
                        months = np.array([date[1] for date in datelist])
                        for month in range(1,12+1,1):
                                xs  = np.where(months==month)[0]
                                Smm = S[:,:,xs,:,:].mean(axis=2)[:,:,np.newaxis,:,:]
                                SMM = np.append(SMM,Smm,axis=2)
                                times.append([datelist[ix] for ix in xs][14])

                        if Source == 'ERAInt': times = [dv.getHours(*time) for time in times]
                        if Source != 'ERAInt': times = [ dv.getDays(*time) for time in times]
                        SMM = np.rollaxis(SMM,2,0)
                        SMM = np.rollaxis(SMM,3,1)
                        SMM = np.rollaxis(SMM,4,2)
                        SMM = np.rollaxis(SMM,4,3)
                        # Make netcdf file   
                        os.system('mkdir -p %s' % (Dir))
                        FileName = '%s/%s' % (Dir,fname)
                        makeNetCDF3(SMM,times,levs,lats,lons,wvns,dv.time_units,FileName,head)
                else:
                        print '\nFile %s/%s already exists. Did nothing...\n' % (Dir,fname)

def resynthMonthlyFiles_NS(Source,case,phase):
	# Use daily amplitude and phase files to resynthesise waves into montly mean files of northward and southward v or vq
	# DataServers and axes from wavenumber file
	if case not in ['v','vqk']:
		print '\ncase is not suitable. Exiting...\n'
		sys.exit()
	if case == 'v'  : head = 'v'
	if case == 'vqk': head = 'vq'
	proj = LambertProjector(boundinglat=30,resolution=200.)
	if Source != 'ERAInt': dv  = cmipDataServer(Field='va',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
	if Source == 'ERAInt': dv  = reDataServer(  Field='V',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(250,500))
	f    = Dataset('../rcp85/BNU-ESM/day/plev/%s_decomp/%s_decomp_%s.nc' % (head,head,2002),'r')
	wvns = f.variables['wavenumber'][0:15]
	lats = f.variables['lat'][:]
	levs = f.variables['lev'][:]
	lons = np.arange(0,360,1)
	nlev,nlat,nlon,nwvn = len(levs),len(lats),len(lons),len(wvns)

	years = range(1980,2005+1,1)
	for year in years:
                if Source == 'ERAInt': Dir,fname = '../%s%s_resynth' % (case,phase),'%s%s_resynth_%s.nc' % (case,phase,year)
                if Source != 'ERAInt': Dir,fname = '../rcp85/%s/mon/plev/%s%s_resynth' % (Source,case,phase),'%s%s_resynth_%s.nc' % (case,phase,year)
		if not os.path.isfile('%s/%s' % (Dir,fname)):
			SMM   = np.zeros((nwvn,nlon,0,nlev,nlat))
			times = []
			print year
			# Decomposition
			if Source != 'ERAInt': f = Dataset('../rcp85/%s/day/plev/%s_decomp/%s_decomp_%s.nc' % (Source,head,head,year),'r')
			if Source == 'ERAInt': f = Dataset('../%s_decomp/%s_decomp_%s.nc' % (head,head,year),'r')
			Amp,Phi = f.variables['Amp'][:,:,:,0:15],f.variables['Phi'][:,:,:,0:15]
			S = []
			# Resynth
			for i in wvns:
				sf = -1
				if i%2 == 0: sf = 1
				S.append([sf*Amp[:,:,:,i-1]*np.cos(j + Phi[:,:,:,i-1]) for j in np.linspace(-i*np.pi,i*np.pi,len(lons)+1)[0:-1]])
			S = np.array(S)	
	
			datelist = dv.getDateList(Year=year,Season='Annual')
			if Source == 'ERAInt':
				datelist = datelist[::4]
			months = np.array([date[1] for date in datelist])	
			for month in range(1,12+1,1):
				xs    = np.where(months==month)[0]
				ndays = len(xs)
				Smm   = S[:,:,xs,:,:]
				if phase == '+': Smm[np.where(Smm<0)] = 0
				if phase == '-': Smm[np.where(Smm>0)] = 0
				Smm   = Smm.sum(axis=2)[:,:,np.newaxis,:,:]/ndays
				SMM   = np.append(SMM,Smm,axis=2)
				times.append([datelist[ix] for ix in xs][14])
	
		        if Source == 'ERAInt': times = [dv.getHours(*time) for time in times]
		        if Source != 'ERAInt': times = [ dv.getDays(*time) for time in times]
		        SMM = np.rollaxis(SMM,2,0)
		        SMM = np.rollaxis(SMM,3,1)
		        SMM = np.rollaxis(SMM,4,2)
		        SMM = np.rollaxis(SMM,4,3)
		        # Make netcdf file   
		        os.system('mkdir -p %s' % (Dir))
		        FileName = '%s/%s' % (Dir,fname)
		        makeNetCDF3(SMM,times,levs,lats,lons,wvns,dv.time_units,FileName,head)
		else:
			print '\nFile %s/%s already exists. Did nothing...\n' % (Dir,fname)

def Amplitude_Phase_DailyFiles(Source,case):
        if case not in ['v','q','T']:
                print '\ncase %s is not suitable. Exiting...\n' % (case)
                sys.exit()
	# Field
	if case == 'v': Field = ['V','va']
	if case == 'q': Field = ['q','hus']
	if case == 'T': Field = ['T','ta']
	# Interpolation axes
	lons    = np.arange(0,360,1)
	lats    = np.arange(0,85+1,1)
	levs    = np.arange(0,1000+100,100)
	# DataServers
	if Source == 'ERAInt':
	        dv = reDataServer(Field=Field[0],LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
	else:
	        dv = cmipDataServer(Field=Field[1],LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
	# Extend lon axis
	lonres = np.abs(dv.lon[1]-dv.lon[0])
	dlon   = np.append(dv.lon,dv.lon[-1]+lonres)
	dlon   = np.append(dv.lon[0]-lonres,dlon)
	# Wavenumber axis and year and month ranges
	N           = 179
	wvn_units   = '1-%s' % (N)
	wavenumbers = np.arange(1,N+1,1)
	years       = range(1980,2005+1,1)
	# Make data file directory
	if Source != 'ERAInt': Dir = '../rcp85/%s/day/plev/%s_decomp' % (Source,case)
	if Source == 'ERAInt': Dir = '../%s_decomp' % (Field[0])
	print Dir
	os.system('mkdir -p %s' % (Dir))
	for year in years:
	        print year
		FileName = '%s/%s_decomp_%s.nc' % (Dir,case,year)
		time     = [dv.getDays(*date) for date in dv.getDateList(Year=year,Season='Annual')]
		# Meridional wind
		F = dv.getDataSnaps(Year=year,Season='Annual')
		if (type(F)==np.ma.core.MaskedArray) and (F.mask.shape == F.data.shape):
		        F = fill(F.data,invalid=F.mask)
		if Source == 'ERAInt':
			print F.shape,len(time)
			F,time = F.reshape(-1,4,dv.nlev,dv.nlat,dv.nlon).mean(axis=1),time[::4]
			print F.shape,len(time)
		F   = extendField(F)	
		F   = interpolateND(F,dlon,lons,axis=3)
		F   = interpolateND(F,dv.lat,lats,axis=2)
		F   = interpolateND(F,dv.lev,levs,axis=1)
		Fzm = F.mean(axis=-1)[:,:,:,np.newaxis]
		Fza = F - Fzm
		Famp,Fang = waveN(Fza,N=N)
		print Famp.shape
		makeNetCDF2(Famp,Fang,time,wavenumbers,lats,levs,wvn_units,dv.time_units,FileName,case)

def Amplitude_Phase_DailyFiles_vq(Source):
        # Field
        # Interpolation axes
        lons    = np.arange(0,360,1)
        lats    = np.arange(0,85+1,1)
        levs    = np.arange(0,1000+100,100)
        # DataServers
        if Source == 'ERAInt':
                dv = reDataServer(Field='V',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
		dq = reDataServer(Field='q',LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
        else:
                dv = cmipDataServer(Field='va',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		dq = cmipDataServer(Field='hus',LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
        # Extend lon axis
        lonres = np.abs(dv.lon[1]-dv.lon[0])
        dlon   = np.append(dv.lon,dv.lon[-1]+lonres)
        dlon   = np.append(dv.lon[0]-lonres,dlon)
        # Wavenumber axis and year and month ranges
        N           = 179
        wvn_units   = '1-%s' % (N)
        wavenumbers = np.arange(1,N+1,1)
        years       = range(1979,2016+1,1)
        # Make data file directory
        if Source != 'ERAInt': Dir,orr = '../rcp85/%s/day/plev/vq_decomp' % (Source), 'vq'
        if Source == 'ERAInt': Dir,orr = '../vq_decomp', 'vq'
        print Dir
        os.system('mkdir -p %s' % (Dir))
        for year in years:
                print year
                FileName = '%s/%s_decomp_%s.nc' % (Dir,orr,year)
		if not os.path.isfile(FileName):
                	time     = [dv.getDays(*date) for date in dv.getDateList(Year=year,Season='Annual')]
                	# Meridional wind and humidity
                	V = dv.getDataSnaps(Year=year,Season='Annual')
			Q = dq.getDataSnaps(Year=year,Season='Annual')	
                	if (type(V)==np.ma.core.MaskedArray) and (V.mask.shape == V.data.shape):
                	        V = fill(V.data,invalid=V.mask)
                        if (type(Q)==np.ma.core.MaskedArray) and (Q.mask.shape == Q.data.shape):
                                Q = fill(Q.data,invalid=Q.mask)
			F = V*Q
                	if Source == 'ERAInt':
                	        print F.shape,len(time)
                	        F,time = F.reshape(-1,4,dv.nlev,dv.nlat,dv.nlon).mean(axis=1),time[::4]
                	        print F.shape,len(time)
                	F   = extendField(F)    
                	F   = interpolateND(F,dlon,lons,axis=3)
                	F   = interpolateND(F,dv.lat,lats,axis=2)
                	F   = interpolateND(F,dv.lev,levs,axis=1)
                	Fzm = F.mean(axis=-1)[:,:,:,np.newaxis]
                	Fza = F - Fzm
                	Famp,Fang = waveN(Fza,N=N)
                	print Famp.shape
                	makeNetCDF2(Famp,Fang,time,wavenumbers,lats,levs,wvn_units,dv.time_units,FileName,'vq')
		else:
			print '\nFile %s already exists, did nothing...\n' % (FileName)

def waveNumberFile(Source,lons,lats,levs,dM,case):
	if Source != ('ERAInt') and (Source[0:3] != 'CAM'):
	        Dir = '/mnt/climstorage/cian/rcp85/%s/mon/surface/%s' % (Source,case)
	        os.system('mkdir -p %s' % (Dir))
	elif Source[0:3] == 'CAM':
		Dir = '/mnt/climstorage/cian/EXPS/Expt_%s/h0/%s' % (Source[3:],case)
		os.system('mkdir -p %s' % (Dir))
	else:
	        Dir = '/mnt/climstorage/cian'
	FileName = '%s/%s.mon.mean.nc' % (Dir,case)
	print FileName
	if not os.path.isfile(FileName):
		# DataServers
		if case == 'vT': scalef,Fieldf =    1005/1e06,['T',  'ta']
		if case == 'vq': scalef,Fieldf = 2260e03/1e06,['q', 'hus']
		if case == 'vv': scalef,Fieldf =       2/1e06,['V',  'va']
		if case == 'vu': scalef,Fieldf =       2/1e06,['U',  'ua']
		if Source == 'ERAInt':
			dv = reDataServer(Field='V'      ,LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
			df = reDataServer(Field=Fieldf[0],LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
		elif Source[0:3] == 'CAM':
			dv = reDataServer(Field='V',Source=Source)
			df = reDataServer(Field=Fieldf[0],Source=Source)
		else:
		        dv = cmipDataServer(Field= 'va'    ,LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		        df = cmipDataServer(Field=Fieldf[1],LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		# Extend lon axis
		lonres = np.abs(dv.lon[1]-dv.lon[0])
		dlon   = np.append(dv.lon,dv.lon[-1]+lonres)
		dlon   = np.append(dv.lon[0]-lonres,dlon)
	
		N             = 179
		#N            = 15
		wavenumbers   = np.arange(0,N+1,1)
		#years,months = range(1901,1919+1,1),range(1,12+1,1)
		#years,months = range(1979,2016+1,1),range(1,12+1,1)
		#years,months = range(1980,1980+1,1),range(1,1+1,1)
		years,months = range(1980,2100+1,1),range(1,12+1,1)
		#years,months  = [1979],range(1,12+1,1)
		VFWN          = np.zeros((0,len(lats),len(wavenumbers)))
		time          = []
		for year in years:
			for month in months:
				print year,month
				time.append(dv.getDays(year,month,1,0))
				# Meridional wind
				V = dv.getDataSnaps(Year=year,Month=month,dailymean=True)
				if (type(V)==np.ma.core.MaskedArray) and (V.mask.shape == V.data.shape):
				        V = fill(V.data,invalid=V.mask)
				V   = extendField(V)
				V   = interpolateND(V,dlon,lons,axis=3)
				V   = interpolateND(V,dv.lat,lats,axis=2)
				V   = interpolateND(V,dv.lev,levs,axis=1)	
				Vzm = V.mean(axis=3)[:,:,:,np.newaxis]
				Vza = V - Vzm
				Vamp,Vang = waveN(Vza,N=N)
				# Specific humidity/Temperature
				F = df.getDataSnaps(Year=year,Month=month,dailymean=True)
				if (type(F)==np.ma.core.MaskedArray) and (F.mask.shape == F.data.shape):
				        F = fill(F.data,invalid=F.mask)
				F   = extendField(F)
				F   = interpolateND(F,dlon,lons,axis=3)
				F   = interpolateND(F,dv.lat,lats,axis=2)
				F   = interpolateND(F,dv.lev,levs,axis=1)
				Fzm = F.mean(axis=3)[:,:,:,np.newaxis]
				Fza = F - Fzm
				Famp,Fang = waveN(Fza,N=N)
				# Meridional moisture flux decompsed into wave numbers
				VF0  = Vzm*Fzm
				VFwn = (Vamp*Famp/2.)*np.cos(Vang-Fang)
				VFwn = np.append(VF0,VFwn,axis=3)
				# Take time mean and vertical integral
				VFwn = VFwn.mean(axis=0)
				#VFwn = (VFwn*dM[:,np.newaxis,np.newaxis]).sum(axis=0)[np.newaxis,:,:]
				VFwn = VFwn.mean(axis=0)[np.newaxis,:,:]*1000e02/9.81
				VFWN = np.append(VFWN,VFwn,axis=0)
		time = np.array(time)
		# Make NetCDF file
		makeNetCDF1(VFWN,time,wavenumbers,lats,dv.time_units,FileName,case)
	else:
		print '\nFile %s already exists, did nothing...\n' % (FileName)

def waveNumberFile_Freq(Source,lons,lats,levs,dM,case,Freq='6hourly'):
	if Freq == '6hourly':   dailymean,step = False,6
	if Freq == 'dailymean': dailymean,step = True,24
	# Directory
	if Source == 'ERAInt':
        	Dir = '/mnt/climstorage/cian/%s_%s' % (case,Freq)
	else:
		Dir = '/mnt/climstorage/cian/rcp85/%s/mon/surface/%s_%s' % (Source,case,Freq)
	os.system('mkdir -p %s' % (Dir))
	# DataServers
	if case == 'vT': scalef,Fieldf =    1005/1e06,['T',  'ta']
	if case == 'vq': scalef,Fieldf = 2260e03/1e06,['q', 'hus']
	if case == 'vv': scalef,Fieldf =       2/1e06,['V',  'va']
	if case == 'vu': scalef,Fieldf =       2/1e06,['U',  'ua']
	if Source == 'ERAInt':
		dv = reDataServer(Field='V'      ,LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
		df = reDataServer(Field=Fieldf[0],LevType='plev',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000),HybToPress=False)
	else:
		dv = cmipDataServer(Field= 'va'    ,LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))
		df = cmipDataServer(Field=Fieldf[1],LevType='plev',Source=Source,ExpType='rcp85',DataFreq='day',LatRange=(-10,90),LonRange=(0,360),LevRange=(0,1000))	
	# Extend lon axis
	lonres = np.abs(dv.lon[1]-dv.lon[0])
	dlon   = np.append(dv.lon,dv.lon[-1]+lonres)
	dlon   = np.append(dv.lon[0]-lonres,dlon)

	N             = 15
	wavenumbers   = np.arange(0,N+1,1)
	years,months  = range(1979,2017+1,1),range(1,12+1,1)
	for year in years:
	        for month in months:
			FileName = '%s/vq_%s_%02d.nc' % (Dir,year,month)
			if not os.path.isfile(FileName):
				try:
		       	        	time = np.array(dv.getHourList(Year=year,Month=month,step=step))	
		       	        	# Meridional wind
		       	        	V = dv.getDataSnaps(Year=year,Month=month,dailymean=dailymean)
		                	if (type(V)==np.ma.core.MaskedArray) and (V.mask.shape == V.data.shape):
		                	        V = fill(V.data,invalid=V.mask)
		                	V   = extendField(V)
		                	V   = interpolateND(V,dlon,lons,axis=3)
		                	V   = interpolateND(V,dv.lat,lats,axis=2)
		                	V   = interpolateND(V,dv.lev,levs,axis=1)
		                	Vzm = V.mean(axis=3)[:,:,:,np.newaxis]
		                	Vza = V - Vzm
		                	Vamp,Vang = waveN(Vza,N=N)
		                	# Specific humidity/Temperature
		                	F = df.getDataSnaps(Year=year,Month=month,dailymean=dailymean)
		                	if (type(F)==np.ma.core.MaskedArray) and (F.mask.shape == F.data.shape):
		                	        F = fill(F.data,invalid=F.mask)
		                	F   = extendField(F)
		                	F   = interpolateND(F,dlon,lons,axis=3)
		                	F   = interpolateND(F,dv.lat,lats,axis=2)
		                	F   = interpolateND(F,dv.lev,levs,axis=1)
		                	Fzm = F.mean(axis=3)[:,:,:,np.newaxis]
		                	Fza = F - Fzm
		                	Famp,Fang = waveN(Fza,N=N)
		                	# Meridional moisture flux decompsed into wave numbers
		                	VF0  = Vzm*Fzm
		                	VFwn = (Vamp*Famp/2.)*np.cos(Vang-Fang)
		                	VFwn = np.append(VF0,VFwn,axis=3)
		                	# Take vertical integral 
		                	VFwn = VFwn.mean(axis=1)*1000e02/9.81
					# Make NetCDF file
					makeNetCDF1(VFwn,time,wavenumbers,lats,dv.time_units,FileName,case)
				except:
					print 'File %s failed to complete ...' % (FileName)
			else:
				print 'File %s already exists. Did nothing ...' % (FileName)

def waveNumberFile_stationary(Source,YearRange,Season,lons,lats,levs,case):
        if Source != ('ERAInt') and (Source[0:3] != 'CAM'):
                Dir = '/mnt/climstorage/cian/rcp85/%s/mon/surface/%s' % (Source,case)
                os.system('mkdir -p %s' % (Dir))
        elif Source[0:3] == 'CAM':
                Dir = '/mnt/climstorage/cian/EXPS/Expt_%s/h0/%s' % (Source[3:],case)
                os.system('mkdir -p %s' % (Dir))
        else:
                Dir = '/mnt/climstorage/cian/%s' % (case)
        FileName = '%s/%s.stat.%s-%s.%s.mon.mean.nc' % (Dir,case,YearRange[0],YearRange[1],Season)
        if not os.path.isfile(FileName):
                # Attributes
                if case == 'vT': scalef,Fieldf =    1005/1e06,['T',  'ta']
                if case == 'vq': scalef,Fieldf = 2260e03/1e06,['q', 'hus']
                if case == 'vv': scalef,Fieldf =       2/1e06,['V',  'va']
                if case == 'vu': scalef,Fieldf =       2/1e06,['U',  'ua']
		# Open stationary data
		levs_rey,lats_rey,lons_rey,vmfm,vafa,vfm,vm,fm=unpick('/mnt/climstorage/cian/Reynolds/%s.%s.0-1000.hPa.%s-%s.%s.p' % (Source,case,YearRange[0],YearRange[1],Season))
                # Extend lon axis and get dP
                lonres = np.abs(lons_rey[1]-lons_rey[0])
                dlon   = np.append(lons_rey,lons_rey[-1]+lonres)
                dlon   = np.append(lons_rey[0]-lonres,dlon)
		dP     = 100*(levs_rey[-1]-levs_rey[0])
		# Extract wavenumber section
                N             = 179 
                wavenumbers   = np.arange(0,N+1,1)
                VFWN          = np.zeros((0,len(lats),len(wavenumbers)))
                time          = [-999]
		# Meridional wind
		V   = vm
		V   = extendField(V)
		V   = interpolateND(V,dlon,lons,axis=2)
		V   = interpolateND(V,lats_rey,lats,axis=1)
		V   = interpolateND(V,levs_rey,levs,axis=0)
		Vzm = V.mean(axis=2)[:,:,np.newaxis]
		Vza = V - Vzm
		Vamp,Vang = waveN(Vza,N=N)
		# Specific humidity/Temperature
		F   = fm
		F   = extendField(F)
		F   = interpolateND(F,dlon,lons,axis=2)
		F   = interpolateND(F,lats_rey,lats,axis=1)
		F   = interpolateND(F,levs_rey,levs,axis=0)
		Fzm = F.mean(axis=2)[:,:,np.newaxis]
		Fza = F - Fzm
		Famp,Fang = waveN(Fza,N=N)
		# Meridional moisture flux decompsed into wave numbers
		VF0  = Vzm*Fzm
		VFwn = (Vamp*Famp/2.)*np.cos(Vang-Fang)
		VFwn = np.append(VF0,VFwn,axis=-1)
		# Take time mean and vertical integral
		VFwn = VFwn.mean(axis=0)[np.newaxis,:,:]*dP/9.81
		VFWN = np.append(VFWN,VFwn,axis=0)
                time = np.array(time)
                # Make NetCDF file
                makeNetCDF1(VFWN,time,wavenumbers,lats,'misc.',FileName,case)
        else:
                print '\nFile %s already exists, did nothing...\n' % (FileName)

Models = [g[28:] for g in glob.glob('/mnt/climstorage/cian/rcp85/*')]
lons    = np.arange(0,360,1)
lats    = np.arange(0,85+1,1)
levs    = np.arange(0,1000+50,50)

# Mass of girdbox on each level
dP_ = []
dP  = np.diff(levs)/2.
dP_.append(dP[0])
for i in range(len(dP)-1):
        dP_.append(dP[i]+dP[i+1])
dP_.append(dP[-1])
dP_ = np.array(dP_)
dM  = 100.*dP_/9.80665

#Model = str(sys.argv[1])
#Year  = int(sys.argv[2])
#Month = int(sys.argv[3])
#test(Model,Year,Month)

#Model = str(sys.argv[1])
#case  = str(sys.argv[2])
#phase = str(sys.argv[3])

# Decompose vq into vq(k)
#Amplitude_Phase_DailyFiles_vq(Model)

# Resynthesise v(k),q(k) or vq(k)
#resynthMonthlyFiles(Model,case)

# Decompose v or q into v(k) or q(k)
#Amplitude_Phase_DailyFiles(Model,case)

# Resynthesise v(k)q(k)
#resynthMontlyFiles_vq(Model)

# Resynthesise v(k),q(k) or vq(k)
#resynthMonthlyFiles_NS(Model,case,phase)

# Decompose vq into monthly mean zonal-mean v(k)q(k)
#waveNumberFile(Model,lons,lats,levs,dM,case)

# Decompose vq into 6xhourly zonal-mean v(k)q(k)
#waveNumberFile_Freq('ERAInt',lons,lats,levs,dM,'vq','dailymean')

# Decompose stationary vq into monthly mean zonal-mean v(k)q(k)
#YearRange = (int(sys.argv[3]),int(sys.argv[4]))
#waveNumberFile_stationary(Model,YearRange,'DJF',lons,lats,levs,case)


Models1 = Models[0:7]
Models2 = Models[7:14]
Models3 = Models[14:21]
Models4 = Models[21:28]

for Source in Models4:
	try:
		waveNumberFile_Freq(Source,lons,lats,levs,dM,'vq','dailymean')
		#waveNumberFile(Source,lons,lats,levs,dM,'vq')
		#Amplitude_Phase_DailyFiles(Source,'T')
		#Amplitude_Phase_DailyFiles_vq(Source)
		#resynthMonthlyFiles(Source,'q')
		#resynthMonthlyFiles_NS(Source,'vqk','-')
	except:
		pass

