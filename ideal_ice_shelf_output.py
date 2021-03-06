#!/usr/bin/env python

# generate the ISOMIP + diagnostics
# Gustavo Marques, Aug. 2016

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os
#from computeOSF import computeOSF

class MyError(Exception):
  """
  Class for error handling
  """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Compute the isomip+ diagnostics and write them into a netCDF file. This script
      assumes that files/variables names are consistent with the default ISOMIP/GFDL
      options (see: https://github.com/NOAA-GFDL/MOM6-examples/tree/dev/master/ocean_on       ly/ISOMIP).
      ''',
  epilog='Written by Gustavo Marques, Aug. 2016.')

  parser.add_argument('-n', type=str, default='M1_exp5_1km_output',
      help='''The name of the experiment. Default is M1_exp5_1km_output.''')

  parser.add_argument('-icfile', type=str, default='IDEAL_IS_IC.nc',
      help='''The name of the ocean initial condition file. Default is IDEAL_IS_IC.nc.''')

  parser.add_argument('-geometry', type=str, default='ocean_geometry.nc',
      help='''The name of the ocean geometry file. Default is ocean_geometry.nc''')

  parser.add_argument('-isfile', type=str, default='MOM_Shelf_IC.nc',
      help='''The name of the ice shelf initial condition file. Default is MOM_Shelf_IC.nc.''')

  parser.add_argument('-month_file', type=str, default='ocean_month.nc',
      help='''The name of monthly mean file. Default is ocean_month.nc''')

  parser.add_argument('-month_z_file', type=str, default='ocean_month_z.nc',
      help='''The name of monthly mean file with z coordinate. Default is ocean_month_z.nc''')

  parser.add_argument('-prog_file', type=str, default='prog.nc',
      help='''The name of prognostic ocean file. Default is prog.nc''')

  parser.add_argument('-nx', type=float, default=250,
      help='''The number of grid points in the zonal direction (for the OUTPUT ncfile). Default is 250.''')

  parser.add_argument('-ny', type=float, default=500,
      help='''The number of grid points in the meridional direction (for the OUTPUT ncfile). Default is 500.''')

  parser.add_argument('-nz', type=float, default=400,
      help='''The number of grid points in the vertical direction (for the OUTPUT ncfile). Default is 400.''')

#  parser.add_argument('--test', action="store_true",
#      help='''Write 4D versions of elevation and overtuningStreamfunction so we can use gplot. ''')

  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., netcdf is generated/populated and functions for each diagnostic are called.
   """
   if name == 'Ocean0_COM_MOM6-LAYER':
      print("WARNING: exp. name has not been defined, using the default name Ocean0_COM_MOM6-LAYER.")


   # load essential variables
   # bedrock
   depth = Dataset(args.geometry).variables['D'][:]
   args.ny, args.nx = depth.shape
   # area under shelf
   shelf_area = Dataset(args.isfile).variables['shelf_area'][0,:,:]
   # base of STATIC ice shelf, which is ssh(t=0); make it positive
   ice_base = -Dataset(args.icfile).variables['ave_ssh'][0,:,:]
   ice_base[ice_base<1e-5] = 0.0
   #mask grounded ice and open ocean
   shelf_area = mask_grounded_ice(shelf_area,depth,ice_base)
   shelf_area = mask_ocean(shelf_area,shelf_area)

   # x,y, at tracer points
   x=Dataset(args.geometry).variables['geolon'][:]*1.0e3 # in m
   y=Dataset(args.geometry).variables['geolat'][:]*1.0e3 # in m

   # read entire time variable from montly mean
   # it might be useful to add the option of passing
   # just part of the data (e.g., last 6 months)
   time = Dataset(args.month_file).variables['time'][:]
   zi = Dataset(args.prog_file).variables['zi'][:]
   args.nz = len(zi)
   print len(zi)

   # create ncfile and zero fields. Each function below will corresponding values
   create_ncfile(name,time,args)

   # load additional variables
   ocean_area = Dataset(args.geometry).variables['Ah'][:]
   h = Dataset(args.month_file).variables['h'][:]
   h[h<1.0e-5] = 0.0
   melt = Dataset(args.month_file).variables['melt'][:]
   mass_flux = Dataset(args.month_file).variables['mass_flux'][:]
   # mask open ocean and grounded ice
   ocean_area = mask_grounded_ice(ocean_area,depth,ice_base)
   ocean_area = np.ma.masked_where(depth<10.,ocean_area)
   melt = mask_grounded_ice(melt,depth,ice_base)
   melt = mask_ocean(melt,shelf_area)
   # tracers
   salt = Dataset(args.month_file).variables['salt'][:]
   temp = Dataset(args.month_file).variables['temp'][:]
   rho = get_rho(salt,temp)

   # compute total volume
   total_volume(ocean_area,h)

   # compute mean sea level
   mean_sea_level(ocean_area, h, args)

   # compute mean melt rate
   mean_melt_rate(melt,shelf_area)

   # compute total melt flux
   total_melt_flux(mass_flux)

   # mean temp
   mean_tracer(ocean_area,h,temp,'meanTemperature','C')

   # mean salt
   mean_tracer(ocean_area,h,salt,'meanSalinity','PSU')

   # horizontal fields

   # bathymetry (negative)
   bathymetry = -mask_grounded_ice(depth,depth,ice_base)
   bathymetry.fill_value=0.0
#   saveXY(bathymetry,'bathymetry')

   # meltRate, already masked above
   melt = melt/(3600.*24*365) # in m/s
   saveXY(melt,'meltRate')

   # frictionVelocity
#   ustar_shelf = Dataset(args.month_file).variables['ustar_shelf'][:]
   # mask open ocean and grounded ice
#   ustar_shelf = mask_grounded_ice(ustar_shelf,depth,ice_base)
#   ustar_shelf = mask_ocean(ustar_shelf,shelf_area)
#   saveXY(ustar_shelf,'frictionVelocity')

   # thermalDriving
   thermal_driving = Dataset(args.month_file).variables['thermal_driving'][:]
   thermal_driving = mask_grounded_ice(thermal_driving,depth,ice_base)
   thermal_driving = mask_ocean(thermal_driving,shelf_area)
   saveXY(thermal_driving,'thermalDriving')

   # halineDriving
   haline_driving = Dataset(args.month_file).variables['haline_driving'][:]
   haline_driving = mask_grounded_ice(haline_driving,depth,ice_base)
   haline_driving = mask_ocean(haline_driving,shelf_area)
   saveXY(haline_driving,'halineDriving')

   # uBoundaryLayer
   u_ml = Dataset(args.month_file).variables['u_ml'][:]
   u_ml = mask_grounded_ice(u_ml,depth,ice_base)
   u_ml = mask_ocean(u_ml,shelf_area)
   saveXY(u_ml,'uBoundaryLayer')

   # vBoundaryLayer
   v_ml = Dataset(args.month_file).variables['v_ml'][:]
   v_ml = mask_grounded_ice(v_ml,depth,ice_base)
   v_ml = mask_ocean(v_ml,shelf_area)
   saveXY(v_ml,'vBoundaryLayer')

   #if (args.type == 'ocean3' or args.type == 'ocean4'):
   # iceDraft
   # will have to works this out when we run these cases
#   iceDraft = ice_base.copy()
#   iceDraft = mask_grounded_ice(iceDraft,depth,ice_base)
   #iceDraft[shelf_area == 0.] = 0.0
#   iceDraft.fill_value = 720.0; iceDraft = iceDraft.filled()
#   saveXY(-iceDraft,'iceDraft')

   # data from ocean_month_z
   temp_z = Dataset(args.month_z_file).variables['temp'][:]
   salt_z = Dataset(args.month_z_file).variables['salt'][:]

   # data at bottom most cell
   bottomTemperature = get_bottom_data(temp_z)
   bottomSalinity = get_bottom_data(salt_z)
   saveXY(bottomTemperature,'bottomTemperature')
   saveXY(bottomSalinity,'bottomSalinity')

   # interpolate tracers to get value at x = 520 km; y = 40 km
   temp_xz = 0.5 * (temp_z[:,:,0:-1,:] + temp_z[:,:,1::,:])
   salt_xz = 0.5 * (salt_z[:,:,0:-1,:] + salt_z[:,:,1::,:])
   temp_yz = 0.5 * (temp_z[:,:,:,0:-1] + temp_z[:,:,:,1::])
   salt_yz = 0.5 * (salt_z[:,:,:,0:-1] + salt_z[:,:,:,1::])

   # XZ y = 40 km (j=19 in the interpolated grid)
   #temperatureXZ = temp_xz[:,:,19,:]
   #salinityXZ = salt_xz[:,:,19,:]
   #saveXY(temperatureXZ,'temperatureXZ')
   #saveXY(salinityXZ,'salinityXZ')

   # YZ x = 520 km (i=99 in the interpolated grid)
   #temperatureYZ = temp_yz[:,:,:,99]
   #salinityYZ = salt_yz[:,:,:,99]
   #saveXY(temperatureYZ,'temperatureYZ')
   #saveXY(salinityYZ,'salinityYZ')

   # barotropic streamfunction
   uhbt = Dataset(args.month_file).variables['uhbt'][:]
   vhbt = Dataset(args.month_file).variables['vhbt'][:]
   # mask grouded region
   #uhbt = mask_grounded_ice(uhbt,depth,ice_base)
   #vhbt = mask_grounded_ice(vhbt,depth,ice_base)
   psi2D = get_psi2D(uhbt,vhbt)

   psi2D = mask_grounded_ice(psi2D,depth,ice_base)
   saveXY(psi2D,'barotropicStreamfunction')

   print('Done!')
   return

def get_psi2D(u,v):
    '''
    Loop in time and compute the barotropic streamfunction psi at h points.
    '''
    NT,NY,NX = u.shape
    uh = np.zeros(u.shape); vh = np.zeros(v.shape); psi = np.zeros(u.shape)
    # u and v at h points
    utmp = 0.5 * (u[:,:,0:-1] + u[:,:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
    vtmp = 0.5 * (v[:,0:-1,:] + v[:,1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    uh[:,:,1::] = utmp; uh[:,:,0] = 0.5*u[:,:,0] #u_i=1 = 0.5*u_(i=3/2)
    vh[:,1::,:] = vtmp; vh[:,0,:] = 0.5*v[:,0,:] #v_j=1 = 0.5*v_(j=3/2)
    for t in range(NT):
        # mask
        u_tmp = np.ma.masked_where(np.abs(uh[t,:]) > 1.0e30, uh[t,:])
        v_tmp = np.ma.masked_where(np.abs(vh[t,:]) > 1.0e30, vh[t,:])
        psi[t,:,:] = (-u_tmp.cumsum(axis=0) + v_tmp.cumsum(axis=1))*0.5

    return psi

def get_bottom_data(data):
    '''
    Return a 3D array (time,y,x) with data at the bottom most cell of each column.
    Data is already masked where there is no ocean.
    '''
    NT,NK,NY,NX = data.shape
    bottom_data = np.zeros((NT,NY,NX))
    for j in range(NY):
      for i in range(NX):
          # ice and bottom are static, so we can get indies at t = 0 since they dont change.
          ind = np.nonzero(data[0,:,j,i]!= data.fill_value)[0]
          if ind.any():
             bottom_data[:,j,i] = data[:,ind[-1],j,i]
          else:
             bottom_data[:,j,i] = np.nan

    return np.ma.masked_invalid(bottom_data)

def saveXY(var,varname):
   '''
   Save 2D (x,y) or 3D array (time,x,y) into the netCDF file.
   '''
   ncwrite(name,varname,var)
   return

def mean_tracer(area,h,var,varname,units):
   '''
   Compute tracer (temp,salt or any other tracer) averaged over the ocean volume.
   '''
   area = np.resize(area,h.shape)
   tracer = np.zeros(h.shape[0])
   for t in range(len(tracer)):
       vol = (area[t,:] * h[t,:])
       tracer[t] = (var[t,:] * vol).sum() /vol.sum()
       print(varname+' ('+units+') at t='+str(t)+' is:'+str(tracer[t]))

   ncwrite(name,varname,tracer)
   return

def total_melt_flux(mass):
    '''
    Compute the total mass flux of freshwater across the ice-oceaninterface, positive
    for melting and negative for freezing.
    '''
    total_mass = np.zeros(mass.shape[0])
    for t in range(len(total_mass)):
        total_mass[t] = mass[t,:].sum()
        print('Total melt flux (kg/s) at t='+str(t)+' is:'+str(total_mass[t]))

    ncwrite(name,'totalMeltFlux',total_mass)
    return

def mean_melt_rate(melt,area):
    '''
    Compute the melt rate, positive for melting negative for freezing, averaged over
    the ice shelf base (area>0).
    '''
    melt = melt/(3600*24.*365) # in m/s
    total_melt = np.zeros(melt.shape[0])
    for t in range(len(total_melt)):
        total_melt[t] = (melt[t,:] * area).sum() /area.sum()
        print('Averaged melt (m/s) at t='+str(t)+' is:'+str(total_melt[t]))

    ncwrite(name,'meanMeltRate',total_melt)
    return

def total_volume(area,h):
    area = np.resize(area,h.shape)
    vol = np.zeros(h.shape[0])
    for t in range(len(vol)):
        vol[t] = (area[t,:] * h[t,:]).sum()
        print('Total vol. (m^3) at t='+str(t)+' is:'+str(vol[t]))

    ncwrite(name,'totalOceanVolume',vol)
    return

def mean_sea_level(area,h,args):
    '''
    Compute mean sea level with respect to initial condition
    '''
    #area = np.resize(area,h.shape)
    h0 = Dataset(args.icfile).variables['h'][0,:]
    vol0 = (area * h0).sum()
    print('Domain area (m^2) is: '+str(area.sum()))
    print('Initial total volume (m^3) is: '+str(vol0))
    vol = np.zeros(h.shape[0])
    sl = np.zeros(h.shape[0])
    for t in range(len(vol)):
        vol[t] = (area * h[t,:]).sum()
        sl[t] = (vol[t] - vol0) / area.sum()
        print('Mean sea level (m) at t='+str(t)+' is:'+str(sl[t]))

    ncwrite(name,'meanSeaLevel',sl)
    return

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset(name+'.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def get_rho(S,T):
    '''
    Compute density with a linear EoS and using isomip+ coeffs.
    '''
    rho0=1027.51; Sref=34.2; Tref=-1.0
    alpha = 3.733e-5; beta = 7.843e-4
    rho = rho0 * (1 - alpha * (T-Tref) + beta * (S-Sref))
    return rho

def mask_ocean(data,area):
   """
   Mask open ocean. Works with 2D or 3D arrays.
   """
   if len(data.shape) == 2: # 2D array
     data = np.ma.masked_where(area==0,data)

   else: # 3D array
     NZ,NY,NX = data.shape
     area=np.resize(area,(NZ,NY,NX))
     data = np.ma.masked_where(area==0,data)

   return  data

def mask_grounded_ice(data,depth,base):
   """
   Mask regions where the ice shelf is grounded (base>=depth). Works with 2D or 3D arrays.
   """
   if len(data.shape) == 2: # 2D array
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here
   else: # 3D array
      NZ,NY,NX = data.shape
      base = np.resize(base,(NZ,NY,NX))
      depth = np.resize(depth,(NZ,NY,NX))
      data = np.ma.masked_where(base+1.0>=depth, data) # need +1 here

   return data

def create_ncfile(exp_name, ocean_time, args): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """

   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'.nc','w',format='NETCDF4')
   # dimensions
   nx = args.nx ; ny = args.ny ; nz = args.nz
   # create dimensions.
   #ncfile.createDimension('time', None)
   ncfile.createDimension('time', len(ocean_time))
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('nz',nz)

   # create variables, assign units and provide decription
   x = ncfile.createVariable('x',np.dtype('float32').char,('nx'))
   x.units = 'm'
   x.description = 'x location of cell centers'
   x.long_name = 'x location of cell centers'
   y = ncfile.createVariable('y',np.dtype('float32').char,('ny'))
   y.units = 'm'
   y.description = 'y location of cell centers'
   y.long_name = 'y location of cell centers'
   z = ncfile.createVariable('z',np.dtype('float32').char,('nz'))
   z.units = 'm'
   z.description = 'z location of cell centers'
   z.long_name = 'z location of cell centers'

   time = ncfile.createVariable('time',np.dtype('float32').char,('time'))
   time.units = 's'
   time.description = 'time since start of simulation'
   time.long_name = 'time since start of simulation'

   meanMeltRate = ncfile.createVariable('meanMeltRate',np.dtype('float32').char,('time'))
   meanMeltRate.units = 'm/s'; meanMeltRate.description = 'mean melt rate averaged over area of floating ice, positive for melting'
   totalMeltFlux = ncfile.createVariable('totalMeltFlux',np.dtype('float32').char,('time'))
   totalMeltFlux.units = 'kg/s'; totalMeltFlux.description = 'total flux of melt water summed over area of floating ice, positive for melting'

   totalOceanVolume = ncfile.createVariable('totalOceanVolume',np.dtype('float32').char,('time'))
   totalOceanVolume.units = 'm^3'; totalOceanVolume.description = 'total volume of ocean'

   meanSeaLevel = ncfile.createVariable('meanSeaLevel',np.dtype('float32').char,('time'))
   meanSeaLevel.units = 'm'; meanSeaLevel.description = 'mean sea level wrt IC'

   meanTemperature = ncfile.createVariable('meanTemperature',np.dtype('float32').char,('time'))
   meanTemperature.units = 'deg C'; meanTemperature.description = 'the potential temperature averaged over the ocean volume'
   meanSalinity = ncfile.createVariable('meanSalinity',np.dtype('float32').char,('time'))
   meanSalinity.units = 'PSU'; meanSalinity.description = 'the salinity averaged over the ocean volume'
   iceDraft = ncfile.createVariable('iceDraft',np.dtype('float32').char,('ny','nx'))
   bathymetry = ncfile.createVariable('bathymetry',np.dtype('float32').char,('ny','nx'))
   bathymetry.units = 'm'; bathymetry.description = 'elevation of the bathymetry'
   meltRate = ncfile.createVariable('meltRate',np.dtype('float32').char,('time','ny','nx'))
   meltRate.units = 'm/s'; meltRate.description = 'melt rate, positive for melting'
   thermalDriving = ncfile.createVariable('thermalDriving',np.dtype('float32').char,('time','ny','nx'))
   thermalDriving.units = 'deg C'; thermalDriving.description = 'thermal driving used in the melt calculation'
   halineDriving = ncfile.createVariable('halineDriving',np.dtype('float32').char,('time','ny','nx'))
   halineDriving.units = 'PSU';  halineDriving.description = 'haline driving used in the melt calculation'
   uBoundaryLayer = ncfile.createVariable('uBoundaryLayer',np.dtype('float32').char,('time','ny','nx'))
   uBoundaryLayer.units = 'm/s'; uBoundaryLayer.description = 'x-velocity in the boundary layer used to compute u*'
   vBoundaryLayer = ncfile.createVariable('vBoundaryLayer',np.dtype('float32').char,('time','ny','nx'))
   vBoundaryLayer.units = 'm/s'; vBoundaryLayer.description = 'y-velocity in the boundary layer used to compute u*'
   barotropicStreamfunction = ncfile.createVariable('barotropicStreamfunction',np.dtype('float32').char,('time','ny','nx'))
   barotropicStreamfunction.units = 'm^3/s'; barotropicStreamfunction.description = 'barotropic streamfunction'

   bottomTemperature = ncfile.createVariable('bottomTemperature',np.dtype('float32').char,('time','ny','nx'))
   bottomTemperature.units = 'deg C'; bottomTemperature.description = 'temperature in the bottom grid cell of each ocean column'
   bottomSalinity = ncfile.createVariable('bottomSalinity',np.dtype('float32').char,('time','ny','nx'))
   bottomSalinity.units = 'PSU'; bottomSalinity.description = 'salinity in the bottom grid cell of each ocean column'

   # write data to coordinate vars.
   lon = Dataset(args.icfile).variables['lonh'][:]
   lat = Dataset(args.icfile).variables['lath'][:]
   D = Dataset(args.geometry).variables['D'][:]
   W = (lon[-1]+(lon[1]-lon[0])*0.5)-(lon[0]-(lon[1]-lon[0])*0.5)
   L = (lat[-1]+(lat[1]-lat[0])*0.5)-(lat[0]-(lat[1]-lat[0])*0.5)
   max_depth = D.max()
   dz = max_depth/nz
   x[:] = lon[:]
   y[:] = lat[:]
   z[:] = -np.arange(0.5*dz,max_depth,dz)
   # time since start of simulation
   time[:] = ocean_time[:]/365.# in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')
   return

def write_ncfile(exp_name, data): # may add exp_type
   """
   This function writes the fields required for the isomip+ experiments into a pre-exiting netcdf file (exp_name). Different fields (given by the structure data) are added based on the type experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   print ('*** SUCCESS writing data into '+exp_name+'.nc!')
   return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
