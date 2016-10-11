#!/usr/bin/env python

# generate files for Idealized Ice Shelf problem.
# Gustavo Marques

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
import os

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Generate files for Idealized Ice Shelf problem.
      ''',
  epilog='Written by Gustavo Marques, Oct. 2016.')

  parser.add_argument('-nx', type=int, default=100,
      help='''The total number of grid points in the x direction (default = 100).''')

  parser.add_argument('-ny', type=int, default=300,
      help='''The total number of grid points in the y direction (default = 300).''')

  parser.add_argument('-nz', type=int, default=63,
      help='''Number of model layers (default = 63).''')
  
  parser.add_argument('-W', type=float, default=500.,
      help='''Domain width in the x direction (km). Default is 40.''')

  parser.add_argument('-L', type=float, default=1500.,
      help='''Domain lenght in the y direction (km). Default is 1.5E3''')

  parser.add_argument('-max_depth', type=float, default=3.0e3,
      help='''Maximum ocean depth (m). Default is 3E3.''')

  parser.add_argument('-heat_flux_polynya', type=float, default=-500.0,
      help='''Sensible heat flux into the ocean in the coastal polynya (W/m^2). Default is -500.''')

  parser.add_argument('-salt_flux_polynya', type=float, default=-4e-5,
      help='''Salt flux into the ocean in the coastal polynya (kg m^-2 s^-1). Default is -4E-5 (equivalent to gorwing 1 m of sea ice per year in the polynya region.)''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   nx = args.nx ; ny = args.ny; nz = args.nz
   L = args.L; W = args.W; D = args.max_depth

   dy = L/ny; dx = W/nx
   y = np.arange(dy/2.,L,dy)
   x = np.arange(dx/2.,W,dx)

   # create topography
   make_topo(x,y,args)

   # initial T/S
   make_ts(x,y,args)

   # create ice shelf
   #make_iceShelf(y,x) #not implemented yet because it needs scipy

   # create forcing
   make_forcing(x,y,args) 
   
   return

def make_ts(x,y,args):
   '''
   Extract T/S from WOA05 for a particulat lat. then interpolate results into ocean grid. 
   '''
   # all values at ~ 50 S (j=40)
   j=40
   temp = Dataset('WOA05_pottemp_salt.nc').variables['PTEMP'][:,:,j,:]
   salt = Dataset('WOA05_pottemp_salt.nc').variables['SALT'][:,:,j,:]
   depth = Dataset('WOA05_pottemp_salt.nc').variables['DEPTH'][:]
 
   # mean (time and space) values
   temp_mean = temp.mean(axis=2).mean(axis=0)
   salt_mean = salt.mean(axis=2).mean(axis=0)
   # model depth
   z = np.linspace(0,args.max_depth,args.nz) # positive downward
   # interpolate, notice that z and depth are normalized
   temp = np.interp(z/z.max(), depth/depth.max(), temp_mean)
   salt = np.interp(z/z.max(), depth/depth.max(), salt_mean)

   temp3D = np.zeros((1,args.nz,len(y),len(x)))
   salt3D = np.zeros((1,args.nz,len(y),len(x)))
   for i in range(len(x)):
      for j in range(len(y)):
          temp3D[0,:,j,i] = temp[:]
          salt3D[0,:,j,i] = salt[:]

   # create ncfiles

   # 1) The name of the z-space input file used to initialize
   #  the layer thicknesses, temperatures and salinities.
   name = 'ic_ts'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('DEPTH',args.nz)
   ncfile.createDimension('LAT',len(y))
   ncfile.createDimension('LON',len(x))
   ncfile.createDimension('TIME',1)

   # create variables
   LAT = ncfile.createVariable('LAT',np.dtype('double').char,('LAT'))
   LAT.units = 'degrees_north'
   LAT.long_name = 'h point nominal latitude'
   LAT.cartesian_axis = 'Y'
   LAT[:] = y[:]

   LON = ncfile.createVariable('LON',np.dtype('double').char,('LON'))
   LON.units = 'degrees_east'
   LON.long_name = 'h point nominal longitude'
   LON.cartesian_axis = 'X'
   LON[:] = x[:]

   TIME = ncfile.createVariable('TIME',np.dtype('double').char,('TIME'))
   TIME.units = 'days since 0001-01-01 00:00:00'
   TIME.calendar = 'noleap'
   TIME.cartesian_axis = 'T'
   TIME[0] = 0.0

   DEPTH = ncfile.createVariable('DEPTH',np.dtype('double').char,('DEPTH'))
   DEPTH.units = 'm'
   DEPTH.direction = -1
   DEPTH.cartesian_axis = 'Z'
   DEPTH[:] = z[:]

   PTEMP = ncfile.createVariable('PTEMP',np.dtype('float32').char,('TIME','DEPTH','LAT','LON'), fill_value = -1.e+34)
   PTEMP.missing_value = -1.e+34
   PTEMP[:] = temp3D[:] 

   SALT = ncfile.createVariable('SALT',np.dtype('float32').char,('TIME','DEPTH','LAT','LON'), fill_value = -1.e+34)  
   SALT.missing_value = -1.e+34
   SALT[:] = salt3D[:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # 2) The file from which the coordinate densities are read
   name = 'ts_ic_profile'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('Layer',args.nz)

   # create variables
   PTEMP = ncfile.createVariable('PTEMP',np.dtype('double').char,('Layer'))
   PTEMP.units = 'Celcius'
   PTEMP.long_name = 'Potential Temperature'
   PTEMP[:] = temp[:]

   SALT = ncfile.createVariable('SALT',np.dtype('double').char,('Layer'))
   SALT.units = 'PSU'
   SALT.long_name = 'Salinity'
   SALT[:] = salt[:]

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

def make_forcing(x,y,args):
   # wind parameters
   name = 'forcing'
   Ly = args.L # domain size km
   Lc = 200.0  # km
   Q0 = 10. # W/m^2
   Yt = 900.0  # km
   Lasf = 500.0 # km
   tau_acc = 0.2 # N/m^2
   tau_asf = -0.075 # N/m^2
   sponge = 100.0 # km
   # salt and salt fluxes
   salt_flux = args.salt_flux_polynya
   heat_flux = args.heat_flux_polynya

   nx = len(x); ny = len(y); nt =1 # time
   tau_x = np.zeros((nt,ny,nx))
   for j in range(ny):
      if y[j] <= (Yt-Lasf):
        tau = 0.0
      elif (y[j] > (Yt-Lasf) and y[j] <= Yt):
         tmp = np.pi*(Yt-y[j])/Lasf
         tau = tau_asf * np.sin(tmp)**2
      else:
         tmp = np.pi*(Yt-y[j])/(Ly-Yt)
         tau = tau_acc * np.sin(tmp)**2

      tau_x[0,j,:] = tau

   # heat 
   heat = np.zeros((nt,ny,nx))
   for j in range(ny):
    if y[j] < Lc:
       heat[0,j,:] = 0.0
    elif (y[j]>= Lc and y[j]<= Ly-sponge):
       tmp = (2*np.pi*(y[j]-Lc))/((Ly-Lc-sponge))
       heat[0,j,:] = -Q0 * np.sin(tmp)
       #heat[0,j,:] = -Q0 * np.cos(tmp)
    else:
       heat[0,j,:] = 0.0
   
   # evap, proxy for brine formation in polynyas
   brine = np.zeros((nt,ny,nx))
   for j in range(ny):
     if (y[j] >= 300.0 and y[j] <= 400.0):
        brine[0,j,:] = salt_flux
        heat[0,j,:] = heat_flux
 
   # create ncfile
   # open a new netCDF file for writing.
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('xh',nx)
   ncfile.createDimension('yh',ny)
   ncfile.createDimension('xq',nx)
   ncfile.createDimension('yq',ny)
   ncfile.createDimension('time',1)
   ncfile.createDimension('nv',2)
   
   # create variables
   xh = ncfile.createVariable('xh',np.dtype('double').char,('xh'))
   xh.units = 'km'
   xh.long_name = 'h point nominal longitude'
   xh.cartesian_axis = 'X'
   xh[:] = x[:]
   
   xq = ncfile.createVariable('xq',np.dtype('double').char,('xq'))
   xq.units = 'km'
   xq.long_name = 'q point nominal longitude'
   xq.cartesian_axis = 'X'
   xq[:] = x[:]

   yh = ncfile.createVariable('yh',np.dtype('double').char,('yh'))
   yh.units = 'km'
   yh.long_name = 'h point nominal latitude'
   yh.cartesian_axis = 'Y'
   yh[:] = y[:]
   
   yq = ncfile.createVariable('yq',np.dtype('double').char,('yq'))
   yq.units = 'km'
   yq.long_name = 'q point nominal latitude'
   yq.cartesian_axis = 'Y'
   yq[:] = y[:]

   time = ncfile.createVariable('time',np.dtype('double').char,('time'))
   time.long_name = 'time'
   time.units = 'days since 0001-01-01 00:00:00'
   time.cartesian_axis = 'T'
   time.calendar_type = 'NOLEAP'
   time.calendar = 'NOLEAP'
   time.bounds = 'time_bounds'
   time[0] = 0

   nv = ncfile.createVariable('nv',np.dtype('double').char,('nv'))   
   nv.long_name = 'vertex number'
   nv.units = 'none'
   nv.cartesian_axis = 'N'
   nv[:] = [1,2]

   SW = ncfile.createVariable('SW',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   SW.units = 'Watt meter-2'
   SW.missing_value = 1.e+20
   SW.long_name = 'surface_net_downward_shortwave_flux'
   SW[:] = 0.0

   LW = ncfile.createVariable('LW',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   LW.units = 'Watt meter-2'
   LW.missing_value = 1.e+20
   LW.long_name = 'surface_net_downward_longwave_flux'
   LW[:] = 0.0

   latent = ncfile.createVariable('latent',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   latent.units = 'Watt meter-2'
   latent.missing_value = 1.e+20
   latent.long_name = 'Latent heat flux into ocean due to fusion and evaporation'
   latent[:] = 0.0

   sensible = ncfile.createVariable('sensible',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   sensible.units = 'Watt meter-2'
   sensible.missing_value = 1.e+20
   sensible.long_name = 'surface_downward_sensible_heat_flux'
   sensible[:] = heat[:]

   evap = ncfile.createVariable('evap',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   evap.units = 'kilogram meter-2 second-1'
   evap.missing_value = 1.e+20
   evap.long_name = 'Evaporation at ocean surface (usually negative)'
   evap[:] = brine[:]

   taux = ncfile.createVariable('taux',np.dtype('float32').char,('time', 'yh', 'xq'), fill_value = 1.e+20)
   taux.units = 'Pascal'
   taux.missing_value = 1.e+20
   taux.long_name = 'Zonal Wind Stress'
   taux[:] = tau_x[:]   

   tauy = ncfile.createVariable('tauy',np.dtype('float32').char,('time', 'yq', 'xh'), fill_value = 1.e+20)
   tauy.units = 'Pascal'
   tauy.missing_value = 1.e+20
   tauy.long_name = 'Meridional Wind Stress'
   tauy[:] = 0.0

   ustar = ncfile.createVariable('ustar',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   ustar.units = 'meter second-1'
   ustar.missing_value = 1.e+20
   ustar.long_name = 'Surface friction velocity'
   ustar[:] = 0.0

   SST = ncfile.createVariable('SST',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = -1.e+34)
   SST.units = 'Celsius'
   SST.missing_value = -1.e+34
   SST.long_name = 'Sea Surface Temperature'
   SST[:] = 0.0

   SSS = ncfile.createVariable('SSS',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = -1.e+34)
   SSS.units = 'PSU'
   SSS.missing_value = -1.e+34
   SSS.long_name = 'Sea Surface Salinity'
   SSS[:] = 0.0

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
   
   return

def make_topo(x,y,args):
   # parameters
   name = 'topog'
   Hs = 500  # max depth at shelf
   Ys = 600.0e3 # width of shelf
   Ws = 70.0e3 # width of slope
   H = args.max_depth # max depth

   [X,Y] = np.meshgrid(x,y) 
   nx = len(x); ny = len(y)
   D = np.zeros((ny,nx))

   for j in range(ny):
      for i in range(nx):
          D[j,i] = Hs + 0.5 * (H-Hs) * (1.0 + np.tanh((Y[j,i]*1.0e3 - Ys)/Ws))

   H1 = 300.
   L1 = 250.0e3
   for j in range(ny):
      for i in range(nx):
          if Y[j,i]<= L1/1.0e3:
             D[j,i] = Hs - (H1/2.0) * (np.tanh((Y[j,i]*1.0e3 - L1/2.)/(L1/10.)) - 1.0) 
             if (i == 0 or i == nx-1): # side walls
                D[j,i] = 0.0 # land

   # open a new netCDF file for writing.
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   
   # create variables
   nx = ncfile.createVariable('nx',np.dtype('double').char,('nx'))
   nx.units = 'km'
   nx.description = 'x location of cell centers'
   nx.long_name = 'x location of cell centers'
   nx[:] = x[:]

   ny = ncfile.createVariable('ny',np.dtype('double').char,('ny'))
   ny.units = 'km'
   ny.description = 'y location of cell centers'
   ny.long_name = 'y location of cell centers'
   ny[:] = y[:]

   depth = ncfile.createVariable('depth',np.dtype('float32').char,('ny','nx'))
   depth.units = 'm'
   depth.description = 'depth at h points'
   depth.long_name = 'depth at h points'
   depth[:,:] = D[:,:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
   return



# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
