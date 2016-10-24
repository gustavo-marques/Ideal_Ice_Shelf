#!/usr/bin/env python

# generate files for Idealized Ice Shelf problem.
# Gustavo Marques

import argparse
from netCDF4 import MFDataset, Dataset
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter
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

  parser.add_argument('-nx', type=int, default=200,
      help='''The total number of grid points in the x direction (default = 200).''')

  parser.add_argument('-ny', type=int, default=400,
      help='''The total number of grid points in the y direction (default = 400).''')

  parser.add_argument('-nz', type=int, default=63,
      help='''Number of model layers (default = 63).''')
  
  parser.add_argument('-W', type=float, default=500.,
      help='''Domain width in the x direction (km). Default is 500.''')

  parser.add_argument('-cshelf_lenght', type=float, default=400.,
      help='''Continental shelf lenght in the y direction (km). Default is 400.''')

  parser.add_argument('-slope_lenght', type=float, default=50.,
      help='''Continental shelf slope lenght in the y direction (km). Default is 50.''')

  parser.add_argument('-ice_shelf_lenght', type=float, default=150.,
      help='''Ice shelf lenght in the x direction (km). Default is 150.''')

  parser.add_argument('-polynya_major', type=float, default=75.,
      help='''Major axis (x direction) of polynya region (km). Default is 75.''')

  parser.add_argument('-tau_acc', type=float, default=0.2,
      help='''Max. wind stress in the ACC (N/m^2). Default is 0.2.''')

  parser.add_argument('-tau_asf', type=float, default=-0.075,
      help='''Max. wind stress in the ASF (N/m^2). Default is -0.075''')

  parser.add_argument('-katabatic_wind', type=float, default=0.05,
      help='''Max. katabatic wind stress (N/m^2). Default is 0.05''')

  parser.add_argument('-polynya_minor', type=float, default=75.,
      help='''Minor axis (y direction) of polynya region (km). Default is 75.''')

  parser.add_argument('-L', type=float, default=1000.,
      help='''Domain lenght in the y direction (km). Default is 1.5E3''')

  parser.add_argument('-max_depth', type=float, default=3.0e3,
      help='''Maximum ocean depth (m). Default is 3E3.''')

  parser.add_argument('-heat_flux_polynya', type=float, default=-200.0,
      help='''Sensible heat flux into the ocean in the coastal polynya (W/m^2). Default is -200.''')

  parser.add_argument('-salt_flux_polynya', type=float, default=1e-5,
      help='''Salt flux into the ocean in the coastal polynya (kg m^-2 s^-1). Default is 1E-5 (equivalent to gorwing ?? cm/s of sea ice in the polynya region.)''')

  parser.add_argument('-latent_heat_flux_polynya', type=float, default=0.0,
      help='''Latent heat flux in the coastal polynya (W/m^2). Default is 0.0)''')

  parser.add_argument('-coupled_run', help='''Generate all the files needed to run an ocean_SIS2 simulation.''', action="store_true")

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
   Ocean_Depth = make_topo(x,y,args)

   # create ice shelf
   make_ice_shelf(x,y,args) 

   if args.coupled_run:

     make_mosaic(x,y,Ocean_Depth,args) 

   # initial T/S
   make_ts(x,y,args)

   # create forcing
   make_forcing(x,y,args) 
   
   return

def make_ice_shelf(x,y,args):
   '''
   Here is a slightly different version but with constants
   H(x) = H0 *(Q0)^(1/4) / [Q0+4*H0^4*C^3*x]^(1/4)
   where H0 is ice thickness at the grounding line Q0 = U0*H0 is ice flux at the grounding line, C = rho*g*(1-rho/rho_w)/4/Bbar, rho is the ice density, rho_w is the sea-water density, Bbar is ice stiffness parameter.
   for the following parameters 
   rho = 917 kg/m^3
   rho_w = 1028 kg/m^3
   Bbar = 1.68e8
   C = 1.4440e-06
   U0= 700 m/yr = 2.2182e-05 m/s
   H0 = 1500 m
   Q0 = 0.0333 m^2/s
   rho = 917. #kg/m^3
   rho_w = 1028. #kg/m^3
   Bbar = 1.68e8
   '''
   x = x * 1.0e3 # im m
   y = y * 1.0e3 # im m
   C = 1.4440e-06
   H0 = 1500.0 #m
   Q0 = 0.03327 #m^2/s
   ISL = args.ice_shelf_lenght * 1.0e3 # ice shelf lenght
   gp = 20.e3 # grouding line position
   dy = y[1]-y[0]
   dx = x[1]-x[0]

   h =  H0 *(Q0)**(1./4.) / (Q0+100*H0**4*C**3*(y-gp))**(1./4.)
   h[y>=ISL] = 0.0
   h[y<gp] = H0

   # smooth
   h_smooth = gaussian_filter(h,4)
   h_smooth[y>=ISL] = 0.0
   h_smooth[y<gp] = H0

   #plt.plot(y,h,'k',y,h_smooth,'r')
   #plt.show()

   area_tmp = np.ones(args.ny) * dx * dy
   area_tmp[h_smooth == 0.0] = 0.0

   # Create a netcdf horizontal ocean-grid file
   name = 'IC_IS'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('nx',args.nx)
   ncfile.createDimension('ny',args.ny)
   thick = ncfile.createVariable('thick','double',('ny','nx',))
   thick.units = 'm'
   thick.standard_name =  'ice shelf thickness'
   area = ncfile.createVariable('area','double',('ny','nx',))
   area.units = 'm2'
   area.standard_name =  'ice shelf area'  

   # write into nc file
   for i in range(args.nx):
        thick[:,i] = h_smooth[:]
        area[:,i] = area_tmp[:]

   ncfile.sync()
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!') 

def make_mosaic(x,y,Ocean_Depth,args):
   '''
   Create files used in coupled ice/ocean runs.
   '''
   ni = args.nx;   nj = args.ny
   snj, sni = 2*nj, 2*ni
   dxo = (x[1] - x[0]) * 1.0e3 # in m
   dyo = (y[1] - y[0]) * 1.0e3 # in m

   # Number of land points used nl GM: There has to be at least ONE land point
   nl=0
   #Define land points:
   for j in range(nj):
      for i in range(ni):
          if Ocean_Depth[j,i]==0:
              nl=nl+1

   # Using a supergrid refinement of 2, the grid lengths, area and coordinates are:
   dx = (0.5 * dxo) * np.ones((snj+1,sni))
   dy = (0.5 * dyo) * np.ones((snj,sni+1))
   area = (0.25 * (dxo * dyo)) * np.ones((snj,sni))
   x = np.zeros((snj+1,sni+1))
   x[:,1:] = np.cumsum( dx, axis=-1 )
   y = np.zeros((snj+1,sni+1))
   y[1:,:] = np.cumsum( dy, axis=-2 )

   # Create a netcdf horizontal ocean-grid file
   name = 'ocean_hgrid'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('nx',sni)
   ncfile.createDimension('ny',snj)
   ncfile.createDimension('nxp',sni+1)
   ncfile.createDimension('nyp',snj+1)
   ncfile.createDimension('string',255)
   dx_h = ncfile.createVariable('dx','f8',('nyp','nx',))
   dx_h.units = 'm'
   dy_h = ncfile.createVariable('dy','f8',('ny','nxp',))
   dy_h.units = 'm'
   area_h = ncfile.createVariable('area','f8',('ny','nx',))
   area_h.units = 'm2'
   x_h = ncfile.createVariable('x','f8',('nyp','nxp',))
   x_h.units = 'm'
   y_h = ncfile.createVariable('y','f8',('nyp','nxp',))
   y_h.units = 'm'
   tile = ncfile.createVariable('tile','c',('string',)) 
   dx_h[:,:] = dx[:,:]
   dy_h[:,:] = dy[:,:]
   area_h[:,:] = area[:,:]
   x_h[:,:] = x[:,:]
   y_h[:,:] = y[:,:]
   tile[:] = '\000' * 255
   tile[:5] = 'tile1'
   set_string(tile,'tile1')
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # create ocean_mask file
   name = 'ocean_mask'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('nx',ni)
   rg.createDimension('ny',nj)
   mask = rg.createVariable('mask','f4',('ny','nx',))
   mask[:,:] = 1.
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # Create a mosaic description for the grid file
   name = 'ocean_mosaic'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('ntiles',1)
   rg.createDimension('string',255)
   mosaic = rg.createVariable('mosaic','c',('string',))
   mosaic.standard_name = 'grid_mosaic_spec'
   mosaic.children = 'contacts'
   mosaic.grid_descriptor = ''
   gridlocation = rg.createVariable('gridlocation','c',('string',))
   gridlocation.standard_name = 'grid_file_location'
   gridfiles = rg.createVariable('gridfiles','c',('ntiles','string',))
   gridtiles = rg.createVariable('gridtiles','c',('ntiles','string',))
   rg.grid_version = '0.2'
   # Fill in data
   mosaic[:] = '\000' * 255
   mosaic[:12] = 'ocean_mosaic'
   gridlocation[:] = '\000' * 255
   gridlocation[:2] = './'
   gridfiles[:] = '\000' * 255
   gridfiles[0,:14] = 'ocean_hgrid.nc'
   gridtiles[:] = '\000' * 255
   gridtiles[0,:5] = 'tile1'
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   name = 'atmos_mosaic_tile1Xland_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',nl)  #It is unclear whether this works when nl=0. It does work for nl>0
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two'))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two'))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells'))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]

   xgrid_area[:] = dxo * dyo
   tile1_distance[:] = 0.
   tile2_distance[:] = 0.
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   name = 'atmos_mosaic_tile1Xocean_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]
   
   xgrid_area[:] = dxo * dyo
   count=-1
   for j in range(nj):
       for i in range(ni):
          if Ocean_Depth[j,i]!=0:
            count=count+1
            tile1_cell[count] = [i,j]
            tile2_cell[count] = [i,j]
            tile1_distance[count] = [0,0]
            tile2_distance[count] = [0,0]
            xgrid_area[count] = dxo * dyo

   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!') 

   name = 'land_mosaic_tile1Xocean_mosaic_tile1'
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
   rg.createDimension('two',2)
   contact = rg.createVariable('contact','c',('string',))
   contact.standard_name = 'grid_contact_spec'
   contact.contact_type = 'exchange'
   contact.parent1_cell = 'tile1_cell'
   contact.parent2_cell = 'tile2_cell'
   contact.xgrid_area_field = 'xgrid_area'
   contact.distant_to_parent1_centroid = 'tile1_distance'
   contact.distant_to_parent2_centroid = 'tile2_distance'
   tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
   tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
   tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
   tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
   xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
   xgrid_area.standard_name = 'exchange_grid_area'
   xgrid_area.units = 'm2'
   tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
   tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
   tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
   tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
   rg.grid_version = '0.2'
   # Fill in data
   contact[:] = '\000' * 255
   contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
   for i in range(nl):
      tile1_cell[i,:] = [ni,nj]
      tile2_cell[i,:] = [ni,nj]
   
   xgrid_area[:] = dxo * dyo
   count=-1
   for j in range(nj):
      for i in range(ni):
        if Ocean_Depth[j,i]!=0:
            count=count+1
            tile1_cell[count] = [i,j]
            tile2_cell[count] = [i,j]
            tile1_distance[count] = [0,0]
            tile2_distance[count] = [0,0]
            xgrid_area[count] = dxo * dyo
  
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')


   name = 'mosaic' # sometimes grid_spec is called mosaic.nc
   rg = Dataset(name+'.nc','w')
   rg.createDimension('string',255)
   rg.createDimension('nfile_aXo',1) # -1 is for a single land point
   rg.createDimension('nfile_aXl',1) # -1 is for a single land point
   rg.createDimension('nfile_lXo',1) # -1 is for a single land point
   atm_mosaic_dir = rg.createVariable('atm_mosaic_dir','c',('string',))
   atm_mosaic_dir.standard_name = 'directory_storing_atmosphere_mosaic'
   atm_mosaic_file = rg.createVariable('atm_mosaic_file','c',('string',))
   atm_mosaic_file.standard_name = 'atmosphere_mosaic_file_name'
   atm_mosaic = rg.createVariable('atm_mosaic','c',('string',))
   atm_mosaic.standard_name = 'atmosphere_mosaic_name'
   lnd_mosaic_dir = rg.createVariable('lnd_mosaic_dir','c',('string',))
   lnd_mosaic_dir.standard_name = 'directory_storing_land_mosaic'
   lnd_mosaic_file = rg.createVariable('lnd_mosaic_file','c',('string',))
   lnd_mosaic_file.standard_name = 'land_mosaic_file_name'
   lnd_mosaic = rg.createVariable('lnd_mosaic','c',('string',))
   lnd_mosaic.standard_name = 'land_mosaic_name'
   ocn_mosaic_dir = rg.createVariable('ocn_mosaic_dir','c',('string',))
   ocn_mosaic_dir.standard_name = 'directory_storing_ocean_mosaic'
   ocn_mosaic_file = rg.createVariable('ocn_mosaic_file','c',('string',))
   ocn_mosaic_file.standard_name = 'ocean_mosaic_file_name'
   ocn_mosaic = rg.createVariable('ocn_mosaic','c',('string',))
   ocn_mosaic.standard_name = 'ocean_mosaic_name'
   ocn_topog_dir = rg.createVariable('ocn_topog_dir','c',('string',))
   ocn_mosaic_dir.standard_name = 'directory_storing_ocean_topog'
   ocn_topog_file = rg.createVariable('ocn_topog_file','c',('string',))
   ocn_topog_file.standard_name = 'ocean_topog_file_name'
   aXo_file = rg.createVariable('aXo_file','c',('nfile_aXo','string',))
   aXo_file.standard_name = 'atmXocn_exchange_grid_file'
   aXl_file = rg.createVariable('aXl_file','c',('nfile_aXl','string',))
   aXl_file.standard_name = 'atmXlnd_exchange_grid_file'
   lXo_file = rg.createVariable('lXo_file','c',('nfile_lXo','string',))
   lXo_file.standard_name = 'lndXocn_exchange_grid_file'
   #Global attributes
   rg.grid_version = '0.2'
   rg.code_version = "$Name:  $"
   rg.history = "/net2/aja/workspace/MOM6-examples/ice_ocean_SIS2/OM4_025/preprocessing/fre_nctools/tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ocean_topog.nc"
   # Fill in data
   atm_mosaic_dir[:] = '\000' * 255
   atm_mosaic_dir[:2] = './'
   atm_mosaic_file[:] = '\000' * 255
   atm_mosaic_file[:15] = 'ocean_mosaic.nc'
   atm_mosaic[:] = '\000' * 255
   atm_mosaic[:12] = 'atmos_mosaic'
   lnd_mosaic_dir[:] = '\000' * 255
   lnd_mosaic_dir[:2] = './'
   lnd_mosaic_file[:] = '\000' * 255
   lnd_mosaic_file[:15] = 'ocean_mosaic.nc'
   lnd_mosaic[:] = '\000' * 255
   lnd_mosaic[:11] = 'land_mosaic'
   ocn_mosaic_dir[:] = '\000' * 255
   ocn_mosaic_dir[:2] = './'
   ocn_mosaic_file[:] = '\000' * 255
   ocn_mosaic_file[:15] = 'ocean_mosaic.nc'
   ocn_mosaic[:] = '\000' * 255
   ocn_mosaic[:12] = 'ocean_mosaic'
   ocn_topog_dir[:] = '\000' * 255
   ocn_topog_dir[:2] = './'
   ocn_topog_file[:] = '\000' * 255
   ocn_topog_file[:14] = 'ocean_topog.nc'
   aXo_file[:,:] = '\000' * 255
   aXo_file[:,:40] = 'atmos_mosaic_tile1Xocean_mosaic_tile1.nc'
   aXl_file[:,:] = '\000' * 255
   aXl_file[:,:39] = 'atmos_mosaic_tile1Xland_mosaic_tile1.nc'
   lXo_file[:,:] = '\000' * 255
   lXo_file[:,:39] = 'land_mosaic_tile1Xocean_mosaic_tile1.nc'
   rg.close()
   print ('*** SUCCESS creating '+name+'.nc!')

def make_ts(x,y,args):
   '''
   Extract T/S from WOA05 for a particulat lat. then interpolate results into ocean grid. 
   '''
   # all values at ~ 58 S (j=27)
   # 185 lon (Off Ross Sea)
   # just winter time
   j=27
   temp = Dataset('WOA05_pottemp_salt.nc').variables['PTEMP'][6,:,j,185]
   salt = Dataset('WOA05_pottemp_salt.nc').variables['SALT'][6,:,j,185]
   depth = Dataset('WOA05_pottemp_salt.nc').variables['DEPTH'][:]
   # remove skin layer
   temp[-3::] = temp[-4]; salt[-3::] = salt[-4]
   # remove mask in the bottom
   temp[depth<=100.] = 0.0; salt[depth<=100.] = 33.94
   # smooth
   salt = gaussian_filter(salt,2)
   temp = gaussian_filter(temp,2)

   # mean (time and space) values
   #temp_mean = temp.mean(axis=2).mean(axis=0)
   #salt_mean = salt.mean(axis=2).mean(axis=0)
   
   # model depth
   z = np.linspace(0,args.max_depth,args.nz) # positive downward
   # interpolate, notice that z and depth are normalized
   # cubic
   f1 = interp1d(depth, temp)
   temp_int = f1(z)
   f1 = interp1d(depth, salt)
   salt_int = f1(z)
   # linear
   #temp_int = np.interp(z/z.max(), depth/depth.max(), temp)
   #salt_int = np.interp(z/z.max(), depth/depth.max(), salt)

   #plt.figure()
   #plt.plot(temp,-depth,temp_int,-z)
   #plt.show()

   #plt.figure()
   #plt.plot(salt,-depth,salt_int,-z)
   #plt.show()

   temp3D = np.zeros((1,args.nz,len(y),len(x)))
   salt3D = np.zeros((1,args.nz,len(y),len(x)))
   for i in range(len(x)):
      for j in range(len(y)):
          temp3D[0,:,j,i] = temp_int[:]
          salt3D[0,:,j,i] = salt_int[:]

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
   PTEMP[:] = temp_int[:]

   SALT = ncfile.createVariable('SALT',np.dtype('double').char,('Layer'))
   SALT.units = 'PSU'
   SALT.long_name = 'Salinity'
   SALT[:] = salt_int[:]

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

def make_forcing(x,y,args):
   # wind parameters
   name = 'forcing'
   Ly = args.L # domain size km
   W = args.W # domain width km
   CSL = args.cshelf_lenght  # km
   Q0 = 10. # W/m^2
   Yt = 550.0  # km
   Lasf = 300.0 # km
   tau_acc = args.tau_acc # default is 0.2 # N/m^2
   tau_asf = args.tau_asf # default is -0.075 # N/m^2
   katabatic_wind = args.katabatic_wind # def is 0.05 N/m2
   sponge = 50.0 # km
   # polynya salt and salt fluxes
   major = args.polynya_major
   minor = args.polynya_minor
   area = np.pi * major * minor * 0.5
   print 'Polynya area is (km^2):', area
   ISL = args.ice_shelf_lenght
   salt_flux = args.salt_flux_polynya 
   heat_flux = args.heat_flux_polynya
   latent_flux = args.latent_heat_flux_polynya
   ice_form = ((np.abs(salt_flux)*3600*24) * area * 1.0e6* 1.0e3/34.) * 1.0e-3 / (area * 1.0e6) 
   print 'Sea ice formation (m/day)', ice_form
   ice_form = ice_form/area
   nx = len(x); ny = len(y); nt =1 # time
   tau_x = np.zeros((nt,ny,nx))
   tau_y = np.zeros((nt,ny,nx))
   evaporation = np.zeros((nt,ny,nx))
   salt = np.zeros((nt,ny,nx))
   heat = np.zeros((nt,ny,nx))
   latent = np.zeros((nt,ny,nx))

   # tau x
   for j in range(ny):
      if y[j] <= (Yt-Lasf):
        tauX = 0.0
      elif (y[j] > (Yt-Lasf) and y[j] <= Yt):
         tmp = np.pi*(Yt-y[j])/Lasf
         tauX = tau_asf * np.sin(tmp)**2
      else:
         tmp = np.pi*(Yt-y[j])/(Ly-Yt)
         tauX = tau_acc * np.sin(tmp)**2
      tau_x[0,j,:] = tauX

   # tau y (decays linearly away from the ice shelf front)
   for j in range(ny):
      if y[j] >= ISL and y[j] <= Yt:
         tau_y[0,j,:] = (-katabatic_wind/(Yt-ISL)) * (y[j]-ISL) + katabatic_wind

   # heat/brine
   tmp1 = ISL + minor
   for j in range(ny):
    if y[j] < tmp1:
       heat[0,j,:] = 0.0
    elif y[j] >= tmp1 and y[j]< Ly-sponge:
       tmp2 = (2*np.pi*(y[j]-tmp1))/((Ly-tmp1-sponge))
       heat[0,j,:] = -Q0 * np.sin(tmp2)
       #heat[0,j,:] = -Q0 * np.cos(tmp)
    else:
       heat[0,j,:] = 0.0

   # evap, proxy for brine formation in polynyas
   for i in range(nx):
     for j in range(ny):
       if (x[i] - W/2. >= -major and x[i] - W/2. <= major):
          if (y[j]-ISL >= 0.) and (y[j]-ISL <= (minor * np.abs((1-(x[i]-W/2.)**2/major**2)**(1/2.)))):
            evaporation[0,j,i] = salt_flux
            salt[0,j,i] = salt_flux
            heat[0,j,i] = heat_flux
            latent[0,j,i] = latent_flux

   # remove heat gain
   for j in range(ny):
       if y[j]>ISL and (heat[0,j,i] > 0.): heat[0,j,:] = 0.0 
   # create ncfile
   # open a new netCDF file for writing.
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('xh',nx)
   ncfile.createDimension('yh',ny)
   ncfile.createDimension('xq',nx)
   ncfile.createDimension('yq',ny)
   ncfile.createDimension('time',None)
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

   if args.coupled_run:
     u_flux = ncfile.createVariable('u_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     u_flux.units = 'Pa'
     u_flux.missing_value = 1.e+20
     u_flux.long_name = 'i-direction wind stress'
     u_flux[:] = -tau_x[:] # change sign in ice

     v_flux = ncfile.createVariable('v_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     v_flux.units = 'Pa'
     v_flux.missing_value = 1.e+20
     v_flux.long_name = 'j-direction wind stress'
     v_flux[:] = -tau_y[:] # change sign in ice

     t_flux = ncfile.createVariable('t_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     t_flux.units = 'Watt meter-2'
     t_flux.missing_value = 1.e+20
     t_flux.long_name = 'Sensible heat flux'
     t_flux[:] = -heat[:] # change sign in ice

     # latent heat
     lt = ncfile.createVariable('latent',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     lt.units = 'Watt meter-2'
     lt.missing_value = 1.e+20
     lt.long_name = 'Latent heat flux'
     lt[:] = -latent # change sign in ice

     salt_flux = ncfile.createVariable('salt_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     salt_flux.units = 'kg/(m^2 s)'
     salt_flux.missing_value = 1.e+20
     salt_flux.long_name = 'salt flux'
     salt_flux[:] = -salt[:] # change sign in ice

   else:
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
     evap[:] = evaporation[:]

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
   name = 'ocean_topog'
   Hs = 500  # max depth at shelf
   Ys = args.cshelf_lenght * 1.0e3 # len of shelf in m
   Ws = args.slope_lenght * 1.0e3 # len of slope in m
   H = args.max_depth # max depth
   W = args.W  # domain width in km
   [X,Y] = np.meshgrid(x,y) 
   nx = len(x); ny = len(y)
   D = np.zeros((ny,nx))

   for j in range(ny):
      for i in range(nx):
          D[j,i] = Hs + 0.5 * (H-Hs) * (1.0 + np.tanh((Y[j,i]*1.0e3 - Ys)/Ws))

   H1 = 300. # # depth increase under ice shelf
   L1 = 200.0e3 # size of ice shelf cavity
   Wl = 50.0 # widht of land region next to ice shelf in km
   for j in range(ny):
      for i in range(nx):
          if Y[j,i]<= L1/1.0e3:
             D[j,i] = Hs - (H1/2.0) * (np.tanh((Y[j,i]*1.0e3 - L1/2.)/(L1/10.)) - 1.0) 
             if (X[j,i] <= Wl or X[j,i] >= W-Wl): # side walls
                D[j,i] = 0.0 # land

   # open a new netCDF file for writing.
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('ntiles',1)
 
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
   return D

# A helper function to use when writing strings in a netcdf file
def set_string(variable, value):
   """Sets "variable" to "value" padded with blanks where
   "variable" is a netcdf variable object and "value" is a string."""
   variable[:] = '\000' * variable.shape[0]
   variable[:len(value)] = value
   return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
