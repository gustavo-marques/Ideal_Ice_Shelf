#!/usr/bin/env python

# generate files for Idealized Ice Shelf problem.
# Gustavo Marques

#import matplotlib
#matplotlib.use('Agg')
import argparse
from netCDF4 import MFDataset, Dataset
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import wright_eos as eos
import gsw
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

  parser.add_argument('-nx', type=int, default=1200,
      help='''The total number of grid points in the x direction (default = 1200).''')

  parser.add_argument('-ny', type=int, default=2400,
      help='''The total number of grid points in the y direction (default = 2400).''')

  parser.add_argument('-nz', type=int, default=63,
      help='''Number of model layers (default = 63).''')
  
  parser.add_argument('-W', type=float, default=600.,
      help='''Domain width in the x direction (km). Default is 600.''')

  parser.add_argument('-L', type=float, default=1200.,
      help='''Domain lenght in the y direction (km). Default is 1.2E3''')

  parser.add_argument('-max_depth', type=float, default=3.0e3,
      help='''Maximum ocean depth (m). Default is 3E3.''')

  parser.add_argument('-min_depth', type=float, default=500.0,
      help='''Minimum ocean depth (m). Default is 500.''')

  parser.add_argument('-debug', help='''Adds prints and plots to help debug.''', action="store_true")

  parser.add_argument('-homogeneous_ts', help='''Make the initial T/S homogeneous in the horizontal.''', action="store_true")

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

   # create dir. to place figures
   os.system('mkdir PNG')

   # create ice shelf
   #make_ice_shelf(x,y,args)

   # create topography
   Ocean_Depth = make_topo(x,y,args)

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
   H0 = 2000.0 #m
   Q0 = 0.03327 #m^2/s
   gp = 20.e3 # grouding line position
   dy = y[1]-y[0]
   dx = x[1]-x[0]
   Lice = 200.0e3
   h =  H0 *(Q0)**(1./4.) / (Q0+100*H0**4*C**3*(y-gp))**(1./4.)
   h[y<gp] = H0
   h[y>Lice] = 0.0 
   # smooth
   h_smooth = gaussian_filter(h,2)
   h_smooth[y<gp] = H0
   h_smooth[h_smooth<50.0] = 0.0
   # find meridional lenght of ice shelf
   tmp = np.nonzero(h_smooth==0.0)[0][0]
   args.ISL = x[tmp] / 1.0e3
   print 'Ice shelf meridional lenght is (km):',x[tmp] / 1.0e3
   
   if args.debug:
      plt.plot(y,h,'k',y,h_smooth,'r')
      plt.savefig('PNG/ice_shelf_profile.png')

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

def get_profile(t,i,j,var,depth,z,args,vname,take_mean = True):
   '''
   Get data (var), remove mask, smooth profile then interpolate to z grid.
   '''
   if take_mean:
      data = np.mean(Dataset('WOA05_pottemp_salt.nc').variables[var][:,:,j,i], axis=0)
      
   else:
      data = Dataset('WOA05_pottemp_salt.nc').variables[var][t,:,j,i]

   # replace mask value with last good value
   tmp = np.nonzero(data.mask==False)[0][-1]
   data[tmp+1::] = data[tmp]
   # smooth  
   data = gaussian_filter(data,2) # filter std 2
   # interpo
   f1 = interp1d(depth, data)
   if args.debug:
      plt.figure()
      plt.plot(data,-depth,'b',f1(z),-z,'rx')
      plt.title(vname)
      plt.savefig('PNG/'+vname+'.png')

   return f1(z)

def set_freezing_temp(T,S,z,y,args):
    '''
    Set the temp. in the upper 50 m to the surface freezing point, then smooth data
    '''
    d = np.nonzero(z<=50.)[0]
    for k in range(len(d)):
        for j in range(args.ny):
           # set T to freezing point just near the coast
           if y[j] <= (args.ISL + args.cshelf_lenght):
              T[0,k,j,:] = eos.tfreeze(S[0,k,j,:],1.0e5) 

    #for j in range(args.ny):
    #  for i in range(args.nx):
    #     T[0,d,j,i] = gaussian_filter(T[0,d,j,i],1) # filter std 1

    # smooth T and S in the cross slope direction
    #for k in range(args.nz):
    #   for j in range(args.ny):
    #       T[0,k,:,i] = gaussian_filter(T[0,k,:,i],2)
    #       S[0,k,:,i] = gaussian_filter(T[0,k,:,i],2)

    return T,S

def make_ts(x,y,args):
   '''
   Extract T/S from WOA05 for a particulat lat. then interpolate results into ocean grid. 
   '''
   # climatology depth
   depth = Dataset('WOA05_pottemp_salt.nc').variables['DEPTH'][:]
   # model depth
   z = np.linspace(0,args.max_depth,args.nz) # positive downward
   i =  176; t = 6; j = 21
   # northern wall conditions (-65)
   temp_north = get_profile(t,i,j,'PTEMP',depth,z,args,'TempNorth')
   salt_north = get_profile(t,i,j,'SALT',depth,z,args,'SaltNorth')

   # 3d fields
   temp3D = np.zeros((1,args.nz,len(y),len(x)))
   salt3D = np.zeros((1,args.nz,len(y),len(x)))

   dist = [0,args.max_depth]
   f1 = interp1d(dist,[0, 2])
   dummy_t = f1(z)
   f1 = interp1d(dist,[33.,35.])
   dummy_s = f1(z)
   sigma2 = eos.wright_eos(dummy_t,dummy_s,2.0e7)
   for i in range(args.nx):
      for j in range(args.ny):
          temp3D[0,:,j,i] = dummy_t[:]
          salt3D[0,:,j,i] = dummy_s[:]

   layers = Dataset('GOLD_IC.2010.11.15.nc').variables['Layer'][:] # used in the global run with sig2

   print 'sigma2 min/max',sigma2.min(),sigma2.max()
   print 'Layers min/max', layers.min(),layers.max()
 
   # create ncfiles
   # vertical coord. (layers)
   name = 'vertical_layers'
   ncfile = Dataset(name+'.nc','w')
   ncfile.createDimension('Layer',args.nz)
   Layer = ncfile.createVariable('Layer',np.dtype('double').char,('Layer'))
   Layer.units = 'kg/m^3'
   Layer[:] = layers[:]

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
   
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

def make_forcing(x,y,args):
   # forcing parameters
   Ly = args.L # domain size km
   W = args.W # domain width km
   tau_acc = 0.1
   sponge = 100.0 # km
   nx = len(x); ny = len(y); nt =1 # time
   tau_x = np.zeros((nt,ny,nx)) # ice forcing
   heat = np.ones((nt,ny,nx))
   # wind and heat
   # x-dir
   for j in range(ny):
      if y[j] >= Ly/4. and y[j] <= 3*Ly/4.:
        tmp2 = (y[j] - Ly/4.)/(Ly/2.)
        tau_x[0,j,:] = tau_acc * np.sin(tmp2 * np.pi)**2
        
   # e-folding scale is 150 km
   for j in range(ny):
      if y[j]<=Ly/4.:
         heat[0,j,:] = -50.
      elif y[j]<=Ly/2.:
         heat[0,j,:] = -50. * np.exp(-(y[j]-Ly/2.)/90.) - 10.0
      else:
         heat[0,j,:] = -10.0
   
   # plots
   if args.debug:
      plt.figure()
      plt.subplot(211)
      plt.plot(y,tau_x[0,:,1])
      plt.title('taux'); plt.xlabel('y [km]'); plt.ylabel('Pa')
      plt.subplot(212)
      plt.contourf(x,y,tau_x[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/taux.png')

      plt.figure()
      plt.subplot(211)
      plt.plot(y,heat[0,:,1])
      plt.title('Sensible Heat'); plt.xlabel('y [km]'); plt.ylabel('W/m^2')
      plt.subplot(212)
      plt.contourf(x,y,heat[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/sensible_heat.png')

   # create ncfile
   # used when forcing is applied in the sea ice
   name = 'forcing'
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

   u_flux = ncfile.createVariable('u_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   u_flux.units = 'Pa'
   u_flux.missing_value = 1.e+20
   u_flux.long_name = 'i-direction wind stress'
   u_flux[:] = -tau_x[:] # change sign in ice

   v_flux = ncfile.createVariable('v_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
   v_flux.units = 'Pa'
   v_flux.missing_value = 1.e+20
   v_flux.long_name = 'j-direction wind stress'
   v_flux[:] = 0.0 # change sign in ice

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
   lt[:] = 0.0 # change sign in ice

   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
   
   return

def make_topo(x,y,args):
   # parameters
   name = 'ocean_topog'
   H = args.max_depth # max depth
   nx = len(x); ny = len(y)
   Hs = 500.0  # shelf's depth
   Ys = 500 * 1.0e3 # len of shelf in m
   Ws = 50 * 1.0e3 # len of slope in m
   [X,Y] = np.meshgrid(x,y)
   nx = len(x); ny = len(y)
   D = np.zeros((ny,nx))
   #np.ones((ny,nx)) *  H
   for j in range(ny):
      for i in range(nx):
          D[j,i] = Hs + 0.5 * (H-Hs) * (1.0 + np.tanh((Y[j,i]*1.0e3 - Ys)/Ws))

   # we need at least one land point,
   # so row 0 is set to land
   D[0,:] = 0.
   Dcoupler = D.copy()
   # 1) topography used in the coupler
   # open a new netCDF file for writing.
   name = 'ocean_topog'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('ntiles',1)
 
   # create variables
   nxx = ncfile.createVariable('nx',np.dtype('double').char,('nx'))
   nxx.units = 'km'
   nxx.description = 'x location of cell centers'
   nxx.long_name = 'x location of cell centers'
   nxx[:] = x[:]

   nyy = ncfile.createVariable('ny',np.dtype('double').char,('ny'))
   nyy.units = 'km'
   nyy.description = 'y location of cell centers'
   nyy.long_name = 'y location of cell centers'
   nyy[:] = y[:]

   depth = ncfile.createVariable('depth',np.dtype('float32').char,('ny','nx'))
   depth.units = 'm'
   depth.description = 'depth at h points'
   depth.long_name = 'depth at h points'
   depth[:,:] = Dcoupler[:,:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')

   # 2) topography used ny the ocean
   # open a new netCDF file for writing.
   name = 'topog'
   ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)
   ncfile.createDimension('ntiles',1)

   # create variables
   nxx = ncfile.createVariable('nx',np.dtype('double').char,('nx'))
   nxx.units = 'km'
   nxx.description = 'x location of cell centers'
   nxx.long_name = 'x location of cell centers'
   nxx[:] = x[:]

   nyy = ncfile.createVariable('ny',np.dtype('double').char,('ny'))
   nyy.units = 'km'
   nyy.description = 'y location of cell centers'
   nyy.long_name = 'y location of cell centers'
   nyy[:] = y[:]

   depth = ncfile.createVariable('depth',np.dtype('float32').char,('ny','nx'))
   depth.units = 'm'
   depth.description = 'depth at h points'
   depth.long_name = 'depth at h points'
   depth[:,:] = D[:,:]
   ncfile.close()
   print ('*** SUCCESS creating '+name+'.nc!')
 
   return Dcoupler

# A helper function to use when writing strings in a netcdf file
def set_string(variable, value):
   """Sets "variable" to "value" padded with blanks where
   "variable" is a netcdf variable object and "value" is a string."""
   variable[:] = '\000' * variable.shape[0]
   variable[:len(value)] = value
   return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
