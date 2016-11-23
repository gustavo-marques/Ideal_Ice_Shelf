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

  parser.add_argument('-cshelf_lenght', type=float, default=500.,
      help='''Continental shelf lenght in the y direction (km). Default is 500.''')

  parser.add_argument('-slope_lenght', type=float, default=50.,
      help='''Continental shelf slope lenght in the y direction (km). Default is 50.''')

  parser.add_argument('-u10_asf', type=float, default=-5.0,
      help='''Max. wind speed at the ASF (m/2). Default is -5.0''')

  parser.add_argument('-L', type=float, default=1200.,
      help='''Domain lenght in the y direction (km). Default is 1.2E3''')

  parser.add_argument('-max_depth', type=float, default=4.0e3,
      help='''Maximum ocean depth (m). Default is 4E3.''')

  parser.add_argument('-min_depth', type=float, default=500.0,
      help='''Minimum ocean depth (m). Default is 500.''')

  parser.add_argument('-Qmax', type=float, default=70.0,
      help='''Maximum heat flux loss (W/m^2). Default is 70.''')

  parser.add_argument('-land_width', type=float, default=100.0,
      help='''Widht of land region next to ice shelf (km)- this also controls the shape
      of the ice shelf. Default is 100.''')

  parser.add_argument('-coupled_run', help='''Generate all the files needed to run an ocean_SIS2 simulation.''', action="store_true")
  
  parser.add_argument('-debug', help='''Adds prints and plots to help debug.''', action="store_true")
  
  parser.add_argument('-ice_shelf', help='''Generate ice shelf ncfile.''', action="store_true")

  parser.add_argument('-trough', help='''Adds a trough cutting the continental shelf.''', action="store_true")
  
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

   if args.ice_shelf:
      # create ice shelf
      make_ice_shelf(x,y,args)
   else:
      args.ISL = 0.0

   # create topography
   Ocean_Depth = make_topo(x,y,args)

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
   os.system('cp mosaic.nc grid_spec.nc')

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
   # all values between -78 to - 60
   # use sw.dist to find distance
   # 176 lon (Ross Sea), July, southern wall conditions (-78)
   i =  176; t = 6; j = 11
   temp_south = get_profile(t,i,j,'PTEMP',depth,z,args,'TempSouth')
   salt_south = get_profile(t,i,j,'SALT',depth,z,args,'SaltSouth')
   # shelf break conditions (-72)
   j = 17
   temp_break = get_profile(t,i,j,'PTEMP',depth,z,args,'TempSlope')
   salt_break = get_profile(t,i,j,'SALT',depth,z,args,'SaltSlope')
   # northern wall conditions (-65)
   j = 21
   temp_north = get_profile(t,i,j,'PTEMP',depth,z,args,'TempNorth')
   salt_north = get_profile(t,i,j,'SALT',depth,z,args,'SaltNorth')

   # 3d fields
   temp3D = np.zeros((1,args.nz,len(y),len(x)))
   salt3D = np.zeros((1,args.nz,len(y),len(x)))

   if args.homogeneous_ts:
    dist = [0,args.max_depth]
    f1 = interp1d(dist,[0, 2])
    dummy_t = f1(z)
    f1 = interp1d(dist,[33.,35.])
    dummy_s = f1(z)
    sigma2 = eos.wright_eos(dummy_t,dummy_s,2.0e7)
    # compute rho0/alpha/beta
    alpha = eos.alpha_wright_eos(dummy_t[0],dummy_s[0],2.0e7)
    beta = eos.beta_wright_eos(dummy_t[0],dummy_s[0],2.0e7)
    Rho_T0_S0 = eos.wright_eos(0.,0.,2.0e7) + 0.# 0.017 is a correction factor
    # compute linear eos
    rho_lin = Rho_T0_S0 + alpha*dummy_t + beta*dummy_t
    print 'sigma2 - rho_lin',sigma2 - rho_lin
    print 'alpha,beta,Rho_T0_S0',alpha.mean(),beta.mean(),Rho_T0_S0
    plt.close('all')
    plt.figure()
    plt.plot(sigma2-1000.,-z,'k',rho_lin-1000,-z,'r')
    plt.show()

    for i in range(args.nx):
       for j in range(args.ny):
           temp3D[0,:,j,i] = dummy_t[:]
           salt3D[0,:,j,i] = dummy_s[:]

   else:
     # distace from southern wall to use in the interp.
     #dist = [0,args.cshelf_lenght+args.slope_lenght*0.5,args.L,args.L]
     #dummy = args.ISL + args.cshelf_lenght + args.slope_lenght*0.5 
     dist = [0, args.ISL, args.L-100.,args.L]
     # horizontal interp
     for k in range(args.nz):
         #f1 = interp1d(dist, [temp_south[k],temp_break[k],temp_north[k], temp_north[k]])
         f1 = interp1d(dist, [temp_south[k],temp_south[k],temp_north[k], temp_north[k]])
         dummy_t = f1(y)
         #f1 = interp1d(dist, [salt_south[k],salt_break[k],salt_north[k], salt_north[k]])
         f1 = interp1d(dist, [salt_south[k],salt_south[k], salt_north[k], salt_north[k]])
         dummy_s = f1(y)
         for i in range(args.nx):
            temp3D[0,k,:,i] = dummy_t[:]
            salt3D[0,k,:,i] = dummy_s[:]

     # set temp below 50 m to freezing point
     #if args.freeze_ic:
     #	temp3D, salt3D = set_freezing_temp(temp3D,salt3D,z,y,args)

     # compute sigma2
     sigma2 = eos.wright_eos(temp3D,salt3D,2.0e7)
     # print min/max sigma2
     print 'sigma2.min(),sigma2.max()',sigma2.min(),sigma2.max()

     # compute rho0/alpha/beta
     alpha = eos.alpha_wright_eos(temp3D,salt3D,2.0e7)
     beta = eos.beta_wright_eos(temp3D,salt3D,2.0e7)
     Rho_T0_S0 = eos.wright_eos(0.,0.,2.0e7) + 0.017 # 0.017 is a correction factor
   
     # compute linear eos
     rho_lin = Rho_T0_S0 + alpha.mean()*temp3D + beta.mean()*salt3D

   layers = Dataset('GOLD_IC.2010.11.15.nc').variables['Layer'][:] # used in the global run with sig2
   if args.debug:
      print 'alpha,beta,Rho_T0_S0',alpha.mean(),beta.mean(),Rho_T0_S0
      print 'sigma2 - rho_lin',sigma2[0,:,0,0] - rho_lin[0,:,0,0]

      plt.figure()
      plt.contourf(y,-z,rho_lin[0,:,:,0])
      plt.colorbar()
      plt.contour(y,-z,rho_lin[0,:,:,0]-1000.,layers-1000,colors='k',linewidths=2)
      plt.title('Density - 1000.')
      plt.savefig('PNG/rho_section.png')
  
 
      plt.figure()
      plt.plot(sigma2[0,:,0,0]-1000.,-z,'b',rho_lin[0,:,0,0]-1000.,-z,'bx')
      plt.plot(sigma2[0,:,args.ny/2.,0]-1000.,-z,'r',rho_lin[0,:,args.ny/2.,0]-1000.,-z,'rx')
      plt.plot(sigma2[0,:,args.ny-1,0]-1000.,-z,'k',rho_lin[0,:,args.ny-1,0]-1000.,-z,'kx')
      plt.title('Linear (dashed) vs. nonlinear EoS at south, center and north (b,r,k)')
      plt.savefig('PNG/rho_profiles.png')
 
      # plot t - tfreeze at surface
      tf = eos.tfreeze(salt3D[0,0,:,:],1.0e5) 
      plt.figure()
      plt.subplot(311)
      plt.contourf(x,y,temp3D[0,0,:,:],np.linspace(temp3D[0,0,:,:].min(),temp3D[0,0,:,:].max(),50))
      plt.colorbar()
      plt.title('SST')
      plt.subplot(312)
      plt.contourf(x,y,tf,np.linspace(tf.min(),tf.max(),50))
      plt.colorbar()
      plt.title('Freezing temp. at surface')
      plt.subplot(313)
      plt.contourf(x,y,temp3D[0,0,:,:] - tf)
      plt.colorbar()
      plt.title('T - TFreeze')
      plt.savefig('PNG/SST_and_Tfreeze.png')

      #(S, T) = np.meshgrid(np.linspace(0.,40.,51), np.linspace(-2.5, 3.0, 51))
      #r1 = eos.wright_eos(T,S,2.0e7)
      #r2 = Rho_T0_S0 + alpha.mean()*T + beta.mean()*S

      #plt.figure()
      #CS1 = plt.contour(S,T,r1-1000,layers-1000,colors='b')
      #plt.clabel(CS1, inline=1, fontsize=10)
      #CS2 = plt.contour(S,T,r2-1000,layers-1000,colors='r')
      #plt.clabel(CS2, inline=1, fontsize=10)
      #plt.xlabel('Salt'); plt.ylabel('Temp')
      #plt.grid()
      #plt.show()

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

   # 2) The file from which the coordinate densities are read
   #name = 'ts_ic_profile'
   #ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   #ncfile.createDimension('Layer',args.nz)

   # create variables
   #PTEMP = ncfile.createVariable('PTEMP',np.dtype('double').char,('Layer'))
   #PTEMP.units = 'Celcius'
   #PTEMP.long_name = 'Potential Temperature'
   #PTEMP[:] = temp_int[:]

   #SALT = ncfile.createVariable('SALT',np.dtype('double').char,('Layer'))
   #SALT.units = 'PSU'
   #SALT.long_name = 'Salinity'
   #SALT[:] = salt_int[:]

   #ncfile.close()
   #print ('*** SUCCESS creating '+name+'.nc!')

def make_forcing(x,y,args):
   # forcing parameters
   Ly = args.L # domain size km
   W = args.W # domain width km
   CSL = args.cshelf_lenght  # km
   Q0 = 20. # W/m^2
   Qmax = args.Qmax
   Yt = 600.0  # km
   Lasf = 300.0 # km
   Lacc = 450.0
   tau_acc = 0.1
   tau_asf = -0.05 # default is -0.075 # N/m^2
   katabatic_wind = 0.05 # def is 0.05 N/m2
   sponge = 100.0 # km
   ISL = args.ISL
   heat_flux = 0.0
   nx = len(x); ny = len(y); nt =1 # time
   tau_x = np.zeros((nt,ny,nx)) # ice forcing
   tau_y = np.zeros((nt,ny,nx))
   heat = np.zeros((nt,ny,nx))

   # wind and heat
   # x-dir
   tmp1 = ISL + 100.
   tmp1 = 300.
   for j in range(ny):
      if y[j] < tmp1:
        tau_x[0,j,:] = 0.0
      elif y[j] >= tmp1 and y[j] <= tmp1 + Lasf:
         tmp = np.pi*(y[j]-tmp1)/Lasf
         tau_x[0,j,:] = tau_asf * np.sin(tmp)**2 
      elif y[j] > tmp1 + Lasf and y[j] <= tmp1 + Lasf + Lacc:
         tmp = np.pi*(y[j]-tmp1-Lasf)/Lacc 
         tau_x[0,j,:] = tau_acc * np.sin(tmp)**2
      else:
         tau_x[0,j,:] = 0.0
        
   # y-dir (decays linearly away from the ice shelf front)
   # has a gaussian shape: tauy = tauy_max * np.exp(((x-W/2.)**2)/(2*W_v10))
   W_v10 = 50000. # km
   dummy = ISL + args.cshelf_lenght #+ args.slope_lenght
   for i in range(nx):
      for j in range(ny):
         if y[j] < ISL:
            tau_y[0,j,i] = katabatic_wind
         elif y[j] <= dummy:
            tmp =  np.exp((-(x[i]-W*0.5)**2)/(2*W_v10))
            tau_y[0,j,i] = ((-katabatic_wind/(dummy)) * (y[j]) + katabatic_wind) * tmp

   # read ocean's initial SST (in K)
   T0 = Dataset('ic_ts.nc').variables['PTEMP'][0,1,:,0] + 273.
   # T10 = T10 - 10 * np.exp(-(y[j]-ISL)/150.)
   # e-folding scale is 150 km
   efold = args.cshelf_lenght + args.slope_lenght
   for j in range(ny):
      if y[j]<=ISL:
         heat[0,j,:] = -Qmax
      elif y[j]<=ISL + efold:
         heat[0,j,:] = -(Qmax - Q0) * np.exp(-(y[j]-ISL)/100.) - Q0
      else:
         heat[0,j,:] = -Q0
   
   # evap, proxy for brine formation in polynyas
   #if t10_polynya > 0.:
   #  for i in range(nx):
   #    for j in range(ny):
   #	 if (x[i] - W/2. >= -major and x[i] - W/2. <= major):
   #	    if (y[j]-ISL >= 0.) and (y[j]-ISL <= (minor * np.abs((1-(x[i]-W/2.)**2/major**2)**(1/2.)))):
   #	      evaporation[0,j,i] = salt_flux
   #	      salt[0,j,i] = salt_flux
   #  	      heat[0,j,i] = heat_flux
   #	      latent[0,j,i] = latent_flux
   #	      t10[0,j,i] = t10_polynya

   # set salt flux such that net buoynacy flux is zero
   # B_heat = B_salt
   #cp = 3974.0 # J/(K kg)
   #rho0 = 1028.0; alpha = 0.11069/rho0; beta = 0.7906/rho0; g = 9.8
   #B_heat = heat[0,:,:] * g * alpha/ (rho0 * cp)
   #salt[0,:,:] = B_heat * rho0 / (beta * g)
   #print 'B_heat',B_heat.min(),B_heat.max()
   #print 'salt',salt.min(),salt.max()
   # plots
   if args.debug:

      #plt.figure()
      #u=u10[0,::5,::5]; v=v10[0,::5,::5]; mag = np.sqrt(u**2 + v**2)
      #plt.quiver(x[::5],y[::5],u, v, mag, cmap = plt.cm.seismic)
      #plt.colorbar()
      #plt.plot(x,np.ones(nx)*ISL,'k')
      #plt.plot(x,np.ones(nx)*(ISL+CSL),'k')
      #plt.title('Wind vel. at 10 m (m/s)')
      #plt.xlabel('x [km]'); plt.ylabel('y [km]')
      #plt.savefig('PNG/wind10m.png')

      #plt.figure()
      #plt.subplot(211)
      #plt.plot(y,u10[0,:,1])
      #plt.title('u10'); plt.xlabel('y [km]'); plt.ylabel('m/s')
      #plt.subplot(212)
      #plt.contourf(x,y,u10[0,:,:]); plt.colorbar()
      #plt.xlabel('x [km]'); plt.ylabel('y [km]')
      #plt.savefig('PNG/u10.png')

      #plt.figure()
      #plt.subplot(211)
      #plt.plot(y,v10[0,:,1])
      #plt.title('v10'); plt.xlabel('y [km]'); plt.ylabel('m/s')
      #plt.subplot(212)
      #plt.contourf(x,y,v10[0,:,:]); plt.colorbar()
      #plt.xlabel('x [km]'); plt.ylabel('y [km]')
      #plt.savefig('PNG/v10.png')

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
      plt.plot(y,tau_y[0,:,args.nx/2.])
      plt.title('tauy'); plt.xlabel('y [km]'); plt.ylabel('Pa')
      plt.subplot(212)
      plt.contourf(x,y,tau_y[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/tauy.png') 

      plt.figure()
      plt.subplot(211)
      plt.plot(y,heat[0,:,1])
      plt.title('Sensible Heat'); plt.xlabel('y [km]'); plt.ylabel('W/m^2')
      plt.subplot(212)
      plt.contourf(x,y,heat[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/sensible_heat.png')

      plt.figure()
      plt.title('Salt flux kg/(m^2 s)')
      plt.contourf(x,y,heat[0,:,:]); plt.colorbar()
      plt.xlabel('x [km]'); plt.ylabel('y [km]')
      plt.savefig('PNG/salt_flux.png')

   # create ncfile
   # open a new netCDF file for writing.
   # # used when forcing is applied in the atm
   #name = 'forcing_10'
   #ncfile = Dataset(name+'.nc','w')
   # create dimensions.
   #ncfile.createDimension('LON',nx)
   #ncfile.createDimension('LAT',ny)
   #ncfile.createDimension('TIME',None)
   # create variables
   #LON = ncfile.createVariable('LON',np.dtype('double').char,('LON'))
   #LON.units = 'km'
   #LON.long_name = 'h point nominal longitude'
   #LON.cartesian_axis = 'X'
   #LON[:] = x[:]

   #LAT = ncfile.createVariable('LAT',np.dtype('double').char,('LAT'))
   #LAT.units = 'km'
   #LAT.long_name = 'h point nominal latitude'
   #LAT.cartesian_axis = 'Y'
   #LAT[:] = y[:]

   #time = ncfile.createVariable('TIME',np.dtype('double').char,('TIME'))
   #time.long_name = 'time'
   #time.units = 'days since 0001-01-01 00:00:00'
   #time.cartesian_axis = 'T'
   #time.calendar_type = 'NOLEAP'
   #time.calendar = 'NOLEAP'
   #time.bounds = 'time_bounds'
   #time[0] = 0

   #t_10 = ncfile.createVariable('T_10',np.dtype('float32').char,('TIME','LAT','LON')) 
   #t_10.long_name = 'Air Temperature'
   #t_10.units = 'Kelvin'
   #t_10[:] = t10[:]

   #u_10 = ncfile.createVariable('U_10',np.dtype('float32').char,('TIME','LAT','LON'))
   #u_10.long_name = 'U wind'
   #u_10.units = 'm/s'
   #u_10[:] = u10[:]

   #v_10 = ncfile.createVariable('V_10',np.dtype('float32').char,('TIME','LAT','LON'))
   #v_10.long_name = 'U wind'
   #v_10.units = 'm/s'
   #v_10[:] = v10[:]

   #ncfile.close()
   #print ('*** SUCCESS creating '+name+'.nc!')

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
     lt[:] = 0.0 # change sign in ice

     salt_flux = ncfile.createVariable('salt_flux',np.dtype('float32').char,('time', 'yh', 'xh'), fill_value = 1.e+20)
     salt_flux.units = 'kg/(m^2 s)'
     salt_flux.missing_value = 1.e+20
     salt_flux.long_name = 'salt flux'
     salt_flux[:] = 0.0 # + adds salt from ocean

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
     evap[:] = 0.0

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
   Hs = args.min_depth  # shelf's depth
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

   # add trough(s)
   trough = np.zeros((ny,nx))
   xx = x - args.W*0.5
   for i in range(nx):
       if x[i]<= (args.W*0.5):
          trough[:,i] = args.min_depth + 0.5 * (700 - args.min_depth) * (1+np.tanh((xx[i]+30.)/10.))
       else:
          trough[:,i] = args.min_depth + 0.5 * (700 - args.min_depth) * (1-np.tanh((xx[i]-30.)/10.))

   if args.trough:
      D[trough>D]=trough[trough>D]

   # remove corners or "right angles" in the topography   
   for i in range(nx):
      D[:,i] = gaussian_filter(D[:,i],2)

 
   if args.ice_shelf:
     H1 = 500. # # depth increase under ice shelf
     ISL = args.ISL  # lenght of ice shelf cavity in m
     L1 = (ISL * 1.0e3) 
     minor = 2*ISL
     major = 200.0
     Wl = args.land_width # widht of land region next to ice shelf in km
     for j in range(ny):
	for i in range(nx):
            #if X[j,i]<= W*0.5:
            #  land = (ISL + 10) - ISL * (np.tanh((Y[j,i]*1.0e3 - L1)/(L1/10.)) - 1.0)
            #else:
            #land = (ISL + 10) - ISL * 
            if X[j,i]<= major*0.5 or X[j,i]>= (W-major*0.5): # add land mask
               if Y[j,i]<= ISL + 10:
                  D[j,i] = 0.0
            else:
	       if Y[j,i]<= ISL: # depth under ice shelf or land
	          D[j,i] = Hs - (H1/2.0) * (np.tanh((Y[j,i]*1.0e3 - L1)/(L1/10.)) - 1.0) 
	          #dummy1 = X[j,i] - W*0.5
                  #dummy2 = Y[j,i] - ISL
                  #if (dummy1 <= -major or dummy1 >= major):
                  #   D[j,i] = 0.0
                  #else: 
                  #   if (dummy2 < -(minor*np.sqrt((1-(abs(dummy1)**2)/major**2)))):
                  #      D[j,i] = 0.0 # land

   # to avoid sea ice formation under ice shelves,
   # two topography files need to be constructed.
   # The coupler topo is where the cavity is masked.

   if not args.ice_shelf:
     D[0,:] = 0.
   else:
     D[0,0] = 0.

   Dcoupler = D.copy()
   for j in range(ny):
     if y[j]<= args.ISL:
       Dcoupler[j,:] = 0.0

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
