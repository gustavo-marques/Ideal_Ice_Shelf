#!/usr/bin/env python

# Compute various diagnostics for the Idealized Ice Shelf problem.
# Gustavo Marques

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import wright_eos as eos
import gsw
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

  parser.add_argument('-n', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-prog_file', type=str, default='ocean_month_z.nc', help='''Name of the prog. file (default is prog.nc).''')

  parser.add_argument('-ice_file', type=str, default='ice_month.nc', help='''Name of the ice file (default is ice_month.nc).''')

  parser.add_argument('-sfc_file', type=str, default='ocean_sfc.nc', help='''Name of the ocean surface fluxes file (default is ocean_sfc.nc).''')

  parser.add_argument('-ice_shelf_file', type=str, default='MOM_Shelf_IC.nc', help='''Name of the file that has the initial conditions for the ice shelf (default is MOM_Shelf_IC.nc).''')

  parser.add_argument('-cshelf_lenght', type=float, default=470.,
     help='''Continental shelf lenght in the y direction (km). Default is 470.''')

  parser.add_argument('-t0', type=int, default=0,
     help='''Initial time indice to start computations. Default is 0.''')

  parser.add_argument('-tf', type=int, default=-1,
     help='''Final time indice for computations. Default is -1 (all starting from t0).''')

  parser.add_argument('-AABW_rho', type=float, default=1037.15,
     help='''Minimum density (sigma2) to compute AABW transport (kg/m3). Default is 1037.15''')

  parser.add_argument('-CDW_rho', type=float, default=1037.05,
     help='''Minimum density (sigma2) to compute CDW transport (kg/m3). Default is 1037.05''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   # load a few variables
   if args.tf == -1:
      time = Dataset(args.prog_file).variables['time'][args.t0::]/365. # in years
   else:
      time = Dataset(args.prog_file).variables['time'][args.t0:args.tf+1]/365.

   x = Dataset(args.prog_file).variables['xh'][:] # in km
   y = Dataset(args.prog_file).variables['yh'][:]
   dy = (y[1]-y[0])*1.0e3 # in m
   zt = Dataset(args.prog_file).variables['zt'][:]
   zw = Dataset(args.prog_file).variables['zw'][:]
   h = zw[1::]-zw[0:-1]
   km, jm, im = len(zt), len(x), len(y)
   depth = np.tile(zt,jm*im).reshape(jm,im,km).transpose()
   h = np.tile(h,jm*im).reshape(jm,im,km).transpose()
   pressure = gsw.p_from_z(-depth,-75.) * 1.0e4

   # ice shelf lenth
   name = args.n
   # create ncfile and zero fields.
   create_ncfile(name,x,y,time,args)

   # lists used to save ncdata
   varname = []; var = []

   # create arrays
   oht1 = np.zeros(len(time)); var.append('oht1'); varname.append('OHT_shelf')
   oht2 = np.zeros(len(time)); var.append('oht2'); varname.append('OHT_ice_shelf')

   # loop in time
   for t in range(args.t0,len(time)+args.t0):
           tt = t - args.t0 # time indice used in the arrays
           print 'Time (years):', time[tt]
	   # load data
	   v = mask_bad_values(Dataset(args.prog_file).variables['v'][t,:])
           vh = v * h * dy
	   temp = mask_bad_values(Dataset(args.prog_file).variables['temp'][t,:])
	   salt = mask_bad_values(Dataset(args.prog_file).variables['salt'][t,:])
           pressure = gsw.p_from_z(depth,-75.) * 1.0e4 # in Pa [1 db = 10e4 Pa]
           # diags functions
           oht1[tt] = get_oht(temp,salt,depth,vh,y,y_loc = 460.)
           print 'oht1',oht1[tt]
           oht2[tt] = get_oht(temp,salt,depth,vh,y,y_loc = 190.)
           print 'oht2',oht2[tt]
           print '\n'

   print 'Saving netcdf data...'

   for i in range(len(varname)):
       s = str("ncwrite(name,'%s',%s)"% (varname[i],var[i]))
       print s
       eval(s)
       #ncwrite(name,varname,var)

   print 'Done!'


def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset(name+'.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def mask_bad_values(data,value = -1.30856803e+26):
         return np.ma.masked_where(data == value, data)

def get_oht(t,s,d,vh,y,y_loc):
         '''
         Compute the onshore heat transport, as defined in St-Laurent et al JPO 2012
         '''
         cp = 3974.0 # heat capacity
         rho0 = 1028.0
         tmp = np.nonzero(y<=y_loc)[0][-1]
         # mask transport. > 0.
         vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])
         t_new = np.ma.masked_where(vh[:,tmp,:]>0.0, t[:,tmp,:])
         s_new = np.ma.masked_where(vh[:,tmp,:]>0.0, s[:,tmp,:])
         d_new = np.ma.masked_where(vh[:,tmp,:]>0.0, d[:,tmp,:])
         p = gsw.p_from_z(d_new,-70.0) * 1.0e4 # in Pascal
         t_freeze = eos.tfreeze(s_new,p)
         dt = (t_new - t_freeze)
         print 't-tf min/max',dt.min(),dt.max()

         oht = - (vhnew*cp*rho0*dt).sum() # watts
         return oht

def create_ncfile(exp_name, xx, yy, ocean_time, args): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'.nc','w',format='NETCDF4')
   # dimensions
   nx = len(xx) ; ny = len(yy)
   # create dimensions.
   ncfile.createDimension('time', len(ocean_time))
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)

   # create variables, assign units and provide decription
   x = ncfile.createVariable('x',np.dtype('float32').char,('nx'))
   x.units = 'km'
   x.description = 'x location of cell centers'
   x.long_name = 'x location of cell centers'
   y = ncfile.createVariable('y',np.dtype('float32').char,('ny'))
   y.units = 'km'
   y.description = 'y location of cell centers'
   y.long_name = 'y location of cell centers'

   time = ncfile.createVariable('time',np.dtype('float32').char,('time'))
   time.units = 'years'
   time.description = 'time since start of simulation'
   time.long_name = 'time since start of simulation'

   OHT_shelf = ncfile.createVariable('OHT_shelf',np.dtype('float32').char,('time'))
   OHT_shelf.units = 'Watts'; OHT_shelf.description = 'Onshore heat transport at shelf break'

   OHT_ice_shelf = ncfile.createVariable('OHT_ice_shelf',np.dtype('float32').char,('time'))
   OHT_ice_shelf.units = 'Watts'; OHT_ice_shelf.description = 'Onshore heat transport at the ice shelf edge'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
