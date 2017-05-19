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

  parser.add_argument('-name', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-n', type=int, default=1, help='''Number of records to be skipped. Detault is 1, where all records are read.''')

  parser.add_argument('-energy', type=str, default='', help='''Extract the total energy (KE + PE) per unit mass (m2 s-2) from the given file and plot it as a function of time.''')

  parser.add_argument('-total_tke', help='''Plot and save volume integreated TKE (m^5/s2).''', action="store_true")

  parser.add_argument('-cshelf_lenght', type=float, default=470.,
     help='''Continental shelf lenght in the y direction (km). Default is 470.''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   if args.energy:
      print("Computing total energy...")
      compute_energy(args)

   if args.total_tke:
      print("Computing turbulent kinetic energy...")
      compute_total_tke(args)

   print 'Done!'

def compute_total_tke(args):
    n = args.n
    time = Dataset('prog.nc').variables['time'][::n]/365.
    y = Dataset('prog.nc').variables['yh'][:]
    ubar = Dataset('ocean_month.nc').variables['u'][:].mean(axis=0)
    vbar = Dataset('ocean_month.nc').variables['v'][:].mean(axis=0)
    area = np.tile(Dataset('ocean_geometry.nc').variables['Ah'][:],(ubar.shape[0],1,1))
    area_ice_shelf = np.tile(Dataset('MOM_Shelf_IC.nc').variables['shelf_area'][0,:,:],(ubar.shape[0],1,1))
    area_cshelf = np.ma.masked_where(area_ice_shelf>0,area)
    tmp1 = np.nonzero(y<=args.cshelf_lenght)[0][-1]
    area_cshelf[:,tmp1::,:] = 0.0
    area_cshelf = np.ma.masked_where(area_cshelf==0.,area_cshelf)

    TKE = np.zeros(len(time))
    TKE_CS = np.zeros(len(time))
    TKE_IS = np.zeros(len(time))

    for t in range(len(time)):
        print 'Time is (days):', time[t]
        y = Dataset('prog.nc').variables['yh'][:]
        u = Dataset('prog.nc').variables['u'][t*n,:]
        v = Dataset('prog.nc').variables['v'][t*n,:]
        h = Dataset('prog.nc').variables['h'][t*n,:]
        u_prime = u - ubar; v_prime = v - vbar
        # u and v at h points
        uh = np.zeros(u.shape); vh = np.zeros(v.shape);
        utmp = 0.5 * (u_prime[:,:,0:-1] + u_prime[:,:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
        vtmp = 0.5 * (v_prime[:,0:-1,:] + v_prime[:,1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
        uh[:,:,1::] = utmp; uh[:,:,0] = 0.5*u_prime[:,:,0] #u_i=1 = 0.5*u_(i=3/2)
        vh[:,1::,:] = vtmp; vh[:,0,:] = 0.5*v_prime[:,0,:] #v_j=1 = 0.5*v_(j=3/2)
        vh = vh*vh; uh = uh*uh # m^2/s
        vel = vh + uh
        vel_tot = 0.5*(vel * h * area) # m^5/s
        TKE[t] = vel_tot.sum()
        TKE_CS[t] = (vel_tot*area_cshelf/area_cshelf).sum()
        TKE_IS[t] = (vel_tot*area_ice_shelf/area_ice_shelf).sum()
        print 'Total TKE is:',TKE[t] 
        print 'TKE on the cont. shelf is:',TKE_CS[t] 
        print 'TKE under ice shelf is:',TKE_IS[t] 
        print '\n'

    # save netcdf file
    # open a new netCDF file for writing.
    ncfile = Dataset(args.name+'_tke.nc','w',format='NETCDF4')
    # create dimensions.
    ncfile.createDimension('time', len(time))

    # create variables, assign units and provide decription
    time1 = ncfile.createVariable('time',np.dtype('float32').char,('time'))
    time1.units = 'years'
    time1.description = 'time since start of simulation'
    time1.long_name = 'time since start of simulation'
    time1[:] = time

    TKE1 = ncfile.createVariable('total_tke',np.dtype('float32').char,('time'))
    TKE1.units = 'm^5/s2'; TKE1.description = 'Volume integrated TKE'
    TKE1[:] = TKE[:]

    TKE2 = ncfile.createVariable('cshelf_tke',np.dtype('float32').char,('time'))
    TKE2.units = 'm^5/s2'; TKE2.description = 'TKE integrated over cont. shelf'
    TKE2[:] = TKE_CS[:]

    TKE3 = ncfile.createVariable('ishelf_tke',np.dtype('float32').char,('time'))
    TKE3.units = 'm^5/s2'; TKE3.description = 'TKE integrated over ice shelf'
    TKE3[:] = TKE_IS[:]

def compute_energy(args):
    file = args.energy
    os.system("awk '/MOM Day/{print $3}' " + file + " > tmp0" )
    os.system("awk '/MOM Day/{print $6}' " + file + " > tmp1" )
    os.system("awk -F',' '{print $1}' tmp1 > tmp2")
    time = np.loadtxt('tmp0')/365. # in yr.
    energy = np.loadtxt('tmp2')
    os.system('rm tmp?')
    plt.figure()
    plt.plot(time,energy,lw=2.5)
    plt.xlabel('Time [years]')
    plt.ylabel('Total energy [m2 s-2]')
    plt.grid()
    plt.savefig(args.name+'_total_energy.png')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
