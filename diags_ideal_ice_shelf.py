#!/usr/bin/env python

# Compute various diagnostics for the Idealized Ice Shelf problem.
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

  parser.add_argument('-exp_name', type=str, default='test', help='''Name of the experiment (default is test).''')
  
  parser.add_argument('-prog_file', type=str, default='prog.nc', help='''Name of the prog. file (default is prog.nc).''')
  
  parser.add_argument('-ice_file', type=str, default='ice_month.nc', help='''Name of the ice file (default is ice_month.nc).''')

  parser.add_argument('-energy', type=str, default='', help='''Extract the total energy (KE + PE) per unit mass (m2 s-2) from the given file and plot it as a function of time.''')
  
  parser.add_argument('-total_tke', type=str, default='', help='''Plot the total domain turbulent kinetic energy per unit mass (m2 s-2) from the given file.''')

  parser.add_argument('-cshelf_lenght', type=float, default=560.,
     help='''Continental shelf lenght in the y direction (km). Default is 560.''')

  parser.add_argument('-savefig', help='''Save the plots as PNG.''', action="store_true")

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   # load a few variables
   time = Dataset(args.ice_file).variables['time'][:]/365. # in years
   x = Dataset(args.prog_file).variables['xh'][:] # in km
   y = Dataset(args.prog_file).variables['yh'][:]
   HI = Dataset('ice_month.nc').variables['HI'][0,:]
   # ice shelf lenth
   args.ISL = y[HI[:,0].mask == True][-1]

   # create arrays
   salt_flux = np.zeros(len(time))
   polynya_area = np.zeros(len(time))
   ice_area = np.zeros(len(time))
   HI_max = np.zeros(len(time))

   # loop in time
   for t in range(len(time)):
           print 'Time (years):', time[t]
	   # read data
	   saltf = mask_bad_values(Dataset(args.ice_file).variables['SALTF'][t,:])
	   CI_tot = mask_bad_values(Dataset(args.ice_file).variables['CI_tot'][t,:])
	   HI = mask_bad_values(Dataset(args.ice_file).variables['HI'][t,:])
           # diags
	   salt_flux[t] = get_saltf(x,y,saltf,CI_tot,args)
	   polynya_area[t] = get_polynya_area(x,y,CI_tot,args)
	   ice_area[t] = get_ice_area(x,y,CI_tot)
	   HI_max[t] = HI.max()
   
   if args.energy:
      print("Plotting total energy...")
      plot_energy(args)

   if args.total_tke:
      print("Plotting turbulent kinetic energy...")
      plot_total_tke(args)

   if args.savefig:
	   print 'Plotting...'
           # salt flux
           plt.figure()
           plt.plot(time,salt_flux,'r',lw=2.5)
           plt.xlabel('Time [years]')
           plt.ylabel('Total (polynya) salt flux [kg / (m2 s)]')
           plt.grid()
           plt.savefig(args.exp_name+'_total_salt_flux.png')
 
           # polynya area
           plt.figure()
           plt.plot(time,polynya_area,'r',lw=2.5)
           plt.xlabel('Time [years]')
           plt.ylabel('Total polynya area [km^2]')
           plt.grid()
           plt.savefig(args.exp_name+'_total_polynya_area.png')

           plt.figure()
           plt.plot(time,ice_area,'r',lw=2.5)
           plt.xlabel('Time [years]')
           plt.ylabel('Total sea ice area [km^2]')
           plt.grid()
           plt.savefig(args.exp_name+'_total_ice_area.png')

 
           plt.figure()
           plt.plot(time,HI_max,'r',lw=2.5)
           plt.xlabel('Time [years]')
           plt.ylabel('Max. sea ice thickness [m]')
           plt.grid()
           plt.savefig(args.exp_name+'_max_hice.png')

   print 'Done!'

def mask_bad_values(data,value = -1.30856803e+26):
         return np.ma.masked_where(data == value, data)

def get_saltf(x,y,saltf,CI_tot,args):
	'''
        Compute the total salt flux into (+) the ocean in the cont. shelf 
	region (cshelf_lenght). Salt flux is just computed when CI_tot < 0.7. 
	'''
	# get indices for region of interest
        tmp = np.nonzero(y<=args.cshelf_lenght)[0]
	# mask points where CI_tot >= 0.7
	CI_tot = np.ma.masked_where(CI_tot >= 0.7, CI_tot)
	saltf = saltf * (CI_tot/CI_tot)
	return saltf[tmp,:].sum()

def get_polynya_area(x,y,CI_tot,args):
	'''
        Compute the total (coastal, i.e., for y <= ISL + 100 km) polynya area.
	'''
	grid_area = np.ones((CI_tot.shape)) * (x[1]-x[0]) * (y[1]-y[0])
        # get indices for region of interest
        tmp = np.nonzero(y<=args.ISL + 100.0)[0]
        # mask points where CI_tot >= 0.7
        CI_tot = np.ma.masked_where(CI_tot >= 0.7, CI_tot)
        grid_area = grid_area * (CI_tot/CI_tot)
	return grid_area[tmp,:].sum()

def get_ice_area(x,y,CI_tot):
	'''
        Compute total area covered by sea ice. Sea ice is present when CI_tot >= 0.7.
	'''
        grid_area = np.ones((CI_tot.shape)) * (x[1]-x[0]) * (y[1]-y[0])
	CI_tot = np.ma.masked_where(CI_tot < 0.7, CI_tot)
        grid_area = grid_area * (CI_tot/CI_tot)
	return grid_area.sum()

def plot_total_tke(args):
    file = args.total_tke
    time = Dataset(file).variables['time'][:]
    TKEu = np.zeros(len(time))
    TKEv = np.zeros(len(time))
    ubar = Dataset(file).variables['u'][0] * 0.; ubar2=0*ubar;
    vbar = Dataset(file).variables['v'][0] * 0.; vbar2=0*vbar;
    # using the following identity
    # <u'^2> = <(u - ubar)^2>, where <> = domain ave; u' = u-ubar
    #        = <u^2> - <ubar^2>
    for t in range(len(time)):
        print 'Time is (days):', time[t]
        u = Dataset(file).variables['u'][t,:]
        v = Dataset(file).variables['v'][t,:]
        ubar = ubar + u
        ubar2 = ubar2 + u*u
        vbar = vbar + v
        vbar2 = vbar2 + v*v
   
    ubar, vbar = ubar/len(time), vbar/len(time)
    print 'Finish computing mean, now tke...'
    for t in range(len(time)):
        print 'Time is (days):', time[t]
        u = Dataset(file).variables['u'][t,:]
        v = Dataset(file).variables['v'][t,:]
        up = (u - ubar)**2; vp = (v - vbar)**2
        TKEu[t] = up.sum() 
        TKEv[t] = vp.sum() 

    if args.savefig:
       # plot
       plt.figure()
       plt.plot(time/365.,TKEu,'r',lw=2.5)
       plt.plot(time/365.,TKEv,'b',lw=2.5)
       plt.xlabel('Time [years]')
       plt.ylabel('Total TKE (red = u, blue = v) [m2 s-2]')
       plt.grid()
       plt.savefig(args.exp_name+'_total_tke.png')

def plot_energy(args):
    file = args.energy
    os.system("awk '/MOM Day/{print $3}' " + file + " > tmp0" ) 
    os.system("awk '/MOM Day/{print $6}' " + file + " > tmp1" ) 
    os.system("awk -F',' '{print $1}' tmp1 > tmp2")
    time = np.loadtxt('tmp0')/365. # in yr.
    energy = np.loadtxt('tmp2')
    os.system('rm tmp?')
    if args.savefig:
       plt.figure()
       plt.plot(time,energy,lw=2.5)
       plt.xlabel('Time [years]')
       plt.ylabel('Total energy [m2 s-2]')
       plt.grid()
       plt.show()

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
