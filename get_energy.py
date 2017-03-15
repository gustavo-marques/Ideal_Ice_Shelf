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

  parser.add_argument('-energy', type=str, default='', help='''Extract the total energy (KE + PE) per unit mass (m2 s-2) from the given file and plot it as a function of time.''')

  parser.add_argument('-total_tke', type=str, default='', help='''Plot the total domain turbulent kinetic energy per unit mass (m2 s-2) from the given file.''')

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

    # plot
    plt.figure()
    plt.plot(time/365.,TKEu,'r',lw=2.5)
    plt.plot(time/365.,TKEv,'b',lw=2.5)
    plt.xlabel('Time [years]')
    plt.ylabel('Total TKE (red = u, blue = v) [m2 s-2]')
    plt.grid()
    plt.savefig(args.n+'_total_tke.png')

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
    plt.savefig(args.n+'_total_energy.png')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
