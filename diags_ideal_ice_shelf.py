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

  parser.add_argument('-energy', type=str, default='', help='''Extract the total energy (KE + PE) per unit mass (m2 s-2) from the given file and plot it as a function of time.''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   if args.energy:
      print("Plotting total energy...")
      plot_energy(args)

def plot_energy(args):
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
    plt.show()

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
