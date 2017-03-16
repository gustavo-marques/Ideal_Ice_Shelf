#!/usr/bin/env python
import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import warnings
import os

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Print some basic statistics from a diag file.
      ''',
  epilog='Written by Gustavo Marques, Mar. 2017.')

  parser.add_argument('-var', type=str, default='', help='''Specifies which variable. If empty (default), stats from all variables will be printed.''')

  parser.add_argument('-fname', type=str, default='M2_exp2_2km.nc', help='''Specifies the name of the netCDF file.''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   if args.var:
      print 'reading variable', args.var + '...'
      data = Dataset(args.fname).variables[args.var][:]
      print 'mean/std:',data.mean(), data.std()

   else:
      print 'I will process all variables ...'
      process_all(args)

   print 'Done!'

def process_all(args):
    file = Dataset(args.fname)
    for f in file.variables:
       print 'reading variable', f + '...'
       data = file.variables[f][:]
       print 'mean/std:',data.mean(), data.std()
       print '######## \n'
    return

if __name__ == '__main__': parseCommandLine()
