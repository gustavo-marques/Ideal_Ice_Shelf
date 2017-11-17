from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import os
import argparse

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

  parser.add_argument('-dx', type=str, default='1',
      help='''Horizontal grid resolution (default = 1).''')

  parser.add_argument('-exp', type=str, default='M2_exp0',
      help='''Experiment name (default = M2_exp0).''')

  parser.add_argument('-out', type=str, default='out6',
      help='''Name of output file (default = out6).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp):
    s=Dataset(exp+'/ocean_month.nc')
    print 'Reading file', exp
    tmp = len(s.variables['time'][:])
    if tmp>24:
      melt=s.variables['mass_flux'][-24::,:,:]
    else:
      melt=s.variables['mass_flux'][:,:,:]
    
    s.close()
    m = np.zeros(melt.shape[0])
    for t in range(melt.shape[0]):
        m[t] = melt[t,:,:].sum() * 1.0e-12 * 3600 * 24 * 365 # in Gt/yr 
    return m

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   if os.path.exists(path+'/out1/'):
     x = Dataset(path+'/out1/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/out1/ocean_geometry.nc').variables['geolat'][:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]

   melt=get_data(path+'/'+args.out)
   print 'melt',melt
   print 'Saving melt rates...'
   os.system('mkdir TXT')
   print 'Melt rate is mean/std (Gt/yr):',melt.mean(), melt.std()
   text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_melt_rate.txt', "w")
   text_file.write("%f \n" % melt.mean())
   text_file.write("%f \n" % melt.std())
   text_file.close()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
