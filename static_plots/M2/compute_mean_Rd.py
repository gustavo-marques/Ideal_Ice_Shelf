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


def get_data(exp, ind1, ind2):
    s=Dataset(exp+'/prog.nc')
    print 'Reading file', exp
    dx = s.variables['xh'][1] - s.variables['xh'][0]
    print 'Dx is:', dx
    tmp = len(s.variables['time'][:])
    if tmp>144:
      rd=s.variables['Rd_dx'][-144::,:,:]*dx
    else:
      rd=s.variables['Rd_dx'][:,:,:]*dx
    
    s.close()
    m = np.zeros(rd.shape[0])
    m1 = np.zeros(rd.shape[0])
    m2 = np.zeros(rd.shape[0])
    for t in range(rd.shape[0]):
        tmp = np.ma.masked_where(rd[t,:] <0.01 , rd[t,:])
        m[t] = tmp.mean() # in km TOTAL
        m1[t] = tmp[0:ind2,:].mean() # in km C. Shelf and I shelf
        m2[t] = tmp[0:ind1,:].mean() # in km Ice Shelf
        
    return m, m1, m2

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
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]

  
   ISL = y[HI[:,0].mask == True,0][-1] 
   print('Ice shelf lenght is (km):',ISL)
   y_ishelf = np.nonzero(y[:,0]<=ISL)[0][-1]
   y_cshelf = np.nonzero(y[:,0]<=400.)[0][-1]
   print 'y_ishelf, y_cshelf', y_ishelf, y_cshelf

   rd_total, rd_ishelf, rd_cshelf=get_data(path+'/'+args.out,y_ishelf, y_cshelf)
   print 'rd_total, rd_ishelf, rd_cshelf:',rd_total, rd_ishelf, rd_cshelf
   print 'Saving Rds ...'
   os.system('mkdir TXT')
   #print ' mean/std (km):',rd.mean(), rd.std()
   text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_Rd.txt', "w")
   text_file.write("%f %f \n" % (rd_total.mean(), rd_total.std()))
   text_file.write("%f %f \n" % (rd_ishelf.mean(), rd_ishelf.std()))
   text_file.write("%f %f \n" % (rd_cshelf.mean(), rd_cshelf.std()))
   text_file.close()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
