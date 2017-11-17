from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
from matplotlib.colors import Normalize
matplotlib.rcParams.update({'font.size': 16})
import sys, os
sys.path.append('../../')
from remapping import mom_remapping
from misc import *
import gsw
import wright_eos as eos
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


def get_data(exp,area,y_loc):
    s=Dataset(exp+'/prog.nc')
    print 'Reading file', exp
    y = s.variables['yh'][:]
    dy = ((y[1]-y[0])*1.0e3)
    area = (dy)**2 # in m^2
    y_tmp = np.nonzero(y<=y_loc)[0][-1]
    tm = len(s.variables['time'][:])
    tmp = 60
    if tm<tmp:
        tmp = tm
    tmp = tm
    Hout = np.zeros(tmp)
    Hin = np.zeros(tmp)
    Heat = np.zeros(tmp)
    transp = np.zeros(tmp)
    time = np.zeros(tmp)
    for t in range(1,tmp):
      time0 = s.variables['time'][t-1] * 3600 * 24. # in sec
      time1 = s.variables['time'][t] * 3600 * 24. # in sec
      dt = time1 - time0
      print 'Time is (secs):',time1, dt
      # at h points
      vh1 = 0.5*(s.variables['vh'][t,:,y_tmp-1,:] + s.variables['vh'][t,:,y_tmp,:])
      temp0 = s.variables['temp'][t-1,:] + 273.
      temp1 = s.variables['temp'][t,:] + 273.
      h0 = s.variables['h'][t-1,:,0:y_tmp,:]
      h1 = s.variables['h'][t,:,0:y_tmp,:]
      transp[t] = ((temp1[:,y_tmp,:])*vh1).sum()
      dTdt = ((((temp1[:,0:y_tmp,:]*h1) - (temp0[:,0:y_tmp,:]*h0))/dt)*area).sum()
      print 'dTdt, transp[t]',dTdt, transp[t]

    s.close()
    return #Hout/1.0e12, Hin/1.0e12, transp/1.0e12, time # in TW

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   #plt.style.use('ggplot')
   color1 = '#8B4513'
   color2 = '#ff6347'
   color3 = '#8470ff'
   color4 = '#3cb371'

   # plot some metrics for runs with varing wind forcing
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   if os.path.exists(path+'/out1/'):
     x = Dataset(path+'/out1/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/out1/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/out1/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/out1/MOM_Shelf_IC.nc').variables['shelf_area'][:]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/MOM_Shelf_IC.nc').variables['shelf_area'][:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]

   e0 = Dataset(path+'/'+args.out+'/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)
   ISL = y[HI[:,0].mask == True,0][-1] - 10.0 # -10 km to make sure it is under cavity
   #ISL = 150.0
   print('Ice shelf lenght is (km):',ISL)

   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   #Hout,Hin, transp, time=get_data(path+'/'+args.out,area,ISL)
   get_data(path+'/'+args.out,area,ISL)
   print 'Hout',Hout.mean(),Hout.std()
   print 'Hin',Hin.mean(),Hin.std()
   plt.figure()
   plt.plot(time,transp)
   plt.show()
   print 'Saving data...'
   os.system('mkdir TXT')
   #text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_heat_budget.txt', "w")
   #text_file.write("%f \n" % (Hout.mean()))
   #text_file.write("%f \n" % (Hout.std()))
   #text_file.write("%f \n" % (Hin.mean()))
   #text_file.write("%f \n" % (Hin.std()))
   #text_file.close()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
