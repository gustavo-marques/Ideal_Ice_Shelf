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

  parser.add_argument('-save', type=str, default='T',
                help='''Save data into a .txt file (default = T).''')

  parser.add_argument('-out', type=str, default='out6',
      help='''Name of output file (default = out6).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp,area):
    s=Dataset(exp+'/prog.nc')
    s1=Dataset(exp+'/ice_month.nc')
    print 'Reading file', exp
    y = s.variables['yh'][:]
    tm = len(s.variables['time'][:])
    tmp = 144
    if tm<tmp:
        tmp = tm
    SSS = np.zeros(tmp-1)
    SST = np.zeros(tmp-1)
    MLD = np.zeros(tmp-1)
    RD = np.zeros(tmp-1)
    HI = np.zeros(tmp-1)
    time = np.zeros(tmp-1)
    for t in range(1,tmp):
      current_time = s.variables['time'][t-tmp]
      time[t-1] = current_time/365
      print 'Time is (years):',time[t-1]
      hi_tmp=s1.variables['HI'][t-tmp,:]
      sst_tmp=s1.variables['SST'][t-tmp,:]
      sss_tmp=s1.variables['SSS'][t-tmp,:]
      mld_tmp=s.variables['ePBL_h_ML'][t-tmp,:]
      rd_tmp=s.variables['Rd_dx'][t-tmp,:]
      tmp2 = np.ma.masked_where(sst_tmp.mask == True,rd_tmp)
      RD[t-1] = (area*tmp2).sum()/area.sum()
      tmp2 = np.ma.masked_where(sst_tmp.mask == True,mld_tmp)
      MLD[t-1] = (area*tmp2).sum()/area.sum()
      HI[t-1] = (area*hi_tmp).sum()/area.sum()
      SSS[t-1] = (area*sss_tmp).sum()/area.sum()
      SST[t-1] = (area*sst_tmp).sum()/area.sum()
      print 'SSS, SST, MLD, RD_DX, HI ',SSS[t-1], SST[t-1], MLD[t-1], RD[t-1], HI[t-1] 
    s.close()
    s1.close()
    return SSS,SST,MLD,RD,HI,time 

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
     area = Dataset(path+'/out1/ocean_geometry.nc').variables['Ah'][:]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/ocean_geometry.nc').variables['Ah'][:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]

   # mask cavity
   area = np.ma.masked_where(HI.mask == True, area)

   # get data
   SSS,SST,MLD,RD,HI,time=get_data(path+'/'+args.out,area)
   print 'SSS, SST, MLD, RD, HI mean/std',SSS.mean(), SSS.std(),SST.mean(), SST.std(), MLD.mean(), MLD.std(), RD.mean(), RD.std(),HI.mean(), HI.std()
   if args.save=='T':
     print 'Saving data...'
     os.system('mkdir ncfiles')
     exp_name = args.exp+'_dx'+args.dx
     create_ncfile(exp_name, x[0,:], y[:,0], time)
     ncwrite(exp_name,'SSS',SSS)
     ncwrite(exp_name,'SST',SST)
     ncwrite(exp_name,'MLD',MLD)
     ncwrite(exp_name,'RD_DX',RD)
     ncwrite(exp_name,'HI',HI)
   print 'Done!'

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset('ncfiles/'+name+'_mean_surface_properties.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def create_ncfile(exp_name, xx, yy, ocean_time): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset('ncfiles/'+exp_name+'_mean_surface_properties.nc','w',format='NETCDF4')
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
   time.description = 'last two years of simulation'
   time.long_name = 'last two years of simulation'

   SSS = ncfile.createVariable('SSS',np.dtype('float32').char,('time'))
   SSS.units = 'psu'; SSS.description = 'Domain-averaged surface salinity'

   SST = ncfile.createVariable('SST',np.dtype('float32').char,('time'))
   SST.units = 'C'; SST.description = 'Domain-averaged surface pot. temperature'

   MLD = ncfile.createVariable('MLD',np.dtype('float32').char,('time'))
   MLD.units = 'm'; MLD.description = 'Domain-average surface boundary layer depth'

   RD_DX = ncfile.createVariable('RD_DX',np.dtype('float32').char,('time'))
   RD_DX.units = 'nondim'; RD_DX.description = 'Domain-average RD_DX'

   HI = ncfile.createVariable('HI',np.dtype('float32').char,('time'))
   HI.units = 'm'; HI.description = 'Domain-average sea ice thickness'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]-ocean_time[0]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'_mean_surface_properties.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
