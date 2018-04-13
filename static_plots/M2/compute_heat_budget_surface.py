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
    s=Dataset(exp+'/ocean_sfc.nc')
    s1=Dataset(exp+'/ice_month.nc')
    print 'Reading file', exp
    y = s.variables['yh'][:]
    tm = len(s.variables['time'][:])
    tmp = 144
    if tm<tmp:
        tmp = tm
    tmp = tm
    Hice = np.zeros(tmp-1)
    Hoff = np.zeros(tmp-1)
    BM = np.zeros(tmp-1)
    BH = np.zeros(tmp-1)
    time = np.zeros(tmp-1)
    for t in range(1,tmp):
      current_time = s.variables['time'][t-tmp]
      time[t-1] = current_time/365
      print 'Time is (years):',time[t-1]
      heat_net=s.variables['net_heat_surface'][t-tmp,:]
      bheat=s1.variables['BHEAT'][t-tmp,:]
      bmelt=s1.variables['BMELT'][t-tmp,:]
      Hice[t-1] = (heat_net*area).sum()
      tmp2 = np.ma.masked_where(area>0,heat_net)
      area_cell = area.max()
      Hoff[t-1] = (area_cell*tmp2).sum()
      BM[t-1] = (area_cell*bmelt).sum()
      BH[t-1] = (area_cell*bheat).sum()
      print 'BM[t-1], BH[t-1]',BM[t-1]/1.0e12, BH[t-1]/1.0e12
    s.close()
    s1.close()
    return Hice/1.0e12, Hoff/1.0e12, BM/1.0e12, BH/1.0e12, time #in TW

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
     area = Dataset(path+'/out1/MOM_Shelf_IC.nc').variables['shelf_area'][0,:]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/MOM_Shelf_IC.nc').variables['shelf_area'][0,:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]


   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   Hice,Hoff,BM,BH,time=get_data(path+'/'+args.out,area)
   print 'net heat from coupler ice_shelf/offshore mean/std',Hice.mean(), Hice.std(),Hoff.mean(), Hoff.std()
   print 'net BMELT/BHEAT mean/std',BM.mean(), BM.std(),BH.mean(), BH.std()
   if args.save=='T':
     print 'Saving data...'
     os.system('mkdir ncfiles')
     exp_name = args.exp+'_dx'+args.dx
     create_ncfile(exp_name, x[0,:], y[:,0], time)
     ncwrite(exp_name,'BMELT',BM)
     ncwrite(exp_name,'BHEAT',BH)
     ncwrite(exp_name,'BCOND',BH-BM)
     ncwrite(exp_name,'HNET',Hice+Hoff)
   print 'Done!'

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset('ncfiles/'+name+'_heat_budget_surface.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def create_ncfile(exp_name, xx, yy, ocean_time): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset('ncfiles/'+exp_name+'_heat_budget_surface.nc','w',format='NETCDF4')
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

   BMELT = ncfile.createVariable('BMELT',np.dtype('float32').char,('time'))
   BMELT.units = 'Watts'; BMELT.description = 'Total bottom surface melting energy flux'

   BHEAT = ncfile.createVariable('BHEAT',np.dtype('float32').char,('time'))
   BHEAT.units = 'Watts'; BHEAT.description = 'Total ocean to ice heat flux'

   BCOND = ncfile.createVariable('BCOND',np.dtype('float32').char,('time'))
   BCOND.units = 'Watts'; BCOND.description = 'Total ice conductive heat flux'

   HNET = ncfile.createVariable('HNET',np.dtype('float32').char,('time'))
   HNET.units = 'Watts'; HNET.description = 'Total surface ocean heat flux'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]-ocean_time[0]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'_heat_budget_surface.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
