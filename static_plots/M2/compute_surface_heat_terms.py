from matplotlib import pyplot as plt
from netCDF4 import Dataset, MFDataset
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
      Surface heat terms for Idealized Ice Shelf problem.
      ''',
  epilog='Written by Gustavo Marques, Apr. 2018.')

  parser.add_argument('-dx', type=str, default='1',
      help='''Horizontal grid resolution (default = 1).''')

  parser.add_argument('-exp', type=str, default='M2_exp0',
      help='''Experiment name (default = M2_exp0).''')

  parser.add_argument('-save', type=str, default='T',
                help='''Save data into a .nc file (default = T).''')

  parser.add_argument('-out', type=str, default='',
      help='''Name of output file (default = '').''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp,area,out):
    if out == '':
      s=Dataset(exp+'/ice_month.nc')
    else:
      s=MFDataset(exp+'/ice_month*.nc')

    print 'Reading file', exp
    tm = len(s.variables['time'][:])
    SW = np.zeros(tm)
    LW = np.zeros(tm)
    SH = np.zeros(tm)
    LH = np.zeros(tm)
    HS = np.zeros(tm) # snow thick
    SV = np.zeros(tm) # snow vol
    SF = np.zeros(tm) # snow fall
    HI = np.zeros(tm) # ice thick
    IV = np.zeros(tm) # ice vol
    FRA = np.zeros(tm)
    BHEAT = np.zeros(tm)
    BMELT = np.zeros(tm)
    SST = np.zeros(tm)
    SSS = np.zeros(tm)
    ALB = np.zeros(tm)
    SALTF = np.zeros(tm)
    time = np.zeros(tm)
    for t in range(tm):
      time[t] = s.variables['time'][t]/365.
      print 'Time (years):',time[t]
      hi_tmp=s.variables['HI'][t,:]
      sst_tmp=s.variables['SST'][t,:]
      sss_tmp=s.variables['SSS'][t,:]
      sw_tmp=s.variables['SW'][t,:]
      lw_tmp=s.variables['LW'][t,:]
      sh_tmp=s.variables['SH'][t,:]
      lh_tmp=s.variables['LH'][t,:]
      hs_tmp=s.variables['HS'][t,:]
      fra_tmp=s.variables['FRAZIL'][t,:]
      bh_tmp=s.variables['BHEAT'][t,:]
      bm_tmp=s.variables['BMELT'][t,:]
      saltf_tmp=s.variables['SALTF'][t,:]
      sf_tmp=s.variables['SNOWFL'][t,:]
      alb_tmp=s.variables['ALB'][t,:]

      HI[t] = (area*hi_tmp).sum()/area.sum()
      SSS[t] = (area*sss_tmp).sum()/area.sum()
      SST[t] = (area*sst_tmp).sum()/area.sum()
      BHEAT[t] = (area*bh_tmp).sum()
      BMELT[t] = (area*bm_tmp).sum()
      HS[t] = (area*hs_tmp).sum()//area.sum()
      SW[t] = (area*sw_tmp).sum()
      LW[t] = (area*lw_tmp).sum()
      LH[t] = (area*lh_tmp).sum()
      SH[t] = (area*sh_tmp).sum()
      FRA[t] = (area*fra_tmp).sum()
      ALB[t] = (area*alb_tmp).sum()/area.sum()
      IV[t] = (area*hi_tmp).sum()
      SV[t] = (area*hs_tmp).sum()
      SALTF[t] = (area*saltf_tmp).sum()

      print 'SSS, SST, HI, ALB, FRAZIL: ',SSS[t], SST[t],HI[t], ALB[t], FRA[t]
    s.close()
    return SSS,SST,HI,BHEAT,BMELT,IV,HS,SV,SH,LH,SW,LW,FRA,ALB,SALTF,time

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
     area = Dataset(path+'/out1/ocean_geometry.nc').variables['Ah'][:] # m^2
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     area = Dataset(path+'/ocean_geometry.nc').variables['Ah'][:] # m^2
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]

   # mask cavity
   area = np.ma.masked_where(HI.mask == True, area) # m^2

   # get data
   SSS,SST,HI,BHEAT,BMELT,IV,HS,SV,SH,LH,SW,LW,FRA,ALB,SALTF,time=get_data(path+'/'+args.out,area,args.out)
   #print 'SSS, SST, MLD, RD, HI mean/std',SSS.mean(), SSS.std(),SST.mean(), SST.std(), MLD.mean(), MLD.std(), RD.mean(), RD.std(),HI.mean(), HI.std()

   if args.save=='T':
     print 'Saving data...'
     os.system('mkdir ncfiles')
     exp_name = args.exp+'_dx'+args.dx
     create_ncfile(exp_name, time)
     ncwrite(exp_name,'SSS',SSS)
     ncwrite(exp_name,'SST',SST)
     ncwrite(exp_name,'BHEAT',BHEAT)
     ncwrite(exp_name,'BMELT',BMELT)
     ncwrite(exp_name,'HI',HI)
     ncwrite(exp_name,'IV',IV)
     ncwrite(exp_name,'HS',HS)
     ncwrite(exp_name,'SV',SV)
     ncwrite(exp_name,'SH',SH)
     ncwrite(exp_name,'LH',LH)
     ncwrite(exp_name,'SW',SW)
     ncwrite(exp_name,'LW',LW)
     ncwrite(exp_name,'FRA',FRA)
     ncwrite(exp_name,'ALB',ALB)
     ncwrite(exp_name,'SALTF',SALTF)
     ncwrite(exp_name,'TOT',(SW+LW+SH+LH-FRA))
      
   print 'Done!'

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset('ncfiles/'+name+'_surface_heat_terms.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def create_ncfile(exp_name, ocean_time): # may add exp_type
   """
   netcdf file with heat terms.
   """
   # open a new netCDF file for writing.
   ncfile = Dataset('ncfiles/'+exp_name+'_surface_heat_terms.nc','w',format='NETCDF4')
   # dimensions
   # create dimensions.
   ncfile.createDimension('time', len(ocean_time))

   # create variables, assign units and provide decription
   time = ncfile.createVariable('time',np.dtype('float32').char,('time'))
   time.units = 'years'
   time.description = 'years of simulation'
   time.long_name = 'years of simulation'

   SSS = ncfile.createVariable('SSS',np.dtype('float32').char,('time'))
   SSS.units = 'psu'; SSS.description = 'Domain-averaged surface salinity'

   SST = ncfile.createVariable('SST',np.dtype('float32').char,('time'))
   SST.units = 'C'; SST.description = 'Domain-averaged surface pot. temperature'

   HI = ncfile.createVariable('HI',np.dtype('float32').char,('time'))
   HI.units = 'm'; HI.description = 'Domain-average sea ice thickness'

   IV = ncfile.createVariable('IV',np.dtype('float32').char,('time'))
   IV.units = 'm^3'; IV.description = 'Sea ice volume'

   HS = ncfile.createVariable('HS',np.dtype('float32').char,('time'))
   HS.units = 'm'; HS.description = 'Domain-average snow thickness'

   SV = ncfile.createVariable('SV',np.dtype('float32').char,('time'))
   SV.units = 'm^3'; SV.description = 'Snow volume'
  
   SW = ncfile.createVariable('SW',np.dtype('float32').char,('time'))
   SW.units = 'W'; SW.description = 'Total shortwave'

   LW = ncfile.createVariable('LW',np.dtype('float32').char,('time'))
   LW.units = 'W'; LW.description = 'Total longtwave'

   SH = ncfile.createVariable('SH',np.dtype('float32').char,('time'))
   SH.units = 'W'; SH.description = 'Total sensible heat'

   LH = ncfile.createVariable('LH',np.dtype('float32').char,('time'))
   LH.units = 'W'; LH.description = 'Total latent heat'

   FRA = ncfile.createVariable('FRA',np.dtype('float32').char,('time'))
   FRA.units = 'W'; FRA.description = 'Total heat from frazil (ocean to ice)'

   TOT = ncfile.createVariable('TOT',np.dtype('float32').char,('time'))
   TOT.units = 'W'; TOT.description = 'Total heat (SW+LW+SH+LH-FRA)'

   ALB = ncfile.createVariable('ALB',np.dtype('float32').char,('time'))
   ALB.units = 'm'; ALB.description = 'Domain-averaged albedo'

   SALTF = ncfile.createVariable('SALTF',np.dtype('float32').char,('time'))
   SALTF.units = 'kg/s'; SALTF.description = 'Total salt flux'
  
   SF = ncfile.createVariable('SF',np.dtype('float32').char,('time'))
   SF.units = 'kg/s'; SF.description = 'Total snow flux'

   BHEAT = ncfile.createVariable('BHEAT',np.dtype('float32').char,('time'))
   BHEAT.units = 'W'; BHEAT.description = 'Total ocean to ice heat flux'

   BMELT = ncfile.createVariable('BMELT',np.dtype('float32').char,('time'))
   BMELT.units = 'W'; BMELT.description = 'Total bottom surface melting energy flux'
 
   # write data to coordinate vars.
   time[:] = ocean_time[:]# in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'_mean_surface_properties.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
