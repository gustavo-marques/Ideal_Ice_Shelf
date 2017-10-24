#!/usr/bin/env python


##!/usr/bin/env /ncrc/home2/Gustavo.Marques/anaconda2/bin/python

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
      Computes heat transports for Idealized Ice Shelf problem.
      ''',
  epilog='Written by Gustavo Marques, Oct. 2017.')

  parser.add_argument('-n', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-prog_file', type=str, default='prog.nc', help='''Name of the prog. file (default is prog.nc).''')

  parser.add_argument('-ice_file', type=str, default='ice_month.nc', help='''Name of the ice file (default is ice_month.nc).''')

  parser.add_argument('-tau_min', type=float, default=0.05,
                help='''Minimum tracer concentration to compute outflow volume flux (m^3/s). Default is 0.05''')

  parser.add_argument('-t0', type=int, default=0,
     help='''Initial time indice to start computations. Default is 0.''')

  parser.add_argument('-tf', type=int, default=-1,
     help='''Final time indice for computations. Default is -1 (all starting from t0).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   # load a few variables
   if args.tf == -1:
      time = Dataset(args.prog_file).variables['time'][args.t0::]/365. # in years
   else:
      time = Dataset(args.prog_file).variables['time'][args.t0:args.tf+1]/365.

   x = Dataset(args.prog_file).variables['xh'][:] # in km
   y = Dataset(args.prog_file).variables['yh'][:]
   HI = Dataset(args.ice_file).variables['HI'][0,:]
   # ice shelf lenth
   args.ISL = y[HI[:,0].mask == True][-1] - 10.0 # -10 km to make sure it is under cavity
   args.ISL_j = np.nonzero(y == args.ISL)
   print('Ice shelf lenght is ~ (km):',args.ISL)
   name = args.n
   # create ncfile and zero fields.
   create_ncfile(name,x,y,time,args)

   # lists used to save ncdata
   varname = []; var = []

   # create arrays
   ht1 = np.zeros(len(time)); var.append('ht1'); varname.append('OnHT')
   ht2 = np.zeros(len(time)); var.append('ht2'); varname.append('OfHT')
   E = np.zeros(len(time)); var.append('E'); varname.append('E')
   Qe = np.zeros(len(time)); var.append('Qe'); varname.append('Qe')
   Qin = np.zeros(len(time)); var.append('Qin'); varname.append('Qin')
   Qout = np.zeros(len(time)); var.append('Qout'); varname.append('Qout')
   # loop in time
   for t in range(args.t0,len(time)+args.t0):
           tt = t - args.t0 # time indice used in the arrays
           print 'Time (years):', time[tt]
	   # load data
	   vh = mask_bad_values(Dataset(args.prog_file).variables['vh'][t,:])
	   temp = mask_bad_values(Dataset(args.prog_file).variables['temp'][t,:])
	   salt = mask_bad_values(Dataset(args.prog_file).variables['salt'][t,:])
	   e = Dataset(args.prog_file).variables['e'][t,:]
           depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
           pressure = gsw.p_from_z(depth,-75.) * 1.0e4 # in Pa [1 db = 10e4 Pa]
           rho = eos.wright_eos(temp,salt,pressure)
           rhopot2 = eos.wright_eos(temp,salt,2.0e7)
           # diags functions
           ht1[tt] = get_ht(temp,salt,depth,vh,y,direction=1,y_loc = args.ISL)
           ht2[tt] = get_ht(temp,salt,depth,vh,y,direction=2,y_loc = args.ISL)
           print 'Heat in, Heat out:', ht1[tt], ht2[tt]
           # entrainment
           tr2 = mask_bad_values(Dataset(args.prog_file).variables['tr2'][t,:])
           # convert FW mass flux into volume flux kg/s to m^3/s
           vol_flux = mask_bad_values(Dataset(args.prog_file).variables['mass_flux'][t,:]) * 1.0e-3
           Qin[tt] = np.ma.masked_where(vol_flux<0.0, vol_flux).sum()
           print 'Total vol. flux of tracer across ice/ocean bnd (m^3/s):',Qin[tt]
           vhnew = np.ma.masked_where(vh[:,args.ISL_j,:]<=0.0, vh[:,args.ISL_j,:])
           vhnew = np.ma.masked_where(tr2[:,args.ISL_j,:]<args.tau_min, vhnew)
           # mask tr2 < tau_min
           tr2new = np.ma.masked_where(tr2[:,args.ISL_j,:]<args.tau_min, tr2[:,args.ISL_j,:])
           Qout[tt] = (vhnew * tr2new).sum() # in m^3/s
           # entrainment vol flux (m^3/s)
           Qe[tt] = Qout[tt] - Qin[tt]
           # Entrainment coeff
           E[tt] = (Qe[tt])/Qin[tt]
           print 'Total vol. flux of tracer across ice shelf front (m^3/s):',Qout[tt]
           print 'Qe, E',Qe[tt], E[tt]
           print '\n'

   print 'Saving netcdf data...'
   print 'Heat efficiency is:', (ht2.mean()-ht1.mean())/ht1.mean()
   for i in range(len(varname)):
       s = str("ncwrite(name,'%s',%s)"% (varname[i],var[i]))
       print s
       eval(s)
       #ncwrite(name,varname,var)

   print 'Done!'

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset(name+'_heat_transports.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def mask_bad_values(data,value = -1.30856803e+26):
         return np.ma.masked_where(data == value, data)

def get_ht(t,s,d,vh,y,direction,y_loc):
         '''
         Compute heat transport, onshore (direction=1) or offshore (direction=2), as defined in St-Laurent et al JPO 2012
         '''
         cp = 3974.0 # heat capacity
         rho0 = 1028.0
         tmp = np.nonzero(y<=y_loc)[0][-1]
         if direction ==1:
           # mask transport. > 0.
           vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])
           t_new = np.ma.masked_where(vh[:,tmp,:]>0.0, t[:,tmp,:])
           s_new = np.ma.masked_where(vh[:,tmp,:]>0.0, s[:,tmp,:])
           d_new = np.ma.masked_where(vh[:,tmp,:]>0.0, d[:,tmp,:])
         else:
           # mask transport. < 0.
           vhnew = np.ma.masked_where(vh[:,tmp,:]<0.0, vh[:,tmp,:])
           t_new = np.ma.masked_where(vh[:,tmp,:]<0.0, t[:,tmp,:])
           s_new = np.ma.masked_where(vh[:,tmp,:]<0.0, s[:,tmp,:])
           d_new = np.ma.masked_where(vh[:,tmp,:]<0.0, d[:,tmp,:])

         p = gsw.p_from_z(d_new,-70.0) * 1.0e4 # in Pascal
         t_freeze = eos.tfreeze(s_new,p)
         dt = (t_new - t_freeze)
         #print 't-tf min/max',dt.min(),dt.max()

         if direction == 1:
           oht = - (vhnew*cp*rho0*dt).sum() # watts
         else:
           oht =  (vhnew*cp*rho0*dt).sum()

         return oht

def create_ncfile(exp_name, xx, yy, ocean_time, args): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'_heat_transports.nc','w',format='NETCDF4')
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
   time.description = 'time since start of simulation'
   time.long_name = 'time since start of simulation'

   OnHT = ncfile.createVariable('OnHT',np.dtype('float32').char,('time'))
   OnHT.units = 'Watts'; OnHT.description = 'Onshore heat transport at ice shelf front'

   OfHT = ncfile.createVariable('OfHT',np.dtype('float32').char,('time'))
   OfHT.units = 'Watts'; OfHT.description = 'Offshore heat transport at ice shelf front'

   E = ncfile.createVariable('E',np.dtype('float32').char,('time'))
   E.units = 'nondim'; E.description = 'Entrainment Coeff'

   Qe = ncfile.createVariable('Qe',np.dtype('float32').char,('time'))
   Qe.units = 'm^3/s'; Qe.description = 'Entrainment volume flux'

   Qin = ncfile.createVariable('Qin',np.dtype('float32').char,('time'))
   Qin.units = 'm^3/s'; Qin.description = 'Total vol. flux in, across ice/ocean bnd'

   Qout = ncfile.createVariable('Qout',np.dtype('float32').char,('time'))
   Qout.units = 'm^3/s'; Qout.description = 'Total vol. flux in, across ice shelf front'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
