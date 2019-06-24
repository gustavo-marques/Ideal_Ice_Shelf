from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import sys, os
sys.path.append('../../')
from nc2VTK import VTKgen
from misc import *
import argparse

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Average T,S,U,V,RHO and Tr2 over last two years of simulation 
      ''',
  epilog='Written by Gustavo Marques, Jun. 2018.')

  parser.add_argument('-dx', type=str, default='1',
      help='''Horizontal grid resolution (default = 1).''')

  parser.add_argument('-exp', type=str, default='M2_exp0',
                help='''Experiment name (default = M2_exp0).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp):
    print 'Reading file', exp
    s=Dataset(exp+'/ocean_month.nc')
    tmp = len(s.variables['time'][:])
    if tmp>24:
        h=s.variables['h'][-24::,:,:,:].mean(axis=0)
        e=s.variables['e'][-24::,:,:,:].mean(axis=0)
        v=s.variables['v'][-24::,:,:,:].mean(axis=0)
        u=s.variables['u'][-24::,:,:,:].mean(axis=0)
        temp=s.variables['temp'][-24::,:,:,:].mean(axis=0)
        salt=s.variables['salt'][-24::,:,:,:].mean(axis=0)
        rho=s.variables['rhopot2'][-24::,:,:,:].mean(axis=0)
        tr2=s.variables['tr2'][-24::,:,:,:].mean(axis=0)
    else:
        h=s.variables['h'][:,:,:,:].mean(axis=0)
        e=s.variables['e'][:,:,:,:].mean(axis=0)
        v=s.variables['v'][:,:,:,:].mean(axis=0)
        u=s.variables['u'][:,:,:,:].mean(axis=0)
        temp=s.variables['temp'][:,:,:,:].mean(axis=0)
        salt=s.variables['salt'][:,:,:,:].mean(axis=0)
        rho=s.variables['rhopot2'][:,:,:,:].mean(axis=0)
        tr2=s.variables['tr2'][:,:,:,:].mean(axis=0) 
    s.close()

    uh = np.zeros(salt.shape); vh = np.zeros(salt.shape) 
    utmp = 0.5 * (u[:,0:-1] + u[:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
    vtmp = 0.5 * (v[0:-1,:] + v[1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    uh[:,1::] = utmp; uh[:,0] = 0.5*u[:,0] #u_i=1 = 0.5*u_(i=3/2)
    vh[1::,:] = vtmp; vh[0,:] = 0.5*v[0,:] #v_j=1 = 0.5*v_(j=3/2)
    
    return temp, salt, rho, tr2, uh, vh, h, e 

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   # plot some metrics for runs with varing wind forcing
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

   lonh = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['lonh'][:]
   lath = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['lath'][:]
   lons, lats = np.meshgrid(lonh,lath)
   D = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['D'][:]
   h=Dataset(path+'/M2_exp0/ocean_month.nc').variables['h'][0,:,:,:]
   e=Dataset(path+'/M2_exp0/ocean_month.nc').variables['e'][0,:,:,:]
   # correct top and bottom, for vis pourposes
   h[0,:,:]=e[0,:,:]; h[-1,:,:]=e[-1,:,:]
   D=np.ma.masked_where(D <= 1, D)
   D.mask = np.ma.array(D); D.mask[:,:]=False
   # create bathymetry
   VTKgen(lats,lons,D.mask,depth=D,h=h,fname='M2_1km')

   # list experiments
   #exps = ['M2_exp0','M2_exp13','M2_exp4','M2_exp14']

   temp, salt, rho, tr2, uh, vh, h, e = get_data(path+'/'+args.exp)
   hz=0.5*(e[0:-1,:,:]+e[1::,:,:])# cell depth

   # vel
   print 'Saving vel... \n'
   VTKgen(lats,lons,D.mask,h=hz,u=uh,v=vh,fname=args.exp,varname='vel',t=9999)
   # temp
   print 'Saving temp... \n'
   VTKgen(lats,lons,D.mask,h=hz,tracer=temp,fname=args.exp,varname='temp',t=9999)
   # salt
   print 'Saving salt... \n'
   VTKgen(lats,lons,D.mask,h=hz,tracer=salt,fname=args.exp,varname='salt',t=9999)
   # rhopot2
   print 'Saving rhopot2... \n'
   VTKgen(lats,lons,D.mask,h=hz,tracer=rho,fname=args.exp,varname='rhopot2',t=9999)
   # tr2
   print 'Saving tr2... \n'
   VTKgen(lats,lons,D.mask,h=hz,tracer=tr2,fname=args.exp,varname='tr2',t=9999)

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
