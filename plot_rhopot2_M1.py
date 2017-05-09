#!/usr/bin/env python

# generate 3D diagnostics
# Gustavo Marques, Sep. 2016

import argparse
from netCDF4 import Dataset
import numpy as np
from remapping import mom_remapping
import numpy
import matplotlib.pyplot as plt
import warnings
import os,sys
import matplotlib
from misc import *

matplotlib.rcParams.update({'font.size': 22})

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Remap data into a z-grid then plot a section based on the parameters specified by the user.
      Ex: plot_rhopot2_M1.py -min 37 -max 37.3 -n M1_exp13_1km -var_name rhopot2 -i 250
      ''',
  epilog='Written by Gustavo Marques, Mar. 2017.')

  parser.add_argument('-n', type=str, default='Test',
      help='''The name of the experiment. Default is Test.''')

  parser.add_argument('-var_name', type=str, default='temp',
      help='''Variable to be plotted. Default is temp.''')

  parser.add_argument('-ti', type=int, default=0,
      help='''The initial time indice to save the data. Default value is 0.''')

  parser.add_argument('-tf', type=int, default=-1,
      help='''The final time indice to save the data. Default value is -1(len(time)).''')

  parser.add_argument('-dt', type=int, default=1,
      help='''The time indice interval to save the data. Default value is 1, which saves the entire dataset.''')

  parser.add_argument('-i', type=int, default=0,
      help='''i-th indice where to take section. Default value is 0, which means all data.''')

  parser.add_argument('-j', type=int, default=0,
      help='''j-th indice where to take section. Default value is 0, which means all data.''')

  parser.add_argument('-min', type=float, default=-2.2,
      help='''Minimum value in the colorbar range. Default value is -2.2.''')

  parser.add_argument('-max', type=float, default=1.2,
      help='''Maximum value in the colorbar range. Default value is 1.2.''')

  parser.add_argument('--oceanfile', type=str, default='prog.nc',
      help='''Name of the netCDF file with the ocean variables. Default is prog.nc.''')

  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver(optCmdLineArgs)

def driver(args):
    dt = args.dt
    ti = args.ti
    if args.tf == -1:
       tf = len(Dataset(args.oceanfile).variables['time'][:])
    else:
       tf = args.tf

    # create dir locally
    os.system('mkdir -p PNG')
    # time
    tind = np.arange(ti,tf,dt)
    tm = len(tind)
    # remapping params
    cs = mom_remapping.Remapping_Cs()
    cs.remapping_scheme = 2
    cs.degree = 2

    # ocean grid
    depth = Dataset('ocean_geometry.nc').variables['D'][:,:]
    #h=Dataset(args.oceanfile).variables['h'][0,:,:,:]
    NZ,NY,NX=(Dataset(args.oceanfile).variables['h'][0,:,:,:]).shape
    z = np.linspace(0,depth.max(),2000)
    h = 0.5* (z[0:-1]+z[1::])

    if args.j > 0:
       x = Dataset('ocean_geometry.nc').variables['lonh'][:]
       i = range(len(x))
       j = args.j
       args.axis = 'x'
       XY, ZZ = np.meshgrid(x,h)
       pos = 'j'+str(j)
    else:
       y = Dataset('ocean_geometry.nc').variables['lath'][:]
       j = range(len(y))
       i = args.i
       XY, ZZ = np.meshgrid(y,h)
       args.axis = 'y'
       pos = 'i'+str(i)

    # depth at i,j
    depth = depth[j,i]
    ssh = Dataset('IDEAL_IS_IC.nc').variables['ave_ssh'][0,j,i]
    for t in range(len(tind)):
       get_data(tind[t],XY,ZZ,h,ssh,i,j,z,cs,depth,pos,args)

    print ' \n' + '==> ' + '  DONE!\n' + ''

    return

def get_data(t, XY, ZZ, h, ssh, i, j, z,cs, depth,pos,args):
       """
       Reads data from netcdf at time indice t, remaps then plots.
       """
    # check if data has been saved
    #path_to_file = str('PNG/%s-%05d.vtk' % (args.n,t))
    #if os.path.isfile(path_to_file):
    #   print ' \n' + '==> ' + 'FILE EXISTS, MOVING TO THE NEXT ONE ...\n' + ''
    #else:
       print 'Saving time indice:', t
       # read data
       print 'Reading data...'
       time=Dataset(args.oceanfile).variables['time'][t]
       e0=Dataset(args.oceanfile).variables['e'][t,:,j,200:400]
       h0=Dataset(args.oceanfile).variables['h'][t,:,j,200:400]
       var0=Dataset(args.oceanfile).variables[args.var_name][t,:,j,200:400]
       rho0=Dataset(args.oceanfile).variables['rhopot2'][t,:,j,200:400]
       units=Dataset(args.oceanfile).variables[args.var_name].units
       # average
       h0 = h0.mean(axis = 2); var0 = var0.mean(axis = 2)
       e0 = e0.mean(axis = 2); rho0 = rho0.mean(axis = 2)
       rho0 = rho0 - 1000.
       colors = [(255,255,255), (0,0,255), (51,255,51), (255,255,51), (255,0,0), (255,0,255)]
       my_cmap = make_cmap(colors, bit=True)
       #remap
       var1=remap(cs,h0,e0,h,var0,z,depth)
       rho1=remap(cs,h0,e0,h,rho0,z,depth)
       print 'Plotting data... \n'
       saltlev = np.linspace(args.min,args.max,50)
       fig = plt.figure(facecolor='white',figsize = (12,10))
       ax = fig.add_subplot(111,axisbg='white')
       cs = ax.contourf(XY,-ZZ,rho1,saltlev,cmap = my_cmap,extend='both')
       cs.set_clim(saltlev.min(),saltlev.max())
       cb = plt.colorbar(cs,orientation='horizontal',ticks=[37,37.05, 37.1, 37.15,37.2,37.25,37.3],extend='both')
       cs = ax.contour(XY,-ZZ,rho1,[37.205, 37.03],colors='k')
       #ax.clabel(cs, inline=1, fontsize=18)
       rho1 = np.ma.masked_where(XY>500., rho1)
       cs = ax.contour(XY,-ZZ,rho1,[37.24],colors='k')
       #ax.clabel(cs, inline=1, fontsize=18)

       #cs = ax.contourf(XY,-ZZ,var1,dyelev,extend='both')
       ax.set_aspect('auto')
       #ax.fill_between(XY[0,:], e0[0,:], 0.0, where=y2 <= y1, facecolor='red', interpolate=True)
       ax.fill_between(XY[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
       ax.fill_between(XY[0,:], e0[-1,:], -4000, facecolor='gray', interpolate=True)
       ax.plot(XY[0,:], ssh,'k-',lw=1.5)
       #ax.plot(XY[0,:],-depth,'k',lw=1.5)
       #ax.plot(XY[0,:],e0[0,:],'k',lw=1.5)
       #cb = plt.colorbar(cs,orientation='horizontal')
       cb.set_label(r'$\sigma_2$ [kg m$^{-3}$]', fontsize=22)
       #ax.plot(xpall,zpall,'k-',lw=0.5)
       #ax.contour(X,Y,depth,[550,650,750,850],colors='k')
       ax.set_xlabel(args.axis+' [km]', fontsize = 22)
       ax.set_ylabel('depth [m]',fontsize = 22)
       ax.set_xlim(XY.min(),XY.max())
       ax.set_ylim(-depth.max(),0)
       #ss = str("ax.set_title('Time: %5.1f days')"% (time))
       ##eval(ss)
       #plt.grid()
       ss = str("plt.savefig('PNG/%s-%s-%s-%03d.png',bbox_inches='tight')"% (args.n,args.var_name,pos,t))
       eval(ss)
       plt.show()
       plt.close()

       return

def remap(cs,h0,e0,h,data0,z,depth):
   """
   remaps data...
   """
   km,im = h0.shape
   data1 = np.ma.zeros((len(h),im))
   for i in range(im):
     h1 = np.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
        h1[h<-e0[0,i]] = 0.0 # top
        h1[h>depth[i]] = 0.0 # bottom
        # need to account for SSH and make sure total thicknesses are
        # the same
        dh = h0[:,i].sum() - h1.sum()
        tmp1 = np.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
          if dh > 0.0:
             # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = np.ma.masked

     # for debugging
        #if h0[:,i].sum() != h1.sum():
        #   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
        data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
        # mask
        data1[h1==0.0,i] = np.ma.masked;
     else: # land/iceshelf
        data1[:,i] = np.ma.masked

   return data1

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
