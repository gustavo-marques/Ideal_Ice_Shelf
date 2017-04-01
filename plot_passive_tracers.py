#!/usr/bin/env python

# Gustavo Marques

from remapping import mom_remapping
import matplotlib.colors as colors
from matplotlib.ticker import LogFormatter
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os
import argparse
import numpy
from misc import *

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Plot maps of passive tracers (tr1 and tr2).
      ''',
  epilog='Written by Gustavo Marques, Mar. 2017.')

  parser.add_argument('-name', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-tracer', type=str, default='tr2', help='''Tracer name in the netcdf file (default is tr2).''')

  parser.add_argument('-n', type=int, default='1', help='''Compute every n values (default is 1).''')

  parser.add_argument('-ti', type=int, default='0', help='''Initial time indice (default is 0).''')
  parser.add_argument('-tf', type=int, default='-1', help='''Final time indice (default is len(time)).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
  n = args.n # plot every n
  ti = args.ti
  if args.tf == -1:
     tf = len(netCDF4.Dataset('prog.nc').variables['time'][:])
  else:
     tf = args.tf

  # create dir locally
  os.system('mkdir -p PNG')

  # static variables
  depth = netCDF4.Dataset('ocean_geometry.nc').variables['D'][:]
  X = netCDF4.Dataset('ocean_geometry.nc').variables['geolon'][:]
  Y = netCDF4.Dataset('ocean_geometry.nc').variables['geolat'][:]
  ave_ssh = netCDF4.Dataset('IDEAL_IS_IC.nc').variables['ave_ssh'][0,:]
  z = np.linspace(0,depth.max(),100)
  h = 0.5* (z[0:-1]+z[1::])

  # time
  tind = np.arange(ti,tf,n)
  tm = len(tind)
  # remapping params
  cs = mom_remapping.Remapping_Cs()
  # PCM
  #cs.remapping_scheme = 0 # PCM
  # PPM-H4
  cs.remapping_scheme = 2
  cs.degree = 2

  #color scale
  #colors = [(255,255,255), (0,0,255), (255,0,0)]
  colors = [(255,255,255), (0,0,255), (51,255,51), (255,255,51), (255,0,0), (255,0,255)]
  my_cmap = make_cmap(colors, bit=True)

  # main loop
  for t in range(tm):
      time = netCDF4.Dataset('prog.nc').variables['time'][tind[t]] # in days
      print 'Time is:',time
      data0 = netCDF4.Dataset('prog.nc').variables[args.tracer][tind[t],:,:,:]
      h0 = netCDF4.Dataset('prog.nc').variables['h'][tind[t],:,:,:]
      e0 = netCDF4.Dataset('prog.nc').variables['e'][tind[t],:,:,:]
      # remap
      #data1 = remap(cs,h0,e0,h,data0,z,depth)
      #
      data_ave = (data0 * h0).sum(axis = 0)/h0.sum(axis=0)
      print 'min/max',data_ave.min(), data_ave.max()
      if args.tracer == 'tr2':
         dyelev = np.linspace(0,50,100)
         dyetic = [0,5,10,25,50]
      else:
         dyelev = np.linspace(0,1.5,100)
         dyetic = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4]

      plt_xy(X,Y,data_ave,args.tracer,depth,time,t,dyelev,dyetic,my_cmap,args.name)

  print 'Done with loop!'
  return

def plt_xy(X,Y,tracer,tracer_name,depth,time,t,dyelev,dyetic,my_cmap,exp):
    fig = plt.figure(facecolor='black')
    ax = fig.add_subplot(111,axisbg='gray')
    cs = ax.contourf(X,Y,tracer,dyelev,cmap=my_cmap,extend='both')
    ax.set_aspect('auto')
    cs.set_clim(dyelev.min(),dyelev.max())
    cb = plt.colorbar(cs,orientation='horizontal',ticks=dyetic,extend='both')
    #cb = plt.colorbar(cs,orientation='horizontal')
    cb.set_label(tracer_name, fontsize=16)
    #ax.plot(xpall,zpall,'k-',lw=0.5)
    ax.contour(X,Y,depth,[550,650,750,850],colors='k')
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_xlim(X.min(),X.max())
    ax.set_ylim(Y.min(),600)
    ss = str("ax.set_title('Time: %5.2f days')"% (time))
    eval(ss)
    #plt.grid()
    ss = str("plt.savefig('PNG/%s-%s-%03d.png')"% (exp,tracer_name,t))
    eval(ss)
    plt.close()

    return
def remap(cs,h0,e0,h,data0,z,depth):
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

