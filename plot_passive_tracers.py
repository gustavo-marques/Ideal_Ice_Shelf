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
import matplotlib

matplotlib.rcParams.update({'font.size': 22})

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

  parser.add_argument('-plot_vel', help='''Plot depth averaged velocity on top of tracer.''', action="store_true")

  parser.add_argument('-nvel', type=int, default='10', help='''Number of grid points to skip when plotting velocity. Default is 10.''')

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
         #dyelev = np.linspace(0,50,100)
         dyelev = np.linspace(0,50,100)
         #dyetic = [0,5,10,25,50]
         dyetic = [0,10.0,20,30.0,40.0,50.0]
      else:
         dyelev = np.linspace(0,1.5,100)
         dyetic = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4]

      if args.plot_vel:
        u = netCDF4.Dataset('prog.nc').variables['u'][tind[t],:,:,:]
        v = netCDF4.Dataset('prog.nc').variables['v'][tind[t],:,:,:]
        # u and v at h points
        uh = np.zeros(u.shape); vh = np.zeros(v.shape);
        utmp = 0.5 * (u[:,:,0:-1] + u[:,:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
        vtmp = 0.5 * (v[:,0:-1,:] + v[:,1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
        uh[:,:,1::] = utmp; uh[:,:,0] = 0.5*u[:,:,0] #u_i=1 = 0.5*u_(i=3/2)
        vh[:,1::,:] = vtmp; vh[:,0,:] = 0.5*v[:,0,:] #v_j=1 = 0.5*v_(j=3/2)
        u_ave = (uh * h0).sum(axis = 0)/h0.sum(axis=0)
        v_ave = (vh * h0).sum(axis = 0)/h0.sum(axis=0)
        # mask bad values
        u_ave = np.ma.masked_where(np.abs(u_ave)>10.,u_ave)
        v_ave = np.ma.masked_where(np.abs(v_ave)>10.,v_ave)
        print 'u max/min, v max/min', u_ave.max(), u_ave.min(), v_ave.max(), v_ave.min()
        plt_xy_uv(X,Y,data_ave,u_ave, v_ave, args.tracer,depth,time,t,dyelev,dyetic,my_cmap,args.name,args.nvel)
      else:
        plt_xy(X,Y,data_ave,args.tracer,depth,time,t,dyelev,dyetic,my_cmap,args.name)

      print '\n'

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

def plt_xy_uv(X,Y,tracer,u,v,tracer_name,depth,time,t,dyelev,dyetic,my_cmap,exp,n):
    ##fig = plt.figure(facecolor='black',figsize=(10,10))
    #fig = plt.figure(facecolor='black')
    ##ax = fig.add_subplot(111,axisbg='gray')
    fig, (ax, cax) = plt.subplots(nrows=2,figsize=(10,10),
                  gridspec_kw={"height_ratios":[1, 0.05]})
    #ax.axis('equal')
    cs = ax.contourf(X,Y,tracer,dyelev,cmap=my_cmap,extend='both')
    ax.set_aspect('equal')
    cs.set_clim(dyelev.min(),dyelev.max())
    #cbaxes = fig.add_axes([0.15, 0.05, 0.75, 0.03])
    cb = plt.colorbar(cs, cax = cax,orientation='horizontal',ticks=dyetic,extend='both')
    #cb = plt.colorbar(cs,orientation='horizontal',ticks=dyetic,extend='both')
    #cb = plt.colorbar(cs,orientation='horizontal')
    cb.set_label('Depth averaged tracer', fontsize=22)
    #ax.plot(xpall,zpall,'k-',lw=0.5)
    ax.contour(X,Y,depth,[550,650,750,850],colors='k')
    Q = ax.quiver(X[::n,::n],Y[::n,::n],u[::n,::n],v[::n,::n],scale=1)
    qk = ax.quiverkey(Q, 0.72, 0.25, 0.05, r'0.05 $\frac{m}{s}$', labelpos='E',
                   coordinates='figure')
    ax.set_xlabel('x [km]', fontsize = 22)
    ax.set_ylabel('y [km]',fontsize = 22)
    ax.set_xlim(0,500)
    ax.set_ylim(Y.min(),500)
    ss = str("ax.set_title('Time: %5.2f days')"% (time))
    eval(ss)
    #forceAspect(ax,aspect=1)
    #plt.grid()
    ss = str("plt.savefig('PNG/%s-%s-%03d.png',bbox_inches='tight')"% (exp,tracer_name,t))
    eval(ss)
    plt.close()

    return

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

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

