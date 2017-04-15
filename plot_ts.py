#!/usr/bin/env python

# Gustavo Marques

from remapping import mom_remapping
import matplotlib.colors as colors
from matplotlib.ticker import LogFormatter
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib
import os
import argparse
import numpy
from misc import *
import wright_eos as eos

matplotlib.rcParams.update({'font.size': 16})

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Plot t/s diagrams using passive tracers (tr1 and tr2) to tracer water masses.
      ''',
  epilog='Written by Gustavo Marques, April 2017.')

  parser.add_argument('-name', type=str, default='test', help='''Name of the experiment (default is test).''')

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
  if not os.path.exists('PNG'):
       os.system('mkdir -p PNG')

  # static variables
  depth = netCDF4.Dataset('ocean_geometry.nc').variables['D'][:]
  X = netCDF4.Dataset('ocean_geometry.nc').variables['geolon'][:]
  Y = netCDF4.Dataset('ocean_geometry.nc').variables['geolat'][:]
  area = netCDF4.Dataset('ocean_geometry.nc').variables['Ah'][:] # m^2

  # time
  tind = np.arange(ti,tf,n)
  tm = len(tind)

  #color scale
  #colors = [(255,255,255), (0,0,255), (255,0,0)]
  colors = [(255,255,255), (0,0,255), (51,255,51), (255,255,51), (255,0,0), (255,0,255)]
  my_cmap = make_cmap(colors, bit=True)

  # main loop
  for t in range(tm):
      time = netCDF4.Dataset('prog.nc').variables['time'][tind[t]] # in days
      print 'Time is:',time
      tracer1 = netCDF4.Dataset('prog.nc').variables['tr1'][tind[t],:,:,:]
      tracer2 = netCDF4.Dataset('prog.nc').variables['tr2'][tind[t],:,:,:]
      h = netCDF4.Dataset('prog.nc').variables['h'][tind[t],:,:,:]
      volume = h * np.tile(area, (h.shape[0],1,1))
      e = netCDF4.Dataset('prog.nc').variables['e'][tind[t],:,:,:]
      depth = -0.5 * (e[0:-1,:,:] + e[1::,:,:])
      temp = netCDF4.Dataset('prog.nc').variables['temp'][tind[t],:,:,:]
      salt = netCDF4.Dataset('prog.nc').variables['salt'][tind[t],:,:,:]
      rhopot2 = netCDF4.Dataset('prog.nc').variables['rhopot2'][tind[t],:,:,:]
      # mask vanished layers
      temp = np.ma.masked_where(h<0.01,temp)
      salt = np.ma.masked_where(h<0.01,salt)
      tracer1 = np.ma.masked_where(h<0.01,tracer1)
      tracer2 = np.ma.masked_where(h<0.01,tracer2)
      volume = np.ma.masked_where(h<0.01,volume)

      TS_diagram(temp,salt,rhopot2,depth,tracer1,tracer2,volume,t,my_cmap,args.name)

      print '\n'

  print 'Done with loop!'
  return

def TS_diagram(temp,salt,rho,depth,tracer1,tracer2,volume,t,cm,exp):
    """ plot a TS diagram
    salt,temp,depth (positive) are obvious
    rho is sigma2 at cell center
    t = counter
    """
    n = 10
    temp = temp[::n,::n,:].mean(axis=2); salt = salt[::n,::n,:].mean(axis=2);
    depth = depth[::n,::n,:].mean(axis=2)
    tracer1 = tracer1[::n,::n,:].mean(axis=2)
    tracer2 = tracer2[::n,::n,:].mean(axis=2)
    volume = volume[::n,::n,:].mean(axis=2)

    # Figure out boudaries (mins and maxs)
    smin = 34.05  #np.min(salts) - (0.01 * np.min(salts))
    smax = 34.8 #np.max(salts) + (0.01 * np.max(salts))
    tmin = -2.5 #np.min(temps) - (0.1 * np.min(temps))
    tmax = 1.6 #np.max(temps) + (0.1 * np.max(temps))
    temps=np.linspace(tmin,tmax,10)
    salts=np.linspace(smin,smax,10)
    tf = eos.tfreeze(salts,0)
    S,T = np.meshgrid(salts,temps)
    sig2=eos.wright_eos(T,S,2e7)-1000

    # normalize area, size of points are proportional to area
    size = (volume/volume.max())* 150.
    size[size<=0.001] = 0.

    rho=rho-1000
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(111)
    ax1.set_aspect('auto')
    rho_vals = [36.9,37,37.1,37.2,37.3]
    CS = plt.contour(S,T,sig2,rho_vals, colors='k',lw=0.25)
    plt.clabel(CS, fontsize=15, inline=1, fmt='%4.2f') # Label every second level
    #sc=ax1.scatter(salt,temp, c=depth, s=1,vmin=0,vmax=depth.max(), cmap=cm,norm=mp.colors.LogNorm())
    sc=ax1.scatter(salt.flatten(),temp.flatten(), c=depth.flatten(), vmin = 5, vmax = depth.max(),s=size, alpha = 1, norm=colors.LogNorm(), marker='o',cmap = cm)#, s=4,vmin=0,vmax=depth.max(), cmap=cm, norm=colors.Normalize())
    cb=fig.colorbar(sc)
    cb.set_label('Depth [m]',fontsize = 18)
    aux(ax1,tmin,tmax,smin,smax,volume,size,salts,tf)
    plt.tight_layout()
    s = str("plt.savefig('PNG/%s-ts-depth-%04d.png',bbox_inches='tight')"% (exp,t))
    eval(s)

    fig = plt.figure(figsize=(10,8))
    ax2 = fig.add_subplot(111)
    ax2.set_aspect('auto')
    CS2 = plt.contour(S,T,sig2,rho_vals, colors='k',lw=0.25)
    plt.clabel(CS2, fontsize=15, inline=1, fmt='%4.2f') # Label every second level
    sc2=ax2.scatter(salt.flatten(),temp.flatten(), c=tracer1.flatten(), cmap=plt.cm.bwr,vmin = 0.0001, vmax = 1.5,s=size, alpha = 1,  norm=colors.LogNorm(), marker='o')#, s=4,vmin=0,vmax=depth.max(), cmap=cm, norm=colors.Normalize())
    cb2=fig.colorbar(sc2)
    cb2.set_label('Tracer (brine rejection)',fontsize = 18)
    aux(ax2,tmin,tmax,smin,smax,volume,size,salts,tf)
    plt.tight_layout()
    s = str("plt.savefig('PNG/%s-ts-tr1-%04d.png',bbox_inches='tight')"% (exp,t))
    eval(s)

    fig = plt.figure(figsize=(10,8))
    ax2 = fig.add_subplot(111)
    ax2.set_aspect('auto')
    CS2 = plt.contour(S,T,sig2,rho_vals, colors='k',lw=0.25)
    plt.clabel(CS2, fontsize=15, inline=1, fmt='%4.2f') # Label every second level
    sc2=ax2.scatter(salt.flatten(),temp.flatten(), c=tracer2.flatten(), cmap=plt.cm.bwr,vmin = 0.0001, vmax = 50,s=size, alpha = 1,  norm=colors.LogNorm(), marker='o')#, s=4,vmin=0,vmax=depth.max(), cmap=cm, norm=colors.Normalize())
    cb2=fig.colorbar(sc2)
    ax2.plot(salts,tf,'gray')
    cb2.set_label('Tracer (melt water)',fontsize = 18)
    aux(ax2,tmin,tmax,smin,smax,volume,size,salts,tf)
    plt.tight_layout()
    s = str("plt.savefig('PNG/%s-ts-tr2-%04d.png',bbox_inches='tight')"% (exp,t))
    eval(s)
    plt.show()
#    plt.close('all')
    return

def aux(ax,tmin,tmax,smin,smax,vol,size,salts,tf):
    ax.set_ylim(tmin,tmax)
    ax.set_xlim(smin,smax)
    ax.set_xlabel('Salinity', fontsize = 18)
    ax.set_ylabel(r'Temperature [$^o$C]', fontsize = 18)
    s = str("ax.text(34.12,0.45,r'= %2.1f  km$^3$')"% (vol.max()/1.0e9))
    eval(s)
    ax.plot(salts,tf,'gray')
    ax.scatter(34.1,0.5,s=size.max(),marker='o',c='k')
    ax.text(34.71,-1.8,r'$T_f$',color='k')

    return

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()

