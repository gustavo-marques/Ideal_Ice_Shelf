from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
from matplotlib.colors import Normalize
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import sys, os
sys.path.append('../../')
from misc import *
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

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp):
    s=Dataset(exp+'/ocean_month.nc')
    print 'Reading file', exp
    tmp = len(s.variables['time'][:])
    if tmp>24:
        v=s.variables['vhbt'][-24::,:,:].mean(axis=0)
        u=s.variables['uhbt'][-24::,:,:].mean(axis=0)
        melt=s.variables['melt'][-24::,:,:].mean(axis=0)
    else:
        v=s.variables['vhbt'][:,:,:].mean(axis=0)
        u=s.variables['uhbt'][:,:,:].mean(axis=0)
        melt=s.variables['melt'][:,:,:].mean(axis=0)
    s.close()

    uh = np.zeros(u.shape); vh = np.zeros(v.shape); psi = np.zeros(u.shape)
    utmp = 0.5 * (u[:,0:-1] + u[:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
    vtmp = 0.5 * (v[0:-1,:] + v[1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    uh[:,1::] = utmp; uh[:,0] = 0.5*u[:,0] #u_i=1 = 0.5*u_(i=3/2)
    vh[1::,:] = vtmp; vh[0,:] = 0.5*v[0,:] #v_j=1 = 0.5*v_(j=3/2)
    psi[:,:] = (-uh[:,:].cumsum(axis=0) + vh[:,:].cumsum(axis=1))*0.5
    
    return psi/1.0e6, melt # in sv

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   # plot some metrics for runs with varing wind forcing
   os.system('mkdir PSI')
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

   if os.path.exists(path+'/M2_exp0/'):
     x = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/M2_exp0/ocean_geometry.nc').variables['D'][:]
     ssh = Dataset(path+'/M2_exp0/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]

   e0 = Dataset(path+'/M2_exp0/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)

   ### Call the function make_cmap which returns your colormap
   colors = [(255,255,255), (0,0,255), (255,0,0)] 
   my_cmap = make_cmap(colors, bit=True)

   # read data
   psi_vals = np.linspace(-0.21,0.1,15)
   melt_vals = np.linspace(-0.1,15,20)
   exps = ['M2_exp0','M2_exp13','M2_exp4','M2_exp14']
   lbs = ['a) $U_{shelf}$ = -5.0 m s$^{-1}$, melting on','b) $U_{shelf}$ = -5.0 m s$^{-1}$, melting off','c) $U_{shelf}$ = 5.0 m s$^{-1}$, melting on','d) $U_{shelf}$ = 5.0 m s$^{-1}$, melting off']

   # loop over experiments
   for i in range(len(exps)):
     # get data  
     psi,melt=get_data(path+'/'+exps[i])
     print 'melt min/max',melt.min(), melt.max()
     #psi = np.ma.masked_where(depth<500,psi)
     tmp1 = np.nonzero(ssh<=-2.0)[0][-1]
     psi_ice_shelf = psi[0:tmp1,:]
     print 'psi min/max',psi_ice_shelf.min(), psi_ice_shelf.max()
     # plot
     fig, ax = plt.subplots()
     cs=ax.contourf(x[0:tmp1,:],y[0:tmp1,:],melt[0:tmp1,:],melt_vals,cmap=my_cmap,extend='min')
     ax.set_axis_bgcolor('gray')
     #ax.set_title(lbs[i],fontsize=18)
     cs.set_clim(melt_vals.min(),melt_vals.max())
     cc=ax.contour(x[0:tmp1,:],y[0:tmp1,:],psi_ice_shelf,psi_vals,colors='k',linewidths=0.5)
     ax.clabel(cc, fontsize=10, inline=1,fmt='%1.2f',)
     if (i == 0 or i == 2):
       y_fill = np.zeros(im)  
       for ii in range(im):
         tmp2 = np.nonzero(melt.mask[:,ii]==True)[0][-1]
         if (y[tmp2,ii] <= 100.):
           tmp3 = np.nonzero(melt[:,ii]>0.0)[0][0]
           y_fill[ii] = y[tmp3,ii]
         else:
           y_fill[ii] = y[tmp2,ii]

     #ax.fill_between(x[0,:],y_fill,0, facecolor='#A9A9A9', interpolate=True)
     ax.fill_between(x[0,:],y_fill,0, facecolor='peru', interpolate=True)
     
     ax.set_xlim([0,500])
     ax.set_ylim([0,y[0:tmp1,:].max()])
     ax.set_ylabel('y [km]',fontsize=16)
     ax.set_xlabel('x [km]', fontsize=16)

     fig.subplots_adjust(bottom=0.28)
     cbar_ax = fig.add_axes([0.1, 0.11, 0.8, 0.05])
     cbar=fig.colorbar(cs, orientation='horizontal', cax=cbar_ax,ticks=[-0.1,2.5,5.0,7.5,10.0,12.5,15.0])
     #cbar=fig.colorbar(cs, orientation='horizontal', cax=cbar_ax)
     cbar.set_label(r'Melt rate [m yr$^{-1}$]', fontsize=16)
     plt.savefig('PSI/'+exps[i]+'_barotropic_melt_dx1.png',format='png',dpi=300,bbox_inches='tight')
     plt.show()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
