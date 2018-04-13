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

  parser.add_argument('-out', type=str, default='out6',
      help='''Name of output file (default = out6).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def get_data(exp):
    s1=Dataset(exp+'/prog.nc')
    s=Dataset(exp+'/ocean_month.nc')
    print 'Reading file', exp
    tmp = len(s.variables['time'][:])
    if tmp>24:
        v=s.variables['vhbt'][-24::,:,:].mean(axis=0)
        u=s.variables['uhbt'][-24::,:,:].mean(axis=0)
    else:
        v=s.variables['vhbt'][:,:,:].mean(axis=0)
        u=s.variables['uhbt'][:,:,:].mean(axis=0)
    #v=s.variables['vhbt'][-12::,:,:].mean(axis=0)
    #u=s.variables['uhbt'][-12::,:,:].mean(axis=0)
    s.close()

    uh = np.zeros(u.shape); vh = np.zeros(v.shape); psi = np.zeros(u.shape)
    utmp = 0.5 * (u[:,0:-1] + u[:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
    vtmp = 0.5 * (v[0:-1,:] + v[1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    uh[:,1::] = utmp; uh[:,0] = 0.5*u[:,0] #u_i=1 = 0.5*u_(i=3/2)
    vh[1::,:] = vtmp; vh[0,:] = 0.5*v[0,:] #v_j=1 = 0.5*v_(j=3/2)
    psi[:,:] = (-uh[:,:].cumsum(axis=0) + vh[:,:].cumsum(axis=1))*0.5
    
    return psi/1.0e6 # in sv


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
   os.system('mkdir PSI')
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   if os.path.exists(path+'/out1/'):
     x = Dataset(path+'/out1/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/out1/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/out1/ocean_geometry.nc').variables['D'][:]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]

   e0 = Dataset(path+'/'+args.out+'/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)

   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   psi_vals = np.linspace(-0.5,40,50)
   norm = MidpointNormalize(midpoint=0)
   psi=get_data(path+'/'+args.out)
   psi = np.ma.masked_where(depth<500,psi)
   tmp1 = np.nonzero(ssh<=-150.0)[0][-1]
   psi_ice_shelf = psi[0:tmp1,:]
 
   fig = plt.figure(figsize=(10,8),facecolor='white')
   ax = fig.add_subplot(111,axisbg='black')
   #ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.bwr, norm=norm)
   #ct=ax.contourf(x,y,psi,psi_vals,cmap=plt.cm.bwr,extend='both')
   ct=ax.contourf(x[0:tmp1,:],y[0:tmp1,:],psi_ice_shelf,cmap=plt.cm.bwr,extend='both')
   cb=plt.colorbar(ct, orientation='horizontal')
   cb.set_label('Barotropic streamfunction [sv]', fontsize=18)
   #ax.contour(Y,-Z,np.abs(psi),5,colors='gray',linewidth=0.2)
   #ax.set_aspect('equal',adjustable = 'box-forced')
   ax.set_xlim([0,500])
   ax.set_ylim([0,200])
   plt.savefig('PSI/'+args.exp+'_dx'+args.dx+'_barotropic_ice_shelf.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()

   fig = plt.figure(figsize=(10,8),facecolor='white')
   ax = fig.add_subplot(111,axisbg='black')
   #ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.bwr, norm=norm)
   ct=ax.contourf(x,y,psi,psi_vals,cmap=plt.cm.bwr,extend='both')
   cb=plt.colorbar(ct, orientation='horizontal')
   cb.set_label('Barotropic streamfunction [sv]', fontsize=18)
   ax.set_xlim([0,500])
   ax.set_ylim([0,1000])
   plt.savefig('PSI/'+args.exp+'_dx'+args.dx+'_barotropic.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()

   print 'Saving streamfunction stats...'
   os.system('mkdir TXT')
   print 'psi amplitude (max-min)',psi_ice_shelf.max()-psi_ice_shelf.min()
   text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_barotropic_amplitude.txt', "w")
   text_file.write("%f \n" % (psi_ice_shelf.max()-psi_ice_shelf.min()))
   text_file.close()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
