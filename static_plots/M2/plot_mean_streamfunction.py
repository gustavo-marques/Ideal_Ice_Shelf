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

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def get_data(cs,h,depth,z,exp,im):
    s1=Dataset(exp+'/prog.nc')
    s=Dataset(exp+'/ocean_month.nc')
    print 'Reading file', exp
    tmp = len(s.variables['time'][:])
    if tmp>24:
        v=s.variables['vh'][-24::,:,:,:].sum(axis=3).mean(axis=0)
    else:
        v=s.variables['vh'][:,:,:,:].sum(axis=3).mean(axis=0)

    temp=s.variables['temp'][-2::,:,:,:].mean(axis=3).mean(axis=0)
    salt=s.variables['salt'][-2::,:,:,:].mean(axis=3).mean(axis=0)
    rho = eos.wright_eos(temp,salt,2.0e7)
    h0=s1.variables['h'][0:1,:,:,im].mean(axis=0)
    e0=s1.variables['e'][0:1,:,:,im].mean(axis=0)
    s.close()

    km, jm = h0.shape
    # new arrays
    vh = np.zeros((km,jm)); psi = np.zeros((km,jm))
    # v at h points
    vtmp = 0.5 * (v[:,0:-1] + v[:,1::]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    vh[:,1::] = vtmp; vh[:,0] = 0.5*v[:,0] #v_j=1 = 0.5*v_(j=3/2)
    # mask and streamfunction

    #vh = np.ma.masked_where(h0<0.001, vh)
    rho = np.ma.masked_where(h0<0.001, rho)
    # integrate from bottom to top
    #psi[:,:] = -np.cumsum(vh[::-1,:],axis=0)[::-1,:]*0.5
    # top to bottom
    psi[:,:] = np.cumsum(vh,axis=0)*0.5
    # interpolate psi (can't remap it, it has been integrated already)
    psi1 = np.zeros((len(h),jm))
    for j in range(jm):
       if h0[:,j].sum() > 0.001:
          tmp = -0.5*(e0[0:-1,j]+e0[1::,j])
          psi1[:,j] = np.interp(h, tmp, psi[:,j])

       else:
          psi1[:,j] = -1.0e32

    #remap rho
    rho1=remap(cs,h0,e0,h,rho,z,depth)
    # mask bad values
    psi1 = np.ma.masked_where(np.ma.getmask(rho1), psi1)

    return rho1,psi1/1.0e6 # in sv


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
   os.system('mkdir OSF')
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   if os.path.exists(path+'/out1/'):
     x = Dataset(path+'/out1/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/out1/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/out1/ocean_geometry.nc').variables['D'][:,im/2]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:,im/2]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]

   e0 = Dataset(path+'/'+args.out+'/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)
   # remapping params
   cs = mom_remapping.Remapping_Cs()
   cs.remapping_scheme = 2
   cs.degree = 2
   z = np.linspace(0,depth.max(),1000)
   h = 0.5* (z[0:-1]+z[1::])
   # grid
   Y,Z = np.meshgrid(y[:,0],h)

   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   psi_vals = np.linspace(-0.1,0.1,20)
   norm = MidpointNormalize(midpoint=0)
   rho,psi=get_data(cs,h,depth,z,path+'/'+args.out,im/2)
   fig = plt.figure(figsize=(10,8),facecolor='white')
   ax = fig.add_subplot(111,axisbg='black')
   #ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.bwr, norm=norm)
   ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.bwr,extend='both')
   cb=plt.colorbar(ct, orientation='horizontal',ticks=[-0.10,-0.05,0.0,0.05,0.10])
   cb.set_label('Overturning streamfunction [sv]', fontsize=18)
   #ax.contour(Y,-Z,np.abs(psi),5,colors='gray',linewidth=0.2)
   #ax.set_aspect('equal',adjustable = 'box-forced')
   ax.set_xlim([0,1000])
   ax.set_ylim([-4000,0])
   #ct.set_clim(dyelev.min(),dyelev.max())
   #psi_vals = np.arange(psi.min(),psi.max(),2)
   s=ax.contour(Y,-Z,rho-1000,[36.90,37.0,37.10,37.2],colors='k',lw=0.5)
   ax.clabel(s, inline=1, fontsize=10,fmt='%4.2f',)
   # fill IS and bottom
   ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
   ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
   ax.plot(Y[0,:], ssh,'k-',lw=1.5)
   ax.set_xlabel('y [km]', fontsize=18)
   ax.set_ylabel('depth [m]', fontsize=18)
   plt.savefig('OSF/'+args.exp+'_dx'+args.dx+'_streamfunction.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()

   #fig = plt.figure(facecolor='white')
   #ax = fig.add_subplot(111,axisbg='black')
   #ct=ax.contourf(Y,-Z,rho-1000,np.linspace(37,37.3,50),)
   #plt.colorbar(ct, orientation='horizontal',ticks=[37,37.05,37.1,37.15,37.2,37.25,37.3])
   #ax.contour(Y,-Z,np.abs(psi),8,colors='gray',linewidth=0.2)
   ##ax.set_aspect('equal',adjustable = 'box-forced')
   #ax.set_xlim([0,1000])
   #ax.set_ylim([-4000,0])
   #s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',lw=0.5)
   #ax.clabel(s, inline=1, fontsize=10,fmt='%4.2f',)
   ## fill IS and bottom
   #ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
   #ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
   #ax.plot(Y[0,:], ssh,'k-',lw=1.5)
   #ax.set_title('Potential density')
   #ax.set_xlabel('y [km]', fontsize=18)
   #ax.set_ylabel('depth [m]', fontsize=18)
   ##plt.savefig('fig_tmp.png',format='png',dpi=300,bbox_inches='tight')
   #plt.close()

   print 'Saving streamfunction stats...'
   os.system('mkdir TXT')
   tmp1 = np.nonzero(ssh<=-200.0)[0][-1]
   print 'Ice shelf ends at y =',y[tmp1,0]
   text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_streamfunction.txt', "w")
   text_file.write("%f \n" % np.abs(psi[Y<=y[tmp1,0]]).max())
   text_file.write("%f \n" % np.abs(psi[Y<=y[tmp1,0]]).mean())
   text_file.close()

   print 'Mean psi:', np.abs(psi[Y<=y[tmp1,0]]).mean()
   print 'Max psi:', np.abs(psi[Y<=y[tmp1,0]]).max()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
