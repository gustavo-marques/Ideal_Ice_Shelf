import numpy
from netCDF4 import Dataset
import argparse
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import sys, os
sys.path.append('../../')
from misc import *
import wright_eos as eos
from remapping import mom_remapping

from computeOSF import computeOSF

parser = argparse.ArgumentParser(description=
      '''
      Generate files for Idealized Ice Shelf problem.
      ''',
epilog='Written by Gustavo Marques, Oct. 2016.')

parser.add_argument('-dx', type=str, default='1',
      help='''Horizontal grid resolution (default = 1).''')

parser.add_argument('-exp', type=str, default='M2_exp0',
      help='''Experiment name (default = M2_exp0).''')
parser.add_argument("--plot", dest="plot", action='store_true',
      help="Perform plotting at each time step.")
args = parser.parse_args()

def remap(cs,h0,e0,h,data0,z,depth):
   """
   remaps data0 into a regular grid
   """
   km,im = h0.shape
   data1 = numpy.ma.zeros((len(h),im))
   for i in range(im):
     h1 = numpy.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
        h1[h<-e0[0,i]] = 0.0 # top
        h1[h>depth[i]] = 0.0 # bottom
        # need to account for SSH and make sure total thicknesses are
        # the same
        dh = h0[:,i].sum() - h1.sum()
        tmp1 = numpy.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
          if dh > 0.0:
             # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = numpy.ma.masked

     # for debugging
        #if h0[:,i].sum() != h1.sum():
        #   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
        data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
        # mask
        data1[h1==0.0,i] = numpy.ma.masked;
     else: # land/iceshelf
        data1[:,i] = numpy.ma.masked

   return data1


path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp

x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]

nxOut = x.shape[1]
nyOut = x.shape[0]
nzOut = 70
dzOut = 4000./nzOut

nzExtra = 0

# the z location of grid-cell corners (or top-bottom inderfaces) on the output
# grid, with extra points at the top to catch anything above 0
zInterfaceOut = -dzOut*numpy.arange(-nzExtra, nzOut+1)

# the z location of grid-cell centers on the output grid
zOut = 0.5*(zInterfaceOut[0:-1] + zInterfaceOut[1:])

ncIn = Dataset(path+'/ocean_month.nc', 'r')

variables = ncIn.variables
nTime = len(ncIn.dimensions['time'])
nx = len(ncIn.dimensions['xh'])
ny = len(ncIn.dimensions['yh'])
nz = len(ncIn.dimensions['zl'])

assert(nx == nxOut)
assert(ny == nyOut)

# 1D x and y of rho points
x = 1e3*variables['xh'][:]
y = 1e3*variables['yh'][:]

dx = x[1]-x[0]
dy = y[1]-y[0]

print "Computing OSF for exp {}".format(path)
zInterface = variables['e'][-24:, :, :, :].mean(axis=0)
h = variables['h'][-24:, :, :, :].mean(axis=0)
# the last u point is on the boundary and can (should?) be ignored
v = variables['v'][-24:, :, 0:-1, :].mean(axis=0)
zInterface_v = 0.5*(zInterface[:, 0:-1, :] + zInterface[:, 1:,:])
h_v = 0.5*(h[:, 0:-1, :] + h[:, 1:, :])
vMask = numpy.ones((nz, ny-1, nx), bool)
osfOut = computeOSF(v, vMask, dx=dx*numpy.ones(v.shape),
                    zInterface=zInterface_v,
                    zInterfaceOut=zInterfaceOut, plot=False,
                    yPlot=y, zPlot=zInterface, tIndex=0)

# skip the extra points we added at the top, since those aren't part of
# standard ISOMIP+ output
osf = \
      osfOut[nzExtra:, :]

km, jm, im = zInterface.shape

# mean ocean properties
temp=Dataset(path+'/ocean_month.nc').variables['temp'][-24:,:,:,:].mean(axis=0).mean(axis=2)
salt=Dataset(path+'/ocean_month.nc').variables['salt'][-24:,:,:,:].mean(axis=0).mean(axis=2)
rho = eos.wright_eos(temp,salt,2.0e7)

# to contour ice shelf and ocean bottom
depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:,im/2]
ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
e0 = Dataset(path+'/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)
h0 = Dataset(path+'/prog.nc').variables['h'][0:2,:,:,im/2].mean(axis=0)
ncIn.close()

# grid
dy = y[1]-y[0]
y0 = y[0]-0.5*dy
(YOut, ZOut) = numpy.meshgrid(y0 + dy*(numpy.arange(ny+1)),
                                      zInterfaceOut,
                                      indexing='xy')
if args.plot:
   os.system('mkdir OSF')
   psi_vals = numpy.linspace(-0.2,0.2,30)
   # remapping params
   cs = mom_remapping.Remapping_Cs()
   cs.remapping_scheme = 2
   cs.degree = 2
   z = numpy.linspace(0,4000.,500)
   h1 = 0.5* (z[0:-1]+z[1::])
   # grid
   Y,Z = numpy.meshgrid(y,h1)
   #remap rho
   rho1=remap(cs,h0,e0,h1,rho,z,depth)
   # interpolate psi (can't remap it, it has been integrated already)
   psi = numpy.zeros((len(h1),jm))
   for j in range(jm):
      if h0[:,j].sum() > 0.001:
         tmp = -0.5*(e0[0:-1,j]+e0[1::,j])
         psi[:,j] = numpy.interp(h1, tmp, osf[:,j])
      else:
         psi[:,j] = -1.0e32

   # mask bad values
   psi = numpy.ma.masked_where(numpy.ma.getmask(rho1), psi)
   psi = numpy.ma.masked_where(psi==-1.0e32, psi)

   fig = plt.figure(figsize=(10,8),facecolor='white')
   ax = fig.add_subplot(111,axisbg='black')
   #ct=ax.contourf(1e-3*Y,-Z,psi*1.0e-6,psi_vals,cmap=plt.cm.bwr,extend='both')
   ct = ax.pcolormesh(1e-3*YOut, ZOut, 1e-6*osf,vmin=-0.2,vmax=0.2,cmap=plt.cm.bwr)
   cb = plt.colorbar(ct,orientation='horizontal',ticks=[-0.2,-0.1,0.0,0.1, 0.2],extend='both')
   cb.set_label('Overturning circulation [sv]', fontsize=18)
   s=ax.contour(1.0e-3*Y,-Z,rho1-1000,[37.0,37.1,37.2],colors='k',lw=0.5)
   ax.clabel(s, inline=1, fontsize=10,fmt='%4.2f')
   # fill IS and bottom
   ax.fill_between(1e-3*y, e0[0,:], 0.0, facecolor='white', interpolate=True)
   ax.fill_between(1e-3*y, e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
   ax.plot(1e-3*y, ssh,'k-',lw=1.5)
   ax.set_xlim([0,1000])
   ax.set_ylim([-4000,0])
   ax.set_xlabel('y [km]', fontsize=18)
   ax.set_ylabel('depth [m]', fontsize=18)
   plt.savefig('OSF/'+args.exp+'_dx'+args.dx+'_overturning_circulation.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()

print 'Saving streamfunction stats...'
os.system('mkdir TXT')
tmp1 = numpy.nonzero(y<=200.0e3)[0][-1]
osf = osf * 1.0e-6 # in sv
# cell averaged osf
h_mean = h.mean(axis=2)
h_mean = numpy.ma.masked_where(h_mean<0.001,h_mean)
asf_mean =  (np.abs(osf[:,0:tmp1])*h_mean[:,0:tmp1]).sum()/h_mean[:,0:tmp1].sum()
print 'Ice shelf ends at y =',y[tmp1]
text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_streamfunction.txt', "w")
text_file.write("%f \n" % np.abs(osf[:,0:tmp1]).max())
text_file.write("%f \n" % asf_mean)
text_file.close()

print 'Mean psi:', asf_mean
print 'Max psi:', np.abs(osf[:,0:tmp1]).max()
print 'Done!'
