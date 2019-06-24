from remapping import mom_remapping
import numpy
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import sys, os
sys.path.append('../../')
from misc import *
import wright_eos as eos
from computeOSF import computeOSF


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
temp=Dataset(path+'/prog.nc').variables['temp'][-24:,:,:,:].mean(axis=0).mean(axis=2)
salt=Dataset(path+'/prog.nc').variables['salt'][-24:,:,:,:].mean(axis=0).mean(axis=2)
rho = eos.wright_eos(temp,salt,2.0e7)

# to contour ice shelf and ocean bottom
depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:,im/2]
ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
e0 = Dataset(path+'/prog.nc').variables['e'][-24:,:,:,im/2].mean(axis=0)
h0 = Dataset(path+'/prog.nc').variables['h'][-24:,:,:,im/2].mean(axis=0)
ncIn.close()

# grid
dy = y[1]-y[0]
y0 = y[0]-0.5*dy
(YOut, ZOut) = numpy.meshgrid(y0 + dy*(numpy.arange(ny+1)),
                                      zInterfaceOut,
                                      indexing='xy')
