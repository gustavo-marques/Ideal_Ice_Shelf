#!/usr/bin/env python

# Gustavo Marques

#import matplotlib
#matplotlib.use('Agg')
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import os

ssh = Dataset('00010101.ocean_sfc.nc').variables['SSH'][:]
nt, ny, nx = ssh.shape
ssh_tend = np.zeros((nt-2,ny,nx))
for t in range(1,nt-1):
   ssh_tend[t-1] = ssh[t+1,:,:] - ssh[t-1,:,:]

ssh_tend = np.ma.masked_where(np.abs(ssh_tend)>10.0, ssh_tend)

name = 'SSH_tend'
ncfile = Dataset(name+'.nc','w', format = 'NETCDF3_CLASSIC')
ncfile.createDimension('nx',nx)
ncfile.createDimension('ny',ny)
ncfile.createDimension('time',nt-2)
s = ncfile.createVariable('SSH_Tend','double',('time','ny','nx'), fill_value = 1.e+20)
s.units = 'm'
s.standard_name =  'ice shelf thickness'
s.missing_value = 1.e+20
s[:] = ssh_tend[:]
ncfile.sync()
ncfile.close()

