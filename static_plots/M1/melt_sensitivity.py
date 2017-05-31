from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 20})

#plt.style.use('ggplot')
color1 = '#6495ed'
color2 = '#ff6347'
color3 = '#8470ff'
color4 = '#3cb371'

# plot some metrics for runs with varing wind forcing

path='/work/gmm/Projects/Ideal_ice_shelf/Mode1/dx1km/Sigma_zstar/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M1_exp5','M1_exp13','M1_exp21','M1_exp22','M1_exp23','M1_exp24']
fnames = ['out3/M1_exp5_1km','out3/M1_exp13_1km','out3/M1_exp21_1km','out3/M1_exp22_1km','out3/M1_exp23_1kmOLD','out4/M1_exp24_1km']

param = 'melting'
# get fw fluxes
FW_mean = np.zeros(len(exps)); FW_std = np.zeros(len(exps));
for i in range(len(exps)):
    print 'path',path+exps[i]+'/'+fnames[i]+'.nc'
    FW_mean[i] = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables['FwIS'][:].mean()/1.0e6
    FW_std[i] = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables['FwIS'][:].std()/1.0e6

vars1 = netCDF4.Dataset(path+exps[0]+'/'+fnames[0]+'.nc').variables
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

width = 0.35
ignore = ['time','x','y']
for var in vars1:
  print 'variable:',var
  if var not in ignore:
   fig, ax = plt.subplots(figsize=(10,8))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+fnames[i]+'.nc'
       data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]
       units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       try:
         description = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].description      
       except:
         description = var

       # plot
       #rects1 = ax.bar(FW_mean[i], data.mean(), width, color=color1, yerr=data.std(),ecolor='k')
       ax.errorbar(FW_mean[i], data.mean(), xerr=FW_std[i], yerr=data.std(),ecolor='b',color='b',marker='o')

   # add some text for labels, title and axes ticks
   ax.set_xlim(-0.1,FW_mean.max()+0.1)
   ax.set_title(description,fontsize=14)
   ax.set_ylabel(r''+var+' ['+units+']', fontsize=20)
   ax.set_xlabel(r'Sub-ice-shelf fresh water flux x 1.0E6 [kg s$^{-1}$]', fontsize=20)
   plt.grid()
   plt.savefig(param+'_sensitivity_'+var+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()
   #plt.show()

vars2 = netCDF4.Dataset(path+exps[0]+'/'+fnames[0]+'_buoyancy.nc').variables

for var in vars2:
  print 'variable:',var
  if var not in ignore:
   fig, ax = plt.subplots(figsize=(10,8))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+fnames[i]+'_buoyancy.nc'
       data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'_buoyancy.nc').variables[var][:]
       units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'_buoyancy.nc').variables[var].units
       # plot
       #rects1 = ax.bar(FW_mean[i], data.mean(), width, color=color1, yerr=data.std(),ecolor='k')
       ax.errorbar(FW_mean[i], data.mean(), xerr=FW_std[i], yerr=data.std(),ecolor='b',color='b',marker='o')

   # add some text for labels, title and axes ticks
   ax.set_xlim(-0.1,FW_mean.max()+0.1)
   ax.set_title(var,fontsize=20)
   ax.set_ylabel(r''+var+' ['+units+']', fontsize=20)
   ax.set_xlabel(r'Sub-ice-shelf fresh water flux x 1.0E6 [kg s$^{-1}$]', fontsize=20)
   plt.grid()
   plt.savefig(param+'_sensitivity_'+var+'_buoyancy.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()
   #plt.show()
