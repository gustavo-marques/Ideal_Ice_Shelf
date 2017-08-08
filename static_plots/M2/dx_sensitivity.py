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

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M2_exp0','M2_exp4','M2_exp13','M2_exp14']
#fnames = ['dx1km/Sigma_zstar/M1_exp5/out3/M1_exp5_1km','dx2.5km/Sigma_zstar/M1_exp5/out2/M1_exp5_2.5km','dx5km/Sigma_zstar/M1_exp5/M1_exp5_5km','dx10km/Sigma_zstar/M1_exp5/out1/M1_exp5_10km']
dx = ['5km','10km']

param = 'dx_and_melting'
labels = ['5','2.5','5','10']
vars1 = netCDF4.Dataset(path+'dx'+dx[0]+'/Sigma_zstar/'+exps[0]+'/'+exps[0]+'_'+dx[0]+'.nc').variables
vars2 = netCDF4.Dataset(path+'dx'+dx[0]+'/Sigma_zstar/'+exps[0]+'/'+exps[0]+'_'+dx[0]+'_buoyancy.nc').variables
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

ignore = ['time','x','y','B0_iceshelf_mean','B0','B0_shelf_mean','B0_mean','Monin_Obukhov_lenght','maxSeaiceThick','seaiceArea','polynyaArea','totalSaltFlux','totalFw','FwIS']
for var in vars2:
  print 'variable:',var
  if var not in ignore:
   fig, ax = plt.subplots(figsize=(10,8))
   width = 1
   for j in range(len(dx)):
     for i in range(len(exps)):
       path_to_file = path+'dx'+dx[j]+'/Sigma_zstar/'+exps[i]+'/'+exps[i]+'_'+dx[j]+'_buoyancy.nc'
       print 'path',path_to_file
       data = netCDF4.Dataset(path_to_file).variables[var][:]
       units = netCDF4.Dataset(path_to_file).variables[var].units
       try:
         description = netCDF4.Dataset(path_to_file).variables[var].description
       except:
         description = var

       # plot
       ax.errorbar(width, data.mean(), yerr=data.std(),ecolor='b',color='b',marker='o')
       width = width+1
     # add some text for labels, title and axes ticks
   #ax.set_xlim(-0.1,len(dx)*len(exps)+1)
   #ax.set_xticklabels((dx))
   ax.set_title(description,fontsize=14)
   ax.set_ylabel(r''+var+' ['+units+']', fontsize=20)
   #ax.set_xlabel(r'Sub-ice-shelf fresh water flux x 1.0E6 [kg s$^{-1}$]', fontsize=20)
   plt.grid()
   plt.savefig(param+'_sensitivity_'+var+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()
