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

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M1_exp5','M1_exp9']
fnames = ['out3/M1_exp5_1km','out3/M1_exp9_1km']

exps1 = ['M1_exp15','M1_exp16']
fnames1 = ['out3/M1_exp15_1km','out3/M1_exp16_1km']
param = 'trough'
labels = ['-5.0','5.0']
#vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt','AABW','CDW','CDW1','polynyaArea']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt','AABW','polynyaArea']
#units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Sea ice volume [m$^3$]','Melt rate  [m yr$^{-1}$]','AABW transport [Sv]','CDW transport [Sv]','CDW transport [Sv]','Coastal polynya area [km$^2$]']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Sea ice volume [m$^3$]','Melt rate  [m yr$^{-1}$]','AABW transport [Sv]','Coastal polynya area [km$^2$]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

width = 0.35
pos2 = [0,1,3,4]
for n in range(len(vars)):
   var = vars[n]

   fig, ax = plt.subplots(figsize=(10,8))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+fnames[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       rects1 = ax.bar(i, np.abs(data.mean()), width, color=color1, yerr=np.abs(data.std()),ecolor='k')

   for i in range(len(exps1)):
       print 'path',path+exps1[i]+'/'+fnames1[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+exps1[i]+'/'+fnames1[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+exps1[i]+'/'+fnames1[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       rects2 = ax.bar(pos2[i]+width, np.abs(data.mean()), width, color=color2, yerr=np.abs(data.std()),ecolor='k')

 #  plt.errorbar(xval, mean_all, yerr=std_all, fmt='o')
  # plt.plot(xval, mean_all,'k-')

   # add some text for labels, title and axes ticks
   ax.set_xlim(-0.1,len(exps))
   ax.set_xticks(np.arange(len(exps)) + width/2 )
   ax.set_xticklabels((labels))
   #ax.set_xticks(len(exps) + width / 2)
   #ax.set_xticklabels(exps)
   ax.legend((rects1[0], rects2[0]), ('with trough', 'without trough'),bbox_to_anchor=(0.8, 1.13),ncol=2,fontsize=16)
   ax.set_ylabel(r''+units[n], fontsize=20)
   ax.set_xlabel(r'Maximum U$_{asf}$ [m s$^{-1}$]', fontsize=20)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   #plt.close()
   plt.show()
