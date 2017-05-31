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

path='/work/gmm/Projects/Ideal_ice_shelf/Mode1/'
icfile = 'IDEAL_IS_IC.nc'
fnames = ['dx1km/Sigma_zstar/M1_exp5/out3/M1_exp5_1km','dx2.5km/Sigma_zstar/M1_exp5/out2/M1_exp5_2.5km','dx5km/Sigma_zstar/M1_exp5/M1_exp5_5km','dx10km/Sigma_zstar/M1_exp5/out1/M1_exp5_10km']
fnames1 = ['dx1km/Sigma_zstar/M1_exp13/out3/M1_exp13_1km','dx2.5km/Sigma_zstar/M1_exp13/out2/M1_exp13_2.5km','dx5km/Sigma_zstar/M1_exp13/M1_exp13_5km','dx10km/Sigma_zstar/M1_exp13/out1/M1_exp13_10km']

param = 'dx_and_melting'
labels = ['1','2.5','5','10']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt','AABW','CDW','CDW1','polynyaArea']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Sea ice volume [m$^3$]','Melt rate  [m yr$^{-1}$]','AABW transport [Sv]','CDW transport [Sv]','CDW transport [Sv]','Coastal polynya area [km$^2$]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

width = 0.35
pos2 = range(len(labels))
for n in range(len(vars)):
   var = vars[n]

   fig, ax = plt.subplots(figsize=(10,8))
   for i in range(len(labels)):
       print 'path',path+'/'+fnames[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+'/'+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+'/'+fnames[i]+'.nc').variables[var][:]
       # plot
       rects1 = ax.bar(i, np.abs(data.mean()), width, color=color1, yerr=np.abs(data.std()),ecolor='k')

   for i in range(len(labels)):
       print 'path',path+'/'+fnames1[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+'/'+fnames1[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+'/'+fnames1[i]+'.nc').variables[var][:]
       # plot
       rects2 = ax.bar(pos2[i]+width, np.abs(data.mean()), width, color=color2, yerr=np.abs(data.std()),ecolor='k')

 #  plt.errorbar(xval, mean_all, yerr=std_all, fmt='o')
  # plt.plot(xval, mean_all,'k-')

   # add some text for labels, title and axes ticks
   ax.set_xlim(-0.1,len(labels))
   ax.set_xticks(np.arange(len(labels)) + width/2 )
   ax.set_xticklabels((labels))
   ax.legend((rects1[0], rects2[0]), ('melting on', 'melting off'),bbox_to_anchor=(0.8, 1.13),ncol=2,fontsize=16)
   ax.set_ylabel(r''+units[n], fontsize=20)
   ax.set_xlabel(r'$\Delta x$ [km]', fontsize=20)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.close()
   #plt.show()
