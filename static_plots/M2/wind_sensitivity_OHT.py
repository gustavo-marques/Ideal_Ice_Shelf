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
dx='2km'
path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx'+dx+'/Sigma_zstar'
icfile = 'IDEAL_IS_IC.nc'
fnames1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4']
fnames2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14']

param = 'oht1'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
labels1 = [-5.0,-2.5,0.0,2.5,5.0]
labels2 = [-5.0,-2.5,0.0,2.5,5.0]
vars = ['OHT_ice_shelf','OHT_shelf']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

width = 0.5
for n in range(len(vars)):
   var = vars[n]

   fig, ax = plt.subplots(figsize=(10,8))
   # Set the ticks and ticklabels for all axes
   plt.setp(ax, xticks=labels1, xticklabels=labels)
   for i in range(len(labels1)):
       path_file = path+'/'+fnames1[i]+'/'+fnames1[i]+'_'+dx+'_oht.nc'

       print 'path',path_file
       data = netCDF4.Dataset(path_file).variables[var][:]/1.0e12 # TW
       # plot
       rects1 = ax.bar(labels1[i]-width, np.abs(data.mean()), width, color=color1, yerr=np.abs(data.std()),ecolor='k')

   for i in range(len(labels2)):
       path_file = path+'/'+fnames2[i]+'/'+fnames2[i]+'_'+dx+'_oht.nc'
       print 'path',path_file
       data = netCDF4.Dataset(path_file).variables[var][:]/1.0e12 # TW
       # plot
       #rects2 = ax.bar(labels2[i], np.abs(data.mean()), width, color=color2, yerr=np.abs(data.std()),ecolor='k')

 #  plt.errorbar(xval, mean_all, yerr=std_all, fmt='o')
  # plt.plot(xval, mean_all,'k-')

   # add some text for labels, title and axes ticks
   ax.set_xlim(-6.5,6.5)
   #ax.set_xticks(np.arange(len(labels)) + width/2 )
   #ax.set_xticklabels((labels))
   #ax.legend((rects1[0], rects2[0]), ('melting on', 'melting off'),bbox_to_anchor=(0.8, 1.13),ncol=2,fontsize=16)
   ax.set_ylabel(r''+units[n], fontsize=20)
   ax.set_xlabel(r'$U_{asf}$ [m s$^{-1}$]', fontsize=20)
   #plt.grid()
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   #plt.close()
   plt.show()
