from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx2km/Sigma_zstar/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M2_exp0','M2_exp2']
fnames = ['M2_exp0_2km','M2_exp2_2km']
param = 'variable_wind'
pos_x = [-3.75,1.25]
labels_var = ['Cooling','Warming']
labels = ['-7.5','-5.0','-2.5','0','2.5','5.0','7.5']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Total Sea ice volume [m$^3$]','<$\dot{m}$>  [m yr$^{-1}$]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(exps))
   std_all = np.zeros(len(exps))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+exps[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       mean_all[i]=data.mean()
       std_all[i]=data.std()
       plt.errorbar(pos_x[i], mean_all[i], yerr=std_all[i], fmt='o',label=labels_var[i])

   #plt.plot(range(len(exps)), mean_all,'k-')
   plt.legend(loc=1)
   plt.xticks(pos_x[2:], labels, rotation='horizontal')
   #plt.xticks(pos_x[i], labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=18)
   plt.xlabel(r'$U_{asf}$ [m s$^{-1}$]', fontsize=16)
   plt.ylabel(r''+units[n], fontsize=16)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()

pos_x = [-3.75,1.25,-7.5,-5.0,-2.5,0,2.5,5.0,7.5]
exps = ['M2_exp0','M2_exp2','M2_exp17','M2_exp3','M2_exp5','M2_exp7','M2_exp6','M2_exp4','M2_exp18']
fnames = ['M2_exp0_2km','M2_exp2_2km','M2_exp17_2km','out1/M2_exp3_2km','M2_exp5_2km','M2_exp7_2km','M2_exp6_2km','out1/M2_exp4_2km','M2_exp18_2km']


for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(exps))
   std_all = np.zeros(len(exps))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+exps[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       mean_all[i]=data.mean()
       std_all[i]=data.std()
       if i==0:
          plt.errorbar(pos_x[i], mean_all[i], yerr=std_all[i], fmt='o',label=labels_var[i],color='blue',ms=8)

       if i==1:
          plt.errorbar(pos_x[i], mean_all[i], yerr=std_all[i], fmt='o',label=labels_var[i],color='red',ms=8)

   plt.errorbar(pos_x[2:], mean_all[2:], yerr=std_all[2:], fmt='o',color='black',ms=8)
   plt.plot(pos_x[2:], mean_all[2:],'k-')
   #plt.plot(range(len(exps)), mean_all,'k-')
   plt.legend(loc='upper left')
   plt.xticks(pos_x[2:], labels, rotation='horizontal')
   #plt.xticks(pos_x[i], labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=18)
   plt.xlabel(r'$U_{asf}$ [m s$^{-1}$]', fontsize=16)
   plt.ylabel(r''+units[n], fontsize=16)
   plt.savefig(param+'_sensitivity_all_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()



