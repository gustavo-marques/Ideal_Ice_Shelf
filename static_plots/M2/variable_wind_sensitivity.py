from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx2km/Sigma_zstar/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4','M2_exp13','M2_exp14']
fnames = ['M2_exp0_2km','M2_exp1_2km','M2_exp2_2km','M2_exp3_2km','M2_exp4_2km','M2_exp13_2km','M2_exp14_2km']
param = 'variable_wind'
pos_x = [-5,-2.5,0,2.5,5,-5,5]
#labels_var = ['Cooling','Warming']
labels_var = ['-5.0','-2.5','0','2.5','5.0','-5 nomelt','5 nomelt']
vars = ['OHT_ice_shelf','OHT_shelf','CDW','CDW1','Melt']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Sv','Sv','<$\dot{m}$>  [m yr$^{-1}$]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(exps))
   std_all = np.zeros(len(exps))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+fnames[i]+'.nc'
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
   plt.legend(loc=1,ncol=2)
   plt.xlim(-6,6)
   #plt.xticks(pos_x[2:], labels, rotation='horizontal')
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
   #plt.show()

vars = ['hml_shelf','Nshelf','Tshelf','Sshelf','seaiceVolume']
units = ['[m]','[s^-1]','[oC]','[psu]','[m^3]']

for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(exps))
   std_all = np.zeros(len(exps))
   for i in range(len(exps)):
       print 'path',path+exps[i]+'/'+fnames[i]+'_buoyancy.nc'
       data = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'_buoyancy.nc').variables[var][:]
       #units = netCDF4.Dataset(path+exps[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       mean_all[i]=data.mean()
       std_all[i]=data.std()
       plt.errorbar(pos_x[i], mean_all[i], yerr=std_all[i], fmt='o',label=labels_var[i])


   #plt.plot(range(len(exps)), mean_all,'k-')
   plt.legend(loc='upper left',ncol=3)
   plt.xlim(-6,6)
   #plt.xticks(pos_x[i], labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   #plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=18)
   plt.xlabel(r'$U_{asf}$ [m s$^{-1}$]', fontsize=16)
   plt.ylabel(r''+units[n], fontsize=16)
   plt.savefig(param+'_sensitivity_all_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()



