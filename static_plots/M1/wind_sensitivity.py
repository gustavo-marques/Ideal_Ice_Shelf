from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 20})

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/'
icfile = 'IDEAL_IS_IC.nc'
exps = ['M1_exp3','M1_exp5','M1_exp6','M1_exp7','M1_exp8','M1_exp4']
fnames = ['out4/M1_exp3_1km','out3/M1_exp5_1km','out3/M1_exp6_1km','out3/M1_exp7_1km','out3_tmp/M1_exp8_1km','out3_tmp/M1_exp4_1km']
param = 'wind'
xval = [-7.5,-5,-2.5,0,2.5,7.5]
labels = ['-7.5','-5','-2.5','0.0','2.5','7.5']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt','AABW','CDW','polynyaArea']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Sea ice volume [m$^3$]','Melt rate  [m yr$^{-1}$]','AABW transport [sv]','CDW transport [sv]','Coastal polynya area [km$^2$]']
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

   plt.errorbar(xval, mean_all, yerr=std_all, fmt='o')
   plt.plot(xval, mean_all,'k-')

   #plt.xticks(range(len(exps)), labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=20)
   plt.xlabel(r'$U_{asf}$ [m s$^{-1}$]', fontsize=20)
   plt.ylabel(r''+units[n], fontsize=20)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()
