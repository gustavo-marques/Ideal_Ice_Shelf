from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/'
dx_path = ['dx0.5km/','dx1km/','dx2km/','dx5km/']
fnames = ['out3/M2_exp2_500m','out1/M2_exp2_1km','M2_exp2_2km','M2_exp2_5km']
param = 'dx'
labels = ['0.5','1.0','2.0','5.0']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Total Sea ice volume [m$^3$]','<$\dot{m}$>  [m yr$^{-1}$]']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(fnames))
   std_all = np.zeros(len(fnames))
   for i in range(len(fnames)):
       print 'path',path+dx_path[i]+'Sigma_zstar/M2_exp2/'+fnames[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+dx_path[i]+'Sigma_zstar/M2_exp2/'+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+dx_path[i]+'Sigma_zstar/M2_exp2/'+fnames[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+fnames[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       mean_all[i]=data.mean()
       std_all[i]=data.std()

   plt.errorbar(range(len(fnames)), mean_all, yerr=std_all, fmt='o',color='black',ms=8)
   plt.plot(range(len(fnames)), mean_all,'k-')

   plt.xticks(range(len(fnames)), labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=18)
   plt.xlabel(r'$\Delta x$ [km]', fontsize=16)
   plt.ylabel(r''+units[n], fontsize=16)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()
