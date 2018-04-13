from matplotlib import pyplot as plt
import netCDF4
import numpy as np
from scipy.optimize import curve_fit
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 16})

#plt.style.use('ggplot')
c1 = '#6495ed'
c2 = '#ff6347'
c3 = 'k'
c4 = '#3cb371'
colors = ['#6495ed','#ff6347','k','#3cb371']
# plot some metrics for runs with varing wind forcing

path='ncfiles/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx2 = ['1 km','2 km','5 km','10 km']
dx = ['dx1','dx2','dx5','dx10']

param = 'heat_budget_dx_and_melting_OSM'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = np.array([-5,-2.5,0.0,2.5,5.0])
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2, ax3)) = plt.subplots(3, 1, figsize=(8,16))

for j in range(len(dx)):
  out1_mean = []; out1_std = []
  out2_mean = []; out2_std = []
  in1_mean = []; in1_std = []
  in2_mean = []; in2_std = []
  lat1_mean = []; lat2_mean = []
  eff1_mean = []
  eff2_mean = []
  for i in range(len(exps1)):
     print path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc'
     print path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc'
     # Hout 
     out1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hout'][:]))
     out1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hout'][:]))
     out2_mean.append(-np.mean(netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hout'][:]))
     out2_std.append(np.std(netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hout'][:]))
     # Hin
     in1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hin'][:]))
     in1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hin'][:]))
     # Latent
     lat1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Lat'][:]))
     lat2_mean.append(-np.mean(netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Lat'][:]))
     
     in2_mean.append(-np.mean(netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hin'][:]))
     in2_std.append(np.std(netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc').variables['Hin'][:]))

     # (Hin + Hout)/Hin
     eff1_mean.append((in1_mean[i]+out1_mean[i])/in1_mean[i])
     #eff1_mean.append(lat1_mean[i]/in1_mean[i])
     eff2_mean.append((in2_mean[i]+out2_mean[i])/in2_mean[i])
     #eff2_mean.append(lat2_mean[i]/in2_mean[i])
     #eff2_mean.append((np.abs(in2_mean[i])-out2_mean[i])/np.abs(in2_mean[i]))

  # plot
  ax1.errorbar(wind, in1_mean, in1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx2[j])
  ax2.errorbar(wind, out1_mean, out1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx2[j])
  ax3.plot(wind, eff1_mean, 'o-', color=colors[j], lw=2)
  ax3.plot([-10,10], [0,0], '-', color='gray', lw=1)

ax2.legend(loc='best', fontsize=14, ncol=2)
ax1.set_ylabel(r'$H_{in}$ [TW]', fontsize=20)
ax2.set_ylabel(r'$H_{out}$ [TW]', fontsize=20)
ax3.set_ylabel(r'$\delta_{th}$', fontsize=20)
#ax1.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
#ax2.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax3.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

#ax1.set_title('Melting on', fontsize=20)
#ax2.set_title('Melting off',fontsize=20)
ax1.set_xlim(-5.5,5.5)
ax2.set_xlim(-5.5,5.5)
ax3.set_xlim(-5.5,5.5)
ax1.set_ylim(-0.1,2.5)
ax2.set_ylim(-2.5,0.1)
ax3.set_ylim(-0.1,0.4)
ax1.set_xticks((wind))
ax2.set_xticks((wind))
ax3.set_xticks((wind))

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
