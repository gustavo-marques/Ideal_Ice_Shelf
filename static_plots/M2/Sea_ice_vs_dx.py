from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 16})

#plt.style.use('ggplot')
color1 = '#6495ed'
color2 = '#ff6347'
color3 = '#8470ff'
color4 = '#3cb371'
colors = ['#6495ed','#ff6347','k','#3cb371']

# plot some metrics for runs with varing wind forcing

path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['1km','2km','5km','10km']

param = 'seaice_dx_and_melting'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2)) = plt.subplots(1, 2, sharey='row', figsize=(12,7))

for j in range(len(dx)):
  ice1_mean = []; ice1_std = []
  ice2_mean = []; ice2_std = []
  for i in range(len(exps1)):
     path_to_file1 = path+'dx'+dx[j]+'/'+exps1[i]+'/'+exps1[i]+'_'+dx[j]+'_buoyancy.nc'
     path_to_file2 = path+'dx'+dx[j]+'/'+exps2[i]+'/'+exps2[i]+'_'+dx[j]+'_buoyancy.nc'
     print 'Path1, path2',path_to_file1, path_to_file2
     time1 = netCDF4.Dataset(path_to_file1).variables['time'][:]
     time2 = netCDF4.Dataset(path_to_file2).variables['time'][:]
     print 'Melt on, path/time:',path_to_file1, time1
     print 'Melt off path/time:',path_to_file2, time2
     print '--------------------------------- \n'
     data1 = netCDF4.Dataset(path_to_file1).variables['seaiceVolume'][:]/1.0e12 # 1000 km^3
     data2 = netCDF4.Dataset(path_to_file2).variables['seaiceVolume'][:]/1.0e12
     ice1_mean.append(data1.mean())
     ice2_mean.append(data2.mean())
     ice1_std.append(data1.std())
     ice2_std.append(data2.std())

  # plot
  ax1.errorbar(wind, ice1_mean, ice1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx[j])
  ax2.errorbar(wind, ice2_mean, ice2_std, linestyle='-', marker='o', color=colors[j], elinewidth=2)

ax1.legend(loc='lower left', fontsize=14, ncol=2)
ax1.set_ylabel(r'Sea ice volume [1000 km$^3$]', fontsize=20)
ax1.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax2.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax1.set_title('a) Melting on', fontsize=20)
ax2.set_title('b) Melting off',fontsize=20)
ax1.set_xlim(-5.5,5.5)
ax2.set_xlim(-5.5,5.5)
ax1.set_ylim(0,1.6)
ax2.set_ylim(0,1.6)
ax1.set_xticks((wind))
ax2.set_xticks((wind))
plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
