from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 16})

#plt.style.use('ggplot')
color1 = '#6495ed'
color2 = '#ff6347'
color3 = '#7063ed'
color4 = '#3cb371'

c=[color1,color2,color3,color4]
# plot some metrics for runs with varing wind forcing

path='ncfiles/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

#dx= ['1km','2km','5km','10km']
dx= ['dx2','dx5','dx10']
dx1= ['2 km','5 km','10 km']

param= 'snow_dx_and_melting'
labels= ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(12,10))

for j in range(len(dx)):
  ice1_mean = []; ice1_std = []
  ice2_mean = []; ice2_std = []
  for i in range(len(exps1)):
     path_to_file1 = path+exps1[i]+'_'+dx[j]+'_surface_heat_terms.nc'
     path_to_file2 = path+exps2[i]+'_'+dx[j]+'_surface_heat_terms.nc'
     print 'Path1, path2',path_to_file1, path_to_file2
     time1 = netCDF4.Dataset(path_to_file1).variables['time'][:]
     time2 = netCDF4.Dataset(path_to_file2).variables['time'][:]
     tmp = np.nonzero((time1>=15) & (time1<20.))[0]
     print 'Melt on, path/time:',path_to_file1, time1[tmp]
     print 'Melt off path/time:',path_to_file2, time2[tmp]
     print '--------------------------------- \n'
     data1 = netCDF4.Dataset(path_to_file1).variables['SV'][tmp]/1.0e12 # 1000 km^3
     data2 = netCDF4.Dataset(path_to_file2).variables['SV'][tmp]/1.0e12
     ice1_mean.append(data1.mean())
     ice2_mean.append(data2.mean())
     ice1_std.append(data1.std())
     ice2_std.append(data2.std())

  ice1_mean = np.array(ice1_mean)
  ice2_mean = np.array(ice2_mean)
  # plot
  if j == 0:
     ax = ax1
  elif j == 1:
     ax = ax2
  elif j == 2:
     ax = ax3
  else:
     ax = ax4

  ax.errorbar(wind, ice1_mean, ice1_std, linestyle='-', marker='o', color=color1, elinewidth=2, label='Melting on')
  ax.errorbar(wind, ice2_mean, ice2_std, linestyle='-', marker='o', color=color2, elinewidth=2, label='Melting off')

  if j == 0:
    ax.legend(loc='upper left', fontsize=14, ncol=2)
  if j == 2:
    ax.set_ylabel(r'Snow volume [1000 km$^3$]', fontsize=20)
    ax.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

  ax.set_title(abcd[j] + r'$\Delta$x = ' + dx1[j])
  ax.set_xlim(-5.5,5.5)
  #ax.set_ylim(0.,2.2)
  ax.set_xticks((wind))


plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
