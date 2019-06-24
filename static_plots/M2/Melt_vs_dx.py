from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 16})

#plt.style.use('ggplot')
color1 = '#6495ed'
color2 = '#ff6347'
color3 = 'k'
color4 = '#3cb371'

# plot some metrics for runs with varing wind forcing

path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['1km','2km','5km','10km']
dx = ['2km','5km','10km']

param = 'Melt_dx'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ax = plt.subplots(1, 1, figsize=(9,8))

for j in range(len(dx)):
  OHT1_mean = []; OHT1_std = []
  OHT2_mean = []; OHT2_std = []
  for i in range(len(exps1)):
     path_to_file1 = path+'dx'+dx[j]+'/'+exps1[i]+'/'+exps1[i]+'_'+dx[j]+'.nc'
     time1 = netCDF4.Dataset(path_to_file1).variables['time'][:]
     print 'Melt on, path/time:',path_to_file1, time1
     print '--------------------------------- \n'
     data1 = netCDF4.Dataset(path_to_file1).variables['Melt'][:]
     OHT1_mean.append(data1.mean())
     OHT1_std.append(data1.std())
     #OHT2_std.append(data2.std())

  if j == 0:
     c = color1
  elif j == 1:
     c = color2
  elif j == 2:
     c = color3
  else:
     c = color4

  ax.errorbar(wind, OHT1_mean, OHT1_std, linestyle='-', marker='o', color=c, elinewidth=3, label=r'$\Delta$x = '+dx[j])

ax.legend(loc='upper left', fontsize=14, ncol=2)
ax.set_ylabel(r'Mean melt rate [m yr$^{-1}$]', fontsize=20)
ax.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax.set_xlim(-5.5,5.5)
#ax.set_ylim(4,10)
ax.set_xticks((wind))
plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
