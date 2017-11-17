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
colors = ['#6495ed','#ff6347','k','#3cb371','y']
# plot some metrics for runs with varing wind forcing

path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['1km','2km','5km','10km']
dx1 =[1.0,2.0,5.0,10.0]
param = 'Rd_dx_and_melting'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(12,12))

for j in range(len(dx)):
  for i in range(len(exps1)):
     path_to_file1 = path+'dx'+dx[j]+'/'+exps1[i]+'/prog.nc'
     path_to_file2 = path+'dx'+dx[j]+'/'+exps2[i]+'/prog.nc'
     print 'Melt on, path/time:',path_to_file1
     print 'Melt off path/time:',path_to_file2
     print '--------------------------------- \n'
     data1 = netCDF4.Dataset(path_to_file1).variables['Rd_dx'][-144::,:].mean(axis=0).mean(axis=1)*dx1[j] # ave last 2 years
     data2 = netCDF4.Dataset(path_to_file2).variables['Rd_dx'][-144::,:].mean(axis=0).mean(axis=1)*dx1[j]
     y = netCDF4.Dataset(path_to_file2).variables['yh'][:]

     # plot
     if j == 0:
       ax = ax1
     elif j == 1:
       ax = ax2
     elif j == 2:
       ax = ax3
     else:
       ax = ax4

     ax.plot(y, data1, '-', color=colors[i], label= labels[i] + 'm/s', lw=2)
     ax.plot(y, data2, '--', color=colors[i], lw=2)
     ax.plot([205,205],[-10,15],'--',color='gray')
     ax.plot([450,450],[-10,15],'--',color='gray')
     ax.plot([650,650],[-10,15],'--',color='gray')
  if j == 0:
     ax.legend(loc='upper left', fontsize=14, ncol=2)
  if j == 2:
     ax.set_ylabel(r'$R_d$ [km]', fontsize=20)
     ax.set_xlabel('y [km]', fontsize=20)

  ax.set_title(abcd[j] + r'$\Delta$x = ' + dx[j])
  ax.set_xlim(0,1000)
  ax.set_ylim(0,11)

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
