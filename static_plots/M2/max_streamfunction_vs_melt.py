from matplotlib import pyplot as plt
import netCDF4
import numpy as np
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

path1='TXT/'
path2='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx1 = ['1km','2km','5km','10km']
param = 'max_streamfunction_vs_melt'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1)) = plt.subplots(1, 1, figsize=(9,7))

for j in range(len(dx)):
  data1_mean = []
  data2_mean = []
  for i in range(len(exps1)):
     data1 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_streamfunction.txt')[0]
     data1_mean.append(data1)
     path_to_file1 = path2+'dx'+dx1[j]+'/'+exps1[i]+'/'+exps1[i]+'_'+dx1[j]+'.nc'
     print path_to_file1
     data2 = netCDF4.Dataset(path_to_file1).variables['Melt'][:]
     data2_mean.append(data2.mean())

  # plot
  ax1.plot(data2_mean, data1_mean, linestyle='', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j])


ax1.legend(loc='upper left', fontsize=14, ncol=2)
ax1.set_ylabel(r'Max. $|\overline{\psi}|$ [sv]', fontsize=20)
ax1.set_xlabel(r'Melt rate [m yr$^{-1}$]', fontsize=20)

ax1.set_xlim(4,10)
ax1.set_ylim(0.06,0.11)

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
