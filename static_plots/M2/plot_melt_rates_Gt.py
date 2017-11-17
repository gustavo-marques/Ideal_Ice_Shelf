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

path='TXT/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['dx1','dx2','dx5','dx10']

param = 'melt_rate_Gt'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1)) = plt.subplots(1, 1, figsize=(12,7))

for j in range(len(dx)):
  data1_mean = []
  data1_std = []
  for i in range(len(exps1)):
     data1 = np.loadtxt('TXT/'+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0]
     data1_mean.append(data1)
     data1 = np.loadtxt('TXT/'+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[1]
     data1_std.append(data1)

  # plot
  ax1.errorbar(wind, data1_mean, data1_std, linestyle='-', marker='o', color=colors[j], elinewidth=3, label=r'$\Delta$x = '+dx[j])

ax1.legend(loc='upper left', fontsize=14, ncol=2)
ax1.set_ylabel(r'Melt rate [Gt yr$^{-1}$]', fontsize=20)
ax1.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

ax1.set_xlim(-5.5,5.5)
#ax1.set_ylim(0,0.7)
ax1.set_xticks((wind))
plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
