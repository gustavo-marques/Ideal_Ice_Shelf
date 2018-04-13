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
dx1 = ['1','2','5','10']

param = 'mean_overturning_transport'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2)) = plt.subplots(1, 2, sharey='row', figsize=(12,7))

for j in range(len(dx)):
  data1_mean = []
  data2_mean = []
  for i in range(len(exps1)):
     data1 = np.loadtxt('TXT/'+exps1[i]+'_'+dx[j]+'_streamfunction.txt')[1]
     data2 = np.loadtxt('TXT/'+exps2[i]+'_'+dx[j]+'_streamfunction.txt')[1]
     data1_mean.append(data1)
     data2_mean.append(data2)

  # plot
  ax1.plot(wind, data1_mean, linestyle='-', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j] +' km')
  ax2.plot(wind, data2_mean, linestyle='-', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j] +' km')


ax1.legend(loc='lower center', fontsize=14, ncol=2)
ax1.set_ylabel(r'Mean overturning transport [Sv]', fontsize=20)
ax1.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax2.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

ax1.set_title('a) Melting on', fontsize=20)
ax2.set_title('b) Melting off',fontsize=20)
ax1.set_xlim(-5.5,5.5)
ax2.set_xlim(-5.5,5.5)
ax1.set_ylim(0,0.014)
ax2.set_ylim(0,0.014)
ax1.set_xticks((wind))
ax2.set_xticks((wind))
  #plt.grid()

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
