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

path='TXT/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx2 = ['1 km','2 km','5 km','10 km']
dx = ['dx1','dx2','dx5','dx10']

param = 'Rd_all'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = np.array([-5,-2.5,0.0,2.5,5.0])
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharey='row', figsize=(12,12))

for j in range(len(dx)):
  rd1_mean = []; rd1_std = []
  rd2_mean = []; rd2_std = []
  rd3_mean = []; rd3_std = []
  for i in range(len(exps1)):
     # Hout 
     rd1_mean.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_Rd.txt')[0,0])
     rd1_std.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_Rd.txt')[0,1])
     rd2_mean.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_Rd.txt')[1,0])
     rd2_std.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_Rd.txt')[1,1])
     rd3_mean.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_Rd.txt')[2,0])
     rd3_std.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_Rd.txt')[2,1])
      

  # plot
  ax1.errorbar(wind, rd1_mean, rd1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx2[j])
  ax2.errorbar(wind, rd2_mean, rd2_std, linestyle='-', marker='o', color=colors[j], elinewidth=2)
  ax3.errorbar(wind, rd3_mean, rd3_std, linestyle='-', marker='o', color=colors[j], elinewidth=2)


ax1.legend(loc='upper left', fontsize=14, ncol=2)
ax3.set_ylabel(r'Rd [km]', fontsize=20)
ax3.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

ax1.set_xlim(-5.5,5.5)
ax2.set_xlim(-5.5,5.5)
ax3.set_xlim(-5.5,5.5)
ax1.set_xticks((wind))
ax2.set_xticks((wind))
ax3.set_xticks((wind))
  #plt.grid()

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
