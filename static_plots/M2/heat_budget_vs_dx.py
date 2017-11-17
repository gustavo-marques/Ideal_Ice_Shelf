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

dx = ['1km','2km','5km','10km']
dx = ['dx1','dx2','dx5','dx10']

param = 'heat_budget_dx_and_melting'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = np.array([-5,-2.5,0.0,2.5,5.0])
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharey='row', figsize=(12,12))

for j in range(len(dx)):
  out1_mean = []; out1_std = []
  out2_mean = []; out2_std = []
  in1_mean = []; in1_std = []
  in2_mean = []; in2_std = []
  eff1_mean = []
  eff2_mean = []
  for i in range(len(exps1)):
     # Hout 
     out1_mean.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_heat_budget.txt')[0])
     out1_std.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_heat_budget.txt')[1])
     out2_mean.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_heat_budget.txt')[0])
     out2_std.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_heat_budget.txt')[1])
     # Hin
     in1_mean.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_heat_budget.txt')[2])
     in1_std.append(np.loadtxt(path+exps1[i]+'_'+dx[j]+'_heat_budget.txt')[3])
     in2_mean.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_heat_budget.txt')[2])
     in2_std.append(np.loadtxt(path+exps2[i]+'_'+dx[j]+'_heat_budget.txt')[3])

     # (Hin + Hout)/Hin
     eff1_mean.append((np.abs(in1_mean[i])-out1_mean[i])/np.abs(in1_mean[i]))
     eff2_mean.append((np.abs(in2_mean[i])-out2_mean[i])/np.abs(in2_mean[i]))

  # plot
  ax1.errorbar(wind, in1_mean, in1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx[j])
  ax2.errorbar(wind, in2_mean, in2_std, linestyle='-', marker='o', color=colors[j], elinewidth=2)
  ax3.errorbar(wind, out1_mean, out1_std, linestyle='-', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx[j])
  ax4.errorbar(wind, out2_mean, out2_std, linestyle='-', marker='o', color=colors[j], elinewidth=2)

  ax5.plot(wind, eff1_mean, '-', color=colors[j], lw=2)
  ax6.plot(wind, eff2_mean, '-', color=colors[j], lw=2)

ax3.legend(loc='upper left', fontsize=14, ncol=2)
ax1.set_ylabel(r'$H_{in}$ [TW]', fontsize=20)
ax3.set_ylabel(r'$H_{out}$ [TW]', fontsize=20)
ax5.set_ylabel(r'$\delta_{th}$', fontsize=20)
ax5.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)
ax6.set_xlabel(r'U$_{shelf}$ [m s$^{-1}$]', fontsize=20)

ax1.set_title('Melting on', fontsize=20)
ax2.set_title('Melting off',fontsize=20)
ax1.set_xlim(-5.5,5.5)
ax2.set_xlim(-5.5,5.5)
ax3.set_xlim(-5.5,5.5)
ax4.set_xlim(-5.5,5.5)
ax5.set_xlim(-5.5,5.5)
ax6.set_xlim(-5.5,5.5)
ax1.set_ylim(-7,0)
ax2.set_ylim(-7,0)
ax3.set_ylim(0,7)
ax4.set_ylim(0,7)
ax1.set_xticks((wind))
ax2.set_xticks((wind))
ax3.set_xticks((wind))
ax4.set_xticks((wind))
ax5.set_xticks((wind))
ax6.set_xticks((wind))
  #plt.grid()

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
