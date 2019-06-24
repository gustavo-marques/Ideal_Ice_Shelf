from matplotlib import pyplot as plt
import netCDF4
import numpy as np
from scipy import stats
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

#plot = False
plot = True
path='TXT/'
path1='ncfiles/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp0_NT','M2_exp1_NT','M2_exp2_NT','M2_exp3_NT','M2_exp4_NT'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx1 = ['1km','2km','5km','10km']

dx = ['dx5']
dx1 = ['5km']

param = 'cavity_transport_vs_Hin'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd)']
tx=[1.4,1.7,2.0,2.3]
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

#f, ((ax1),(ax2)) = plt.subplots(1, 2, figsize=(14,8))
f, ((ax1)) = plt.subplots(1, 1, figsize=(10,8))

for j in range(len(dx)):
  tr1_mean = []
  tr2_mean = []
  Hi1_mean = []
  Hi2_mean = []
  print 'dx = '+dx[j]
  for i in range(len(exps1)):
     # melting on
     tr1 = np.loadtxt(path+exps1[i]+'_'+dx[j]+'_cavity_transports.txt')[2]
     tr1_mean.append(np.abs(tr1)/1.0e6)
     Hi1 = np.loadtxt(path+exps1[i]+'_'+dx[j]+'_heat_budget_monthly.txt')[2] 
     Hi1_mean.append(np.abs(Hi1))
     # melting off
     tr2 = np.loadtxt(path+exps2[i]+'_'+dx[j]+'_cavity_transports.txt')[2]
     tr2_mean.append(np.abs(tr2)/1.0e6)
     Hi2 = np.loadtxt(path+exps2[i]+'_'+dx[j]+'_heat_budget_monthly.txt')[2]         
     Hi2_mean.append(np.abs(Hi2))


  # plot
  ax1.plot(tr1_mean, Hi1_mean, linestyle='', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j], markersize=10)
  ax1.plot(tr2_mean, Hi2_mean, linestyle='', marker='o', color=colors[j+1], lw=2, label=r'$\Delta x$ = '+dx1[j], markersize=10)
  print 
  # fit line 
  yy=np.array(Hi1_mean)
  xx=np.array(tr1_mean)
  #m,b = np.polyfit(xx,yy, 1)
  m, b, r_value, p_value, std_err = stats.linregress(xx, yy)
  print m, b, r_value,  p_value, std_err
  transp = np.linspace(np.min(xx),np.max(xx),10)
  y=m*transp+b
  m = format(m, '.2f')
  r_value = format(r_value, '.2f')
  if plot:
    ax1.plot(transp,y,linestyle='-', color=colors[j], lw=2)
    # TEXT
    ax1.text(0.025,tx[j],r'$TF$ = '+str(m)+' (r = '+str(r_value)+')', color=colors[j], fontsize=20)


ax1.legend(loc='upper left', fontsize=18, ncol=4)
ax1.set_ylabel(r'$H_{in}$ [TW]', fontsize=20)
ax1.set_xlabel(r'$\Phi_{in}$ [Sv]', fontsize=20)
ax1.set_xlim(0,0.4)

if plot:
  plt.savefig(param+'_2.png',format='png',dpi=300,bbox_inches='tight')
else:
  plt.savefig(param+'_1.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
