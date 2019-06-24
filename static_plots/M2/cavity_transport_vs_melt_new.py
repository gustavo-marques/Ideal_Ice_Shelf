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
plot = False
plot = True
path1='TXT/'
path2='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx1 = ['1km','2km','5km','10km']

param = 'cavity_transport_vs_melt'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
symbs = ['s','s','*','o','o']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
tx=[0.10,0.12,0.14,0.16]
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1)) = plt.subplots(1, 1, figsize=(10,8))

for j in range(len(dx)):
  data1_mean = []
  data2_mean = []
  print 'dx = ',dx[j]
  for i in range(len(exps1)):
     data1 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_cavity_transports.txt')[2]
     data1_mean.append(np.abs(data1)/1.0e6)
     data2 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0]
     data2_mean.append(data2) # Gt

  # plot
  xx=np.array(data2_mean)*1.0e3/(365*3600*24.)# in Sv
  yy=np.array(data1_mean)
  #m,b = np.polyfit(xx,yy, 1) # sv and sv to print
  m, b, r_value, p_value, std_err = stats.linregress(xx, yy)
  print m, b, r_value,  p_value, std_err
  m = int(m)
  r_value = format(r_value, '.2f')
  melt = np.linspace(np.min(data2_mean),np.max(data2_mean),10)
  m1, b1, r_value1, p_value, std_err = stats.linregress(np.array(data2_mean), yy) # for plotting
  print 'background transport: ',b1
  melt = np.linspace(np.min(data2_mean),np.max(data2_mean),10)
  y=m1*melt+b1

  for i in range(len(exps1)):
      ax1.plot(data2_mean[i], data1_mean[i], linestyle='', marker=symbs[i], color=colors[j], markersize=10)

  ax1.plot(data2_mean, data1_mean, linestyle='-', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j])
 
  #ax1.plot(data2_mean, data1_mean, linestyle='', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j], markersize=10)
  if plot: 
    ax1.plot(melt,y,linestyle='--', color=colors[j], lw=2)
    # TEXT
    ax1.text(41,tx[j],r'$\eta_{dyn}$ = '+str(m)+' (r = '+str(r_value)+')', color=colors[j], fontsize=20)

# wind legened
w1=[2,4,6,8]
w2=[0.30,0.32,0.34]
symbs = ['s','*','o']
ms=[10,15,10]
for i in range(len(w1)):
  for j in range(len(w2)):
    ax1.plot(w1[i], w2[j], linestyle='', marker=symbs[j], color=colors[i], markersize=ms[j])


ax1.legend(loc='upper left', fontsize=18, ncol=2)
ax1.set_ylabel(r'$\Phi_{in}$ [Sv]', fontsize=20)
ax1.set_xlabel(r'Total melt flux [Gt/yr]', fontsize=20)

ax1.set_xlim(0,68)
ax1.set_ylim(0.075,0.4)
if plot:
  plt.savefig(param+'_2.png',format='png',dpi=300,bbox_inches='tight')
else:
  plt.savefig(param+'_1.png',format='png',dpi=300,bbox_inches='tight')

plt.show()
