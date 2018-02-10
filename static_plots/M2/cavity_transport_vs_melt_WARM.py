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
exps1 = ['M2_exp4','M2_exp3','M2_exp2'] # melt on
exps2 = ['M2_exp14','M2_exp17','M2_exp16'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx1 = ['1km','2km','5km','10km']

param = 'cavity_transport_vs_melt_WARM'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
tx=[0.25,0.2,0.15,0.1]
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1)) = plt.subplots(1, 1, figsize=(9,7))

for j in range(len(dx)):
  data1_mean = []
  data2_mean = []
  for i in range(len(exps1)):
     data1 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_cavity_transports.txt')[2]
     data1_mean.append(np.abs(data1)/1.0e6)
     data2 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0]
     data2_mean.append(data2) # Gt
  # plot
  xx=np.array(data2_mean)*1.0e3/(365*3600*24.)# in Sv
  yy=np.array(data1_mean)
  m,b = np.polyfit(xx,yy, 1) # sv and sv to print
  m = int(m)
  m1,b1 = np.polyfit(np.array(data2_mean),yy, 1) # gt/yr and sv to plot
  melt = np.linspace(np.min(data2_mean),np.max(data2_mean),10)
  y=m1*melt+b1
  ax1.plot(data2_mean, data1_mean, linestyle='', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j])
  # FIT
  ax1.plot(melt,y,linestyle='-', color=colors[j], lw=1)
  # TEXT
  ax1.text(2.5,tx[j],r'$\eta_{dyn}$ = '+str(m), color=colors[j], fontsize=20)

#ax1.legend(loc='upper center', fontsize=20, ncol=2)
ax1.set_ylabel(r'Cavity inflow [Sv]', fontsize=20)
ax1.set_xlabel(r'Melt rate [Gt yr$^{-1}$]', fontsize=20)

ax1.set_xlim(0,45)
ax1.set_ylim(0.0,0.35)

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
