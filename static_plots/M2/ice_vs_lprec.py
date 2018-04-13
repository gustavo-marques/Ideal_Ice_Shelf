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
path3='ncfiles/'
icfile = 'IDEAL_IS_IC.nc'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx1 = ['1km','2km','5km','10km']
dx = ['dx2','dx5','dx10']
dx1 = ['2km','5km','10km']

param = 'ice_vs_lprec'
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1)) = plt.subplots(1, 1, figsize=(9,7))

for j in range(len(dx)):
  lprec_ratio = []
  ice_ratio = []
  for i in range(len(exps1)):
     path_to_file1 = path2+dx[j]+'km/'+exps1[i]+'/'+exps1[i]+'_'+dx1[j]+'_buoyancy.nc'
     path_to_file2 = path2+dx[j]+'km/'+exps2[i]+'/'+exps2[i]+'_'+dx1[j]+'_buoyancy.nc'
     print 'Path1, path2',path_to_file1, path_to_file2
     print '--------------------------------- \n'
     data1 = netCDF4.Dataset(path_to_file1).variables['seaiceVolume'][:]/1.0e12 # 1000 km^3
     data2 = netCDF4.Dataset(path_to_file2).variables['seaiceVolume'][:]/1.0e12
     ice_ratio.append(data1.mean()/data2.mean())
 
     data1 = netCDF4.Dataset(path3+exps1[i]+'_'+dx[j]+'_net_surface_water_flux.nc').variables['NET'][:]
     data2 = np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0] # Gt
     #lprec_ratio.append(data2/np.abs(data1.mean())) # Gt
     lprec_ratio.append(63./data2) # Gt/Gt

  # plot
  ax1.plot(lprec_ratio, ice_ratio, linestyle='', marker='o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx1[j])


ax1.legend(loc='upper left', fontsize=14, ncol=2)
ax1.set_ylabel(r'H* ', fontsize=20)
ax1.set_xlabel(r'LPREC*', fontsize=20)

#ax1.set_xlim(4,10)
#ax1.set_ylim(0.06,0.11)

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
