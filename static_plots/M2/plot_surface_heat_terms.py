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

path='ncfiles/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx = ['dx1','dx2','dx5','dx10']
dx2 = ['1 km','2 km','5 km','10 km']

dx = ['dx2','dx5','dx10']
dx2 = ['2 km','5 km','10 km']

param = 'surface_heat_terms_'
labels = ['-5.0','-2.5','0.0','2.5','5.0']
wind = [-5,-2.5,0,2.5,5]
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
colors = ['#6495ed','#ff6347','k','#3cb371','m']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

area = netCDF4.Dataset('/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx10km/M2_exp0/ocean_geometry.nc').variables['Ah'][:]
HI = netCDF4.Dataset('/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx10km/M2_exp0/ice_month.nc').variables['HI'][0,:]
area = np.ma.masked_where(HI.mask == True, area) # m^2
width = 0.35
variables = ['SW','LW','SH','LH']
variables = ['SW','LW']
for j in range(len(dx)):
  for v in range(len(variables)):
      # melt on
      txt = str('%s_1 = np.zeros(%d)' %(variables[v],len(wind))) # each variable with len = 5
      exec(txt)
      # melt off
      txt = str('%s_2 = np.zeros(%d)' %(variables[v],len(wind))) # each variable with len = 5
      exec(txt)
      
  for i in range(len(exps1)):
     file1 = netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_surface_heat_terms.nc') 
     file2 = netCDF4.Dataset(path+exps2[i]+'_'+dx[j]+'_surface_heat_terms.nc') 
     time = file2.variables['time'][:] # in yr
     tmp = np.nonzero((time>= 15.) & (time<20.))[0]
     # first melting on
     for v in range(len(variables)):
         if variables[v] == 'FRA': # negative
           txt1 = str("%s_1[%d] = -file1.variables['%s'][tmp].mean()/area.sum()" %(variables[v],i,variables[v])) # W/m^2
           txt2 = str("%s_2[%d] = -file2.variables['%s'][tmp].mean()/area.sum()" %(variables[v],i,variables[v])) # W/m^2
           exec(txt1); exec(txt2)
         else:
           txt1 = str("%s_1[%d] = file1.variables['%s'][tmp].mean()/area.sum()" %(variables[v],i,variables[v])) # W/m^2
           txt2 = str("%s_2[%d] = file2.variables['%s'][tmp].mean()/area.sum()" %(variables[v],i,variables[v])) # W/m^2
           exec(txt1); exec(txt2)
         
  # compute total and plot
  TOT_1 = np.zeros(len(wind)); TOT_2 = np.zeros(len(wind));
  f, ax = plt.subplots(len(variables)+1, sharex=True, figsize=(12,12))
  for v in range(len(variables)):
      if v == 0:
        ax[v].set_title(r'$\Delta x$ = %s' %(dx2[j]), fontsize=18)
        txt1 = str('ax[v].plot(wind, %s_1, "-o", color=colors[0], linewidth=2, label="Melting on")' %(variables[v]))
        txt2 = str('ax[v].plot(wind, %s_2, "-*", color=colors[1], linewidth=2, label="Melting off")' %(variables[v]))
        exec(txt1); exec(txt2)
        ax[v].legend(framealpha=0, fancybox=True)
      else:
        # plot
        txt = str('ax[v].plot(wind, %s_1, "-o", color=colors[0], linewidth=2)' %(variables[v]))
        exec(txt)
        txt = str('ax[v].plot(wind, %s_2, "-*", color=colors[1], linewidth=2)' %(variables[v]))
        exec(txt)

      # y label
      ax[v].set_ylabel('%s [W m$^{-2}$]' %(variables[v]))
      # total
      txt = str("TOT_1 = TOT_1[:] + %s_1[:]" %(variables[v]))
      exec(txt)
      txt = str("TOT_2 = TOT_2[:] + %s_2[:]" %(variables[v]))
      exec(txt)

  # plot total
  ax[v+1].plot(wind, TOT_1, '-o', color=colors[0], linewidth=3)
  ax[v+1].plot(wind, TOT_2, '-*', color=colors[1], linewidth=3)
  #ax[v+1].plot(wind, np.abs(TOT_1)-np.abs(TOT_2), '-*', color='k', linewidth=2)
  #ax[v+1].plot([-10,10], np.zeros(2), '-', color='k', linewidth=1)
  ax[v+1].set_xlabel(r'$U_{asf}$ [m s$^{-1}$]',fontsize=18)
  ax[v+1].set_ylabel(r'TOTAL [W m$^{-2}$]')
  ax[v+1].set_xlim(-5.5,5.5)
  plt.savefig(param+dx[j]+'.png',format='png',dpi=300,bbox_inches='tight')
  plt.show()

