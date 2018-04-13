from matplotlib import pyplot as plt
import netCDF4
from scipy import stats
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

path='ncfiles/'
path1='TXT/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

#dx2 = ['1 km','2 km','5 km','10 km']
dx2 = ['2 km','5 km','10 km']
dx = ['dx2','dx5','dx10']
#dx = ['dx1']

rho_ice = 918.
rho0 = 1028.0
cp = 3974.0
cp = 9.0e-06
#cp = 3.2145778682161476e-08
#Lf = 334.0
Lf = 4.1e-05 #3.84e-05 
#Lf = 3.34e5

TF = np.array([0.06,0.09,0.07])
eta = [51.,85.0,96.0]

TR0 = np.array([0.02,0.01,0.01]) #Sv 

param = 'model_vs_theory'
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

f, ((ax1, ax2, ax3)) = plt.subplots(3, 1, figsize=(8,16))

for j in range(len(dx)):
  out1_mean = []; out1_std = []
  in1_mean = []; in1_std = []
  melt1_mean = []; 
  eff1_mean = []; 
  for i in range(len(exps1)):
     #print path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity.nc'
     #print path+exps2[i]+'_'+dx[j]+'_heat_budget_cavity.nc'
     # Hout 
     out1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hout'][:]))
     out1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hout'][:]))
     # Hin
     in1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hin'][:]))
     in1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hin'][:]))
     # melt
     melt1_mean.append(np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0])
     
     # (Hin + Hout)/Hin
     eff1_mean.append((in1_mean[i]+out1_mean[i])/in1_mean[i])

  melt1_mean = np.array(melt1_mean)*rho_ice # in kg/s
  # plot
  ax1.errorbar(melt1_mean, in1_mean, in1_std, linestyle='', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx2[j])
  ax2.errorbar(melt1_mean, out1_mean, out1_std, linestyle='', marker='o', color=colors[j], elinewidth=2)
  ax3.plot(melt1_mean, eff1_mean, 'o', color=colors[j], lw=2)
  # Fit
  M = np.linspace(melt1_mean.min(),melt1_mean.max(),10)
  print 'Hin'
  a, b, r_value, p_value, std_err = stats.linregress(melt1_mean, in1_mean)
  print 'a, b, r_value,  p_value, std_err', a, b, r_value,  p_value, std_err

  print 'a1, b1',(TF[j]*eta[j]), (TF[j])
  cccp = a/(TF[j]*eta[j]); print 'cccp',cccp
  #TR0[j] = b/(TF[j]); print 'TR0[j]',TR0[j]
  #cp = a/(TF[j]*eta[j]); print 'cp', cp

  print 'Hout'
  a, b, r_value, p_value, std_err = stats.linregress(melt1_mean, out1_mean)
  print 'a, b, r_value,  p_value, std_err', a, b, r_value,  p_value, std_err
  print 'a1, b1',-((TF[j]/rho0)*eta[j]-Lf), (TF[j]*TR0[j])
  print 'Lf',((cp*TF[j]*eta[j])-a)
  Hin = (TF[j] * TR0[j] + cp*TF[j]*eta[j]*M)
  #Hin = (TF[j]*eta[j]*M/rho0)
  Hout = (-TF[j]  * TR0[j]) + (cp*TF[j]*eta[j]-Lf)*M
  print Hout
  #Hout = -((TF[j]*eta[j]/rho0)-Lf)*M

  delta = (Hin-np.abs(Hout))/Hin #((TR_warm[j] * TR0[j]) + (cp*TR_warm[j]*eta_warm[j]*M_warm))

  ax1.plot(M,Hin, '-', color=colors[j], lw=2)
  ax2.plot(M,Hout, '-', color=colors[j], lw=2)
  ax3.plot(M,delta, '-', color=colors[j], lw=2)

ax1.legend(loc='best', fontsize=18, ncol=1)
ax1.set_ylabel(r'$H_{in}$ [TW]', fontsize=20)
ax2.set_ylabel(r'$H_{out}$ [TW]', fontsize=20)
ax3.set_ylabel(r'$\delta_{th}$', fontsize=20)
ax3.set_xlabel(r'Mass flux [kg s$^{-1}$]', fontsize=20)


ax1.set_xlim(0,70e3)
ax2.set_xlim(0,70e3)
ax3.set_xlim(0,70e3)

ax1.set_ylim(-0.1,2.6)
ax2.set_ylim(-2.6,0.1)
#ax3.set_ylim(-0.1,0.4)

plt.savefig(param+'.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
