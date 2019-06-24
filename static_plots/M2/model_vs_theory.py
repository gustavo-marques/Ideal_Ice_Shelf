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
ty1=[1.5,1.8,2.1,2.4]
ty2=[-2.4,-2.1,-1.8,-1.5]
# plot some metrics for runs with varing wind forcing
plot = False
plot = True
path='ncfiles/'
path1='TXT/'
exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off

dx2 = ['1 km','2 km','5 km','10 km']
dx = ['dx1','dx2','dx5','dx10']

rho0 = 1028.0
Lf = 3.34e5

TF = np.array([9.30, 12.06,11.11,11.79]) # kg/(m s^2) ?
eta = np.array([60.99,51.68,75.75,106.03]) # nondim

TR0 = np.array([0.142e6,0.135e6,0.121e6,0.167e6])/1.0e12 # e6 to convert to m^3/2; 1.0e12 to TW 

param = 'model_vs_theory'
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

f, ((ax1, ax2, ax3)) = plt.subplots(3, 1, figsize=(8,16))


for j in range(len(dx)):
  out1_mean = []; out1_std = []
  in1_mean = []; in1_std = []
  melt1_mean = []; 
  eff1_mean = []; 
  for i in range(len(exps1)):
     # Hout 
     out1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hout'][:]))
     out1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hout'][:]))
     # Hin
     in1_mean.append(-np.mean(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hin'][:]))
     in1_std.append(np.std(netCDF4.Dataset(path+exps1[i]+'_'+dx[j]+'_heat_budget_cavity_monthly.nc').variables['Hin'][:]))
     # melt
     melt1_mean.append(np.loadtxt(path1+exps1[i]+'_'+dx[j]+'_melt_rate.txt')[0]) # Gt/yr
     
     # (Hin + Hout)/Hin
     eff1_mean.append((in1_mean[i]+out1_mean[i])/in1_mean[i])

  melt1_mean = np.array(melt1_mean)
  melt1 = (melt1_mean)/(1.0e-12 * 3600. * 24 * 365) # in kg/s

  # plot
  ax1.errorbar(melt1_mean, in1_mean, in1_std, linestyle='', marker='o', color=colors[j], elinewidth=2, label=r'$\Delta x$ = '+dx2[j])
  ax2.errorbar(melt1_mean, out1_mean, out1_std, linestyle='', marker='o', color=colors[j], elinewidth=2)
  ax3.plot(melt1_mean, eff1_mean, 'o', color=colors[j], lw=2, label=r'$\Delta x$ = '+dx2[j])
  # Fit
  M = np.linspace(melt1.min(),melt1.max(),10)
  print 'Hin #######'
  a, b, r_value, p_value, std_err = stats.linregress(melt1, in1_mean)
  print 'a, b, r_value,  p_value, std_err', a, b, r_value,  p_value, std_err
  print 'a1, b1',(TF[j]*eta[j]/(rho0*1.0e6)), (TF[j]*TR0[j])

  print 'Hout #######'
  a, b, r_value, p_value, std_err = stats.linregress(melt1, out1_mean)
  print 'a, b, r_value,  p_value, std_err', a, b, r_value,  p_value, std_err
  print 'a1, b1',-((TF[j])*eta[j]/(rho0*1.0e6)-Lf), (TF[j]*TR0[j])

  Hin = ((TF[j] * TR0[j]) + (TF[j]*eta[j]*M/(rho0*1.0e6)))
  print 'Hin', Hin
  print 'in1_mean', in1_mean
  Hout = (-TF[j]  * TR0[j]) - ((TF[j]*eta[j]/(rho0*1.0e6)-Lf/1.0e12)*M)
  print 'Hout',Hout
  print 'out1_mean', out1_mean
  # to be used in the RMSE
  Hin1 = ((TF[j] * TR0[j]) + (TF[j]*eta[j]*melt1/(rho0*1.0e6)))
  Hout1 = (-TF[j]  * TR0[j]) - ((TF[j]*eta[j]/(rho0*1.0e6)-Lf/1.0e12)*melt1)

  delta = (Hin-np.abs(Hout))/Hin #((TR_warm[j] * TR0[j]) + (cp*TR_warm[j]*eta_warm[j]*M_warm))
  M = M * 1.0e-12 * 3600. * 24 * 365 # Gt/yr

  #print 'shapes', np.array(in1_mean).shape(), np.array(Hin).shape()
  # RMSE
  rms_in = format(rmse(in1_mean, Hin1), '.2f')
  rms_out = format(rmse(out1_mean, Hout1), '.2f')
  print 'RMSE in/out',rms_in, rms_out

  if plot:
    # plot
    ax1.plot(M,Hin, '-', color=colors[j], lw=2)
    ax2.plot(M,Hout, '-', color=colors[j], lw=2)
    ax3.plot(M,delta, '-', color=colors[j], lw=2)
    # text
    ax1.text(2.5,ty1[j],r'RMSE = '+str(rms_in), color=colors[j], fontsize=20 )
    ax2.text(2.5,ty2[j],r'RMSE = '+str(rms_out), color=colors[j], fontsize=20 )

ax3.legend(loc='best', fontsize=16, ncol=1)
ax1.set_ylabel(r'$H_{in}$ [TW]', fontsize=20)
ax2.set_ylabel(r'$H_{out}$ [TW]', fontsize=20)
ax3.set_ylabel(r'$\delta_{th}$', fontsize=20)
ax3.set_xlabel(r'Mass flux [Gt yr$^{-1}$]', fontsize=20)


ax1.set_xlim(0,70)
ax2.set_xlim(0,70)
ax3.set_xlim(0,70)
ax1.set_ylim(-0.1,2.7)
ax2.set_ylim(-2.7,0.1)

if plot:
  plt.savefig(param+'_2.png',format='png',dpi=300,bbox_inches='tight')
else:
  plt.savefig(param+'_1.png',format='png',dpi=300,bbox_inches='tight')

plt.show()
