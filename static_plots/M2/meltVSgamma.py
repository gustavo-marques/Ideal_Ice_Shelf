import numpy as np
import matplotlib.pyplot as plt

#gamma_h10 = [0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01, 0.015]
gamma_h10 = [0.001, 0.002, 0.003, 0.004, 0.005, 0.0075]
#melt_h10 = np.array([842081.0, 1.07734e+06, 1.23636e+06, 1.30643e+06, 1.36153e+06, 1.46125e+06, 1.4963e+06, 1.54164e+06]) * 1.0e-12 * 3600. * 24 * 365 # in Gt/yr
melt_h10 = np.array([1.35282e+06, 1.66366e+06, 1.79328e+06, 1.88344e+06, 1.92375e+06, 2.01395e+06]) * 1.0e-12 * 3600. * 24 * 365 # in Gt/yr

print melt_h10
#coeffs = np.polyfit(gamma_h10,melt_h10,3) # third-order poly
coeffs = np.polyfit(gamma_h10,melt_h10,4) # fourth-order poly
#coeffs = np.polyfit(gamma_h10,melt_h10,2) # second-order poly
new_gamma = np.arange(0.001,0.0075,0.0001)
new_melt = coeffs[0]*new_gamma**4 + coeffs[1]*new_gamma**3 + coeffs[2]*new_gamma**2 + coeffs[3]*new_gamma + coeffs[4]
#new_melt = coeffs[0]*new_gamma**3 + coeffs[1]*new_gamma**2  + coeffs[2]*new_gamma + coeffs[3]
#new_melt = coeffs[0]*new_gamma**2 + coeffs[1]*new_gamma  + coeffs[2]
# find values that yields melt~31
#tmp = np.nonzero(new_melt<=34.7)[-1][-1]
tmp = np.nonzero(new_melt<=61.3)[-1][-1]
print('Gamma is: '+str(new_gamma[tmp]))

plt.figure()
#plt.plot(gamma_com,melt_com,'ko')
plt.plot(gamma_h10,melt_h10,'ko')
plt.plot(new_gamma,new_melt,'k-',lw=1.0)
plt.plot(new_gamma[tmp],new_melt[tmp],'b*',markersize=12)
plt.xlabel(r'$\Gamma_T$',fontsize=18)
plt.ylabel('Mean melt rate [GT yr$^{-1}$]',fontsize=18)
plt.xlim(0.0,0.008)
plt.savefig('melt_vs_gamma.png',format='png',dpi=300,bbox_inches='tight')
plt.show()

#plt.figure()
#plt.plot(gamma_com,melt_com,'k-o',label='COM',lw=1)
#plt.plot(gamma_tpy,melt_tpy,'r-o',label='TPY',lw=1)
#plt.legend(loc='upper left', shadow=True)
#plt.xlabel(r'$\Gamma_T$')
#plt.ylabel('mean melt rate (m a$^{-1}$)')
#plt.xlim(0,3.01)
#plt.savefig('melt_gamma_COMvsTPY.png')

