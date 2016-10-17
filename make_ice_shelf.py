# setup initial T/S

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# setup params
#Gustavo,
#
#Here is a slightly different version but with constants
#
#H(x) = H0 *(Q0)^(1/4) / [Q0+4*H0^4*C^3*x]^(1/4)
#
#where H0 is ice thickness at the grounding line Q0 = U0*H0 is ice flux at the grounding line, C = rho*g*(1-rho/rho_w)/4/Bbar, rho is the ice density, rho_w is the sea-water density, Bbar is ice stiffness parameter.
#
#for the following parameters 
#
#rho = 917 kg/m^3
#rho_w = 1028 kg/m^3
#Bbar = 1.68e8
#C = 1.4440e-06
#
#U0= 700 m/yr = 2.2182e-05 m/s
#H0 = 1500 m
#Q0 = 0.0333 m^2/s
#
#rho = 917. #kg/m^3
#rho_w = 1028. #kg/m^3
#Bbar = 1.68e8

C = 1.4440e-06
H0 = 1500.0 #m
Q0 = 0.03327 #m^2/s
ISL = 250e3 # ice shelf lenght
gp = 20.e3 # grouding line position
# IS file
file = Dataset('IC_IS.nc', 'r+')
y = Dataset('ic_ts.nc').variables['LAT'][:] * 1.0E3 # in m
dy = y[1]-y[0]
x = Dataset('ic_ts.nc').variables['LON'][:] * 1.0E3 # in m
dx = x[1]-x[0]
thick = file.variables['thick'][:]

h =  H0 *(Q0)**(1./4.) / (Q0+100*H0**4*C**3*(y-gp))**(1./4.)
h[y>=ISL] = 0.0
h[y<gp] = H0

# smooth
h_smooth = gaussian_filter(h,1)
h_smooth[y>=ISL] = 0.0
h_smooth[y<gp] = H0

plt.plot(y,h,'k',y,h_smooth,'r')
plt.show()

NY, NX = thick.shape
area = np.ones(NY) * dx * dy
area[h_smooth == 0.0] = 0.0
# write into nc file
for i in range(NX):
        file.variables['thick'][:,i] = h_smooth[:]
        file.variables['area'][:,i] = area[:]
# close
file.sync()
file.close()
