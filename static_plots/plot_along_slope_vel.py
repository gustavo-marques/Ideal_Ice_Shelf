import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

path = '/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M2_exp2/out1/'
u = Dataset(path+'ocean_month_z.nc').variables['u'][-1]
y = Dataset(path+'ocean_month_z.nc').variables['yh'][:]
x = Dataset(path+'ocean_month_z.nc').variables['xh'][:]
z = Dataset(path+'ocean_month_z.nc').variables['zt'][:]
depth = Dataset(path+'ocean_geometry.nc').variables['D'][:]
ssh = Dataset(path+'ocean_sfc.nc').variables['SSH'][-1,:]
[Y,Z] = np.meshgrid(y,z)
[XX,YY] = np.meshgrid(x,y)

fig, ax = plt.subplots()
ax.contour(XX,YY,depth,[0,500,600,1000,2000,3000],colors = 'black')
ax.plot(XX[:,200]*np.ones(len(y)),YY[:,100],'g',lw=1.5)
ax.set_xlabel('x [km]',fontsize=15)
ax.set_ylabel('y [km]',fontsize=15)
plt.savefig('vel_sec_x100.png')

fig, ax = plt.subplots()
ax.contour(XX,YY,depth,[0,500,600,1000,2000,3000],colors = 'black')
ax.plot(XX[:,300]*np.ones(len(y)),YY[:,150],'g',lw=1.5)
ax.set_xlabel('x [km]',fontsize=15)
ax.set_ylabel('y [km]',fontsize=15)
plt.savefig('vel_sec_x150.png')

vel = np.linspace(-0.5,0.5,100)
fig, ax = plt.subplots()
ct = ax.contourf(Y,-Z,u[:,:,200],vel, cmap=plt.cm.bwr,extend="both")
ct.cmap.set_under('b')
ct.cmap.set_over('r')
cbar = fig.colorbar(ct, ticks=[-0.4, -0.2, 0, 0.2, 0.4], orientation = 'horizontal')
cbar.set_label('m/s')
ax.contour(Y,-Z,u[:,:,200],[0],linestyles='dashed',colors = 'gray')
ax.set_xlim(0,800)
ax.set_ylim(-2000,0)
ax.set_title('Along-slope vel. [m/s] at x = 200 km')
ax.set_xlabel('y [km]',fontsize=15)
ax.set_ylabel('depth [m]',fontsize=15)
ax.plot(y,-depth[:,200],'k',lw=1)
ax.plot(y,ssh[:,200],'k',lw=1)
plt.savefig('M2_exp2_1km_x100.png')
plt.show()

fig, ax = plt.subplots()
ct = ax.contourf(Y,-Z,u[:,:,300],vel, cmap=plt.cm.bwr,extend="both")
ct.cmap.set_under('b')
ct.cmap.set_over('r')
cbar = fig.colorbar(ct, ticks=[-0.4, -0.2, 0, 0.2, 0.4], orientation = 'horizontal')
cbar.set_label('m/s')
ax.set_xlim(0,800)
ax.contour(Y,-Z,u[:,:,300],[0],linestyles='dashed',colors = 'gray')
ax.set_ylim(-2000,0)
ax.set_title('Along-slope vel. [m/s] at x = 300 km')
ax.set_xlabel('y [km]',fontsize=15)
ax.set_ylabel('depth [m]',fontsize=15)
ax.plot(y,-depth[:,300],'k',lw=1)
ax.plot(y,ssh[:,300],'k',lw=1)
plt.savefig('M2_exp2_1km_x150.png')
plt.show()
