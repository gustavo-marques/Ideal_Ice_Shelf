from matplotlib import pyplot as plt
import netCDF4
import numpy as np

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp5/out3/'
icfile = 'IDEAL_IS_IC.nc'

lath = netCDF4.Dataset(path+icfile).variables['lath'][:] # km
lonh = netCDF4.Dataset(path+icfile).variables['lonh'][:] # km
Rd   = netCDF4.Dataset(path+'prog.nc').variables['Rd_dx'][700:,:].mean(axis=0)
Rd_m   = netCDF4.Dataset(path+'prog.nc').variables['Rd_dx'][700:,:].mean(axis=0).mean(axis=1)

plt.figure(figsize=(10, 10))
plt.subplot(2, 1, 1)
plt.contourf(lonh,lath,Rd,50)
plt.colorbar(orientation = 'horizontal')
plt.title('First baroclinic radius of deformation ($R_d$)', fontsize=18)
plt.xlabel(r' x [km]', fontsize=16)
plt.ylabel(r' y [km]', fontsize=16)
plt.subplot(2, 1, 2)
plt.plot(lath,Rd_m,'k-')
plt.xlabel(r' y [km]', fontsize=16)
plt.ylabel(r'$R_d$ [km]', fontsize=16)
plt.savefig('Rd.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
