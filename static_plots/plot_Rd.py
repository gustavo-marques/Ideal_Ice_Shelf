from matplotlib import pyplot as plt
import netCDF4
import numpy as np

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx2km/Sigma_zstar/M2_exp2/'
icfile = 'IDEAL_IS_IC.nc'

lath = netCDF4.Dataset(path+icfile).variables['lath'][:] # km
Rd   = netCDF4.Dataset(path+'prog.nc').variables['Rd_dx'][1000:,:].mean(axis=0).mean(axis=1) * 2.0

plt.figure(figsize=(10, 4))
plt.plot(lath,Rd,'k-')
plt.title('First baroclinic radius of deformation ($R_d$)', fontsize=18)
plt.xlabel(r' y [km]', fontsize=16)
plt.ylabel(r'$R_d$ [km]', fontsize=16)
plt.savefig('Rd.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
