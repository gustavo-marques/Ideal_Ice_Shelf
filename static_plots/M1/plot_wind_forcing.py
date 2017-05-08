from matplotlib import pyplot as plt
import netCDF4
import numpy as np
import sys
import matplotlib
sys.path.append('../../')
matplotlib.rcParams.update({'font.size': 18})

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp5/INPUT/forcing_10.nc'

lat = netCDF4.Dataset(path).variables['LAT'][:] # km
lon = netCDF4.Dataset(path).variables['LON'][:] # km
v_10   = netCDF4.Dataset(path).variables['V_10'][:]
u_10   = netCDF4.Dataset(path).variables['U_10'][0,:]
u2   = netCDF4.Dataset('/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp9/INPUT/forcing_10.nc').variables['U_10'][0,:]

wet = netCDF4.Dataset('/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp5/out1/ocean_geometry.nc').variables['wet'][:]

ice_shelf = netCDF4.Dataset('/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp5/out1/MOM_Shelf_IC.nc').variables['shelf_area'][:]

wind = np.ma.masked_where(wet == 0, v_10[0,:])
wind = np.ma.masked_where(ice_shelf[0,:] > 0, wind)

colors = [(255,255,255), (0,0,255), (51,255,51), (255,255,51), (255,0,0), (255,0,255)]
my_cmap = make_cmap(colors, bit=True)

lim = np.linspace(0,wind.max(),100)
#fig = plt.figure(facecolor='white',figsize=(8, 10))
f, (ax2, ax) = plt.subplots(1, 2, sharey=True, figsize=(10, 8))
#ax = plt.subplot(211,axisbg='gray')
#ax = fig.add_subplot(111,axisbg='gray')
cs = ax.contourf(lon,lat,wind,lim,cmap = my_cmap)
ax.set_aspect('equal')
ax.set_axis_bgcolor('gray')
cbar = plt.colorbar(cs,orientation = 'vertical',ticks=[0.0,5.0,10.0,15.0])
ax.set_title(r'max($\tau_y$) [m s$^{-1}$]',fontsize=18)
ax.fill([0,0,500,500],[900,1000,1000,900], fill=False, hatch='\\')
ax.plot(lon,np.ones(len(lon))*460,'k--')
ax.plot(lon,np.ones(len(lon))*720,'k--')
#plt.title('First baroclinic radius of deformation ($R_d$)', fontsize=18)
ax.set_xlabel(r' x [km]', fontsize=18)
#ax.set_ylabel(r' y [km]', fontsize=18)

#ax2 = plt.subplot(212,axisbg='gray',sharey=ax)
ax2.plot(u_10[:,0],lat,'k')
ax2.plot(u2[:,0],lat,'k--')
#ax2.set_xlabel(r' x [km]', fontsize=18)
ax2.set_ylabel(r' y [km]', fontsize=18)
labels = [item.get_text() for item in ax2.get_xticklabels()]
labels[0] = r'$\tau_x$'
labels[-1] = r'-$\tau_x$'
ax2.set_xticklabels(labels)
plt.savefig('wind_forcing.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
