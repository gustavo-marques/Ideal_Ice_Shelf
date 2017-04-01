from matplotlib import pyplot as plt
import netCDF4
import numpy as np
import sys
sys.path.append('/ncrc/home2/Gustavo.Marques/python/pyGVtools/')
import m6toolbox

# Define a function to plot a section
def plot_section(file_handle, record, xq, i=125, variable='temp',eta='e',yvar='yh',clim=(0,30), plot_grid=True, rep='pcm', xlim=(0,1000), ylim=(-4000,0), cmap=plt.cm.jet):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.

    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

    e = file_handle.variables[eta][record,:,:,i] # Vertical grid positions
    s = file_handle.variables[variable][record,:,:,i] # Scalar field to color
    y = file_handle.variables[yvar][:]
    x = file_handle.variables['xh'][:]
    #x,z,q = m6toolbox.section2quadmesh(xq, e, s, representation=rep) # This yields three areas at twice the model resolution
    ##cs = plt.pcolormesh(x, z, q, cmap=cmap);
    #plt.colorbar()
    ##cs.set_clim(clim)
    [Y,TMP] = np.meshgrid(y,e[:,0])
    print e.shape, Y.shape
    if plot_grid:
       plt.plot(Y.T, e.T, 'k', lw=0.5, hold=True)
       plt.plot(y,e[0,:],'k',lw=1); plt.plot(y,e[-1,:],'k',lw=1);
    if variable == 'v':
       plt.text(5,-500, r'max$(|v|)$:'+str(np.max(np.abs(s))), fontdict=font)
    plt.ylim(ylim)
    plt.xlim(xlim)
    return #cs

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx2km/Sigma_zstar/M2_exp2/'

lath = netCDF4.Dataset(path+'IDEAL_IS_IC.nc').variables['lath'][:] # km
yh = netCDF4.Dataset(path+'prog.nc').variables['yh'][:] # km
file   = netCDF4.Dataset(path+'prog.nc')
e   = netCDF4.Dataset(path+'prog.nc').variables['e'][-1,:,:,125]
h   = netCDF4.Dataset(path+'prog.nc').variables['h'][-1,:,:,125]
# This is the coordinate of the cell corners (u-points in 1D)
yq = netCDF4.Dataset(path+'IDEAL_IS_IC.nc').variables['latq'][:]
yq = np.concatenate(([0],yq)) # Inserts left most edge of domain in to coordinate

fig = plt.figure(figsize=(12,8))
ax1 = plt.subplot(111)
#cs1 = plot_section(file, 0, yq, variable='h', eta='e',yvar='yh')
plot_section(file, 0, yq, variable='h', eta='e',yvar='yh')
#plt.title(coords[n] + '(final)', fontsize=18)
#plt.tick_params(axis='x',labelbottom='off')
plt.xlabel('y [km]',fontsize=15); plt.ylabel('depth [m]',fontsize=15)
#cbar_ax = fig.add_axes([0.16, 0.025, 0.7, 0.05])
#cbar=fig.colorbar(cs1, orientation='horizontal', cax=cbar_ax,ticks=[-1.0,-0.5,0.0])
#cbar=fig.colorbar(cs1, orientation='horizontal')
#cbar.set_label(r'$\Delta z$ [m]', fontsize=16)
plt.savefig('sigma_shelf_zstar.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
