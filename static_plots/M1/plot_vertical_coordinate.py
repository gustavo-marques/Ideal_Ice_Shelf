from matplotlib import pyplot as plt
import netCDF4
import numpy as np
import sys
sys.path.append('/work/gmm/pyGVtools/')
import m6toolbox
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

# Define a function to plot a section
def plot_section(file_handle, record, xq, ax, i=125, variable='temp',eta='e',yvar='yh',clim=(0,30), plot_grid=True, rep='pcm', xlim=(0,1000), ylim=(-4000,0), cmap=plt.cm.jet, hide_ticks=False):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.

    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

    import matplotlib
    matplotlib.rcParams.update({'font.size': 16})
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
       #ax.plot(Y.T, e.T, 'k', lw=0.5, hold=True)
       ax.plot(Y.T, e.T, 'k', lw=0.5)
       ax.plot(y,e[0,:],'k',lw=1); plt.plot(y,e[-1,:],'k',lw=1);
    if variable == 'v':
       ax.text(5,-500, r'max$(|v|)$:'+str(np.max(np.abs(s))), fontdict=font)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    if hide_ticks:
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])

    return #cs

path='/archive/gmm/Ideal_ice_shelf/Mode1/dx2km/Sigma_zstar/M1_exp5/out1/'

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
plot_section(file, 0, yq, ax1, variable='h', eta='e',yvar='yh')
# box for zoom
ax1.plot([220,5],[-920,-920],'b',linewidth=2)
ax1.plot([220,5],[-5,-5],'b',linewidth=2)
ax1.plot([220,220],[-920,-5],'b',linewidth=2)
ax1.plot([5,5],[-920,-5],'b',linewidth=2)
# zoom
ax2 = add_subplot_axes(ax1,[.05, .15, .45, .5])
plot_section(file, 0, yq, ax2, variable='h', eta='e',yvar='yh',xlim=(0,220), ylim=(-920,0),hide_ticks=True)
# lines
ax1.plot([5,50],[-920,-3412],'b',linewidth=2)
ax1.plot([220,500],[-5,-1380],'b',linewidth=2)
ax1.set_xlabel('y [km]',fontsize=18); ax1.set_ylabel('depth [m]',fontsize=18)
#fig.subplots_adjust(left=0.05)
# resolution
ax1.annotate(r'$\Delta z$ = 10 m', xy=(0.8, 1.05), xytext=(0.8, 1.05), xycoords='axes fraction', 
            fontsize=15, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'))

ax1.annotate(r'$\Delta z$ = 25 m', xy=(0.8, 0.85), xytext=(0.8, 0.85), xycoords='axes fraction', 
            fontsize=15, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'))

ax1.annotate(r'$\Delta z$ = 110 m', xy=(0.8, 0.4), xytext=(0.8, 0.4), xycoords='axes fraction', 
            fontsize=15, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'))

plt.savefig('sigma_shelf_zstar.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
