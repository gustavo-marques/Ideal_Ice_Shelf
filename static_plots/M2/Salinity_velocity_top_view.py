from matplotlib import pyplot as plt
import netCDF4
import numpy as np
import sys,os
sys.path.append('/work/gmm/pyGVtools/')
sys.path.append('colormap/')
import colormaps as cmaps
import m6toolbox
#plt.style.use('ggplot')
import matplotlib
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams.update({'font.size': 18})

class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""

    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels-self.levels.min()
        self._x/= self._x.max()
        self._y = np.linspace(0, 1, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi, alpha)

# Define a function to plot a xy map
def plot_xy(x, y, ax, data1, data2, u, v, depth, j=0, eta='e',yvar='yh', xvar='xh', xlim=(0,500), ylim=(0,500), cmap=cmaps.viridis, hidex=False, hidey=False, m1=-2, m2=2):
    """Plots a xy map reading scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.

    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 18,
        }

    #levels = np.linspace(m1,m2,30)
    levels = [33.78, 33.785, 33.79, 33.795, 33.80, 33.805, 33.81, 33.82, 33.83, 33.84, 33.85, 33.86, 33.87, 33.88, 
              33.89, 33.9, 33.91, 33.92, 33.93, 33.94, 33.95, 33.96, 33.97,
              33.98, 33.99, 34.0, 34.1, 34.15, 34.2, 34.25, 34.3, 34.4, 34.5, 
              34.6, 34.7, 34.8, 35]
    cmap_nonlin = nlcmap(cmap, levels)
    
    cs = ax.contourf(x,y,data1,levels, cmap=cmap_nonlin)
    ax.contour(x,y,depth,[550, 730], colors='white',lw=5)
    n = 10
    q = ax.quiver(x[::n,::n], y[::n,::n], u[::n,::n], v[::n,::n], scale=0.45, scale_units='inches', color='black')
    if j == 0: # 
      ax.quiverkey(q, X=410, Y=50, U=0.1, label=r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='data')
      #qk = ax.quiverkey(q, 400, 50, 1, r'$1 \frac{m}{s}$', labelpos='E')

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    if not hidex:
      #ax.axes.get_xaxis().set_ticks([])
      ax.set_xlabel('y [km]')

    if not hidey:
      #ax.axes.get_yaxis().set_ticks([])
      ax.set_ylabel('Depth [m]')

    return cs

path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx1km/'
dx = 1.0
icfile = '/IDEAL_IS_IC.nc'
exps = ['M2_exp4','M2_exp14','M2_exp0','M2_exp13'] 
labels = ['on','off','on','off']
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
lonh = netCDF4.Dataset(path+'M2_exp4/ocean_geometry.nc').variables['geolon'][:] # km
lath = netCDF4.Dataset(path+'M2_exp4/ocean_geometry.nc').variables['geolat'][:] # km
depth = netCDF4.Dataset(path+'M2_exp4/ocean_geometry.nc').variables['D'][:] # km
# base of STATIC ice shelf, which is ssh(t=0); make it positive
ice_base = -netCDF4.Dataset(path+exps[0]+icfile).variables['ave_ssh'][0,:,:]
ice_base[ice_base<1e-5] = 0.0
u = np.zeros(depth.shape); v = np.zeros(depth.shape)
f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10,12))

for j in range(len(exps)):
  if dx == 1.0:  
    path_to_file = path+exps[j]+'/out/prog*.nc'
  else:
    path_to_file = path+exps[j]+'/prog.nc'
  print 'path:',path_to_file
  print '--------------------------------- \n'
  # plot
  if j == 0:
    ax = ax1; var1='salt'; var2='Rd_dx'; axm=1; hidex=True; hidey=False; time=3650
  elif j == 1:
    ax = ax2; var1='salt'; var2='Rd_dx'; axm=2; hidex=True; hidey=True; time=3650
  elif j == 2:
    ax = ax3; var1='salt'; var2='Rd_dx'; axm=1; hidex=False; hidey=False; time=3650
  else:
    ax = ax4; var1='salt'; var2='Rd_dx'; axm=2; hidex=False; hidey=True; time=3650

  file_handle = netCDF4.MFDataset(path_to_file)

  t_index = np.nonzero(file_handle.variables['time'][:] <= time)[0][-1]
  data1 = file_handle.variables[var1][t_index,1,:] # surface
  u1 = file_handle.variables['u'][t_index,1,:] # surface
  v1 = file_handle.variables['v'][t_index,1,:] # surface
  # mask grounded ice
  data1 = np.ma.masked_where(ice_base+1.0>=depth, data1) # need +1 here
  # u and v at h points
  utmp = 0.5 * (u1[:,0:-1] + u1[:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
  vtmp = 0.5 * (v1[0:-1,:] + v1[1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
  u[:,1::] = utmp; u[:,0] = 0.5*u1[:,0] #u_i=1 = 0.5*u_(i=3/2)
  v[1::,:] = vtmp; v[0,:] = 0.5*v1[0,:] #v_j=1 = 0.5*v_(j=3/2)
  # mask
  u = np.ma.masked_where(data1.mask==True, u)
  v = np.ma.masked_where(data1.mask==True, v)

  data2 = file_handle.variables[var2][t_index,:] * dx # Rd in km
  if j == 0 or j == 2 :
    cs1 = plot_xy(lonh, lath, ax, data1, data2, u, v, depth, j, hidex=hidex, hidey=hidey, m1=33.6, m2=34.7)
  else:
    cs2 = plot_xy(lonh, lath, ax, data1, data2, u, v, depth, j, hidex=hidex, hidey=hidey, m1=33.6, m2=34.7)

  if j == 3: plt.tight_layout(w_pad=0.5, h_pad=0.5) ;f.subplots_adjust(bottom=0.11) 
  
  # aspect ratio
  ratio = 1.0
  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()
  ax.set_aspect(np.abs((xmax-xmin)/(ymax-ymin))*ratio, adjustable='box-forced')

  # text
  ax.text(10,190,abcd[j], color='b', fontsize=18)
  ax.text(8,50,'Melting', color='b', fontsize=18)
  ax.text(40,10,labels[j], color='b', fontsize=18)
  # fill IS and plot bathymetry
  #ax.fill_between(lonh[0,:], 0.0, 50.0, facecolor='white', interpolate=True)
  #ax.fill_between(yh, e[-1,:], -4000, facecolor='lightgray', interpolate=True)
  #ax.plot(yh, ssh,'k-',lw=1.5)

f.text(0.3,0.94,r'Upwelling (U$_{asf}$ = 5.0 m s$^{-1}$)', color='k', fontsize=24)
f.text(0.3,0.53,r'Downwelling (U$_{asf}$ = -5.0 m s$^{-1}$)', color='k', fontsize=24)

cbar_ax1 = f.add_axes([0.128, 0.07, 0.8, 0.03])
#cbar1=f.colorbar(cs1, orientation='horizontal', cax=cbar_ax1, ticks=[33.8,34.0,34.2,34.4,34.6])
cbar1=f.colorbar(cs1, orientation='horizontal', cax=cbar_ax1)
cbar1.set_label(r'Salinity', fontsize=18)
plt.savefig('salinity_velocity_top_view.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
