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

matplotlib.rcParams.update({'font.size': 18})

# Define a function to plot a section
def plot_section(file_handle, xq, ax, data1, data2,im=250, eta='e',yvar='yh', rep='pcm', xlim=(0,650), ylim=(-2000,0), cmap=plt.cm.bwr, hidex=False, hidey=False, m1=-2, m2=2):
    """Plots a section of by reading vertical grid and scalar variable and super-sampling
    both in order to plot vertical and horizontal reconstructions.

    Optional arguments have defaults for plotting salinity and overlaying the grid.
    """
    font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 18,
        }

    e = file_handle.variables[eta][-24::,:,:,im].mean(axis=0) # Vertical grid positions
    y = file_handle.variables[yvar][:]
    x,z,q = m6toolbox.section2quadmesh(xq, e, data1, representation=rep) # This yields three areas at twice the model resolution
    #cs = ax.pcolormesh(x, z, q,norm=LogNorm(vmin=1, vmax=110), cmap=cmap)
    cs = ax.pcolormesh(x, z, q, vmin=m1, vmax=m2, cmap=cmap)
    if (len(data2.shape)>1):
      z = 0.5*(e[0:-1,:]+e[1:,:])
      [Y,TMP] = np.meshgrid(y,z[:,0])
      s = ax.contour(Y,z,data2-1000,[37.1],colors='gray',lw=10)
      ax.clabel(s, inline=1, fontsize=16,fmt='%4.2f', manual=[(500,-500)])
    else:
      print 'data2 will not be plotted!'
      #ax.plot(y,data2,'gray',lw=2);

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
icfile = 'IDEAL_IS_IC.nc'
exps = ['M2_exp0','M2_exp0','M2_exp13','M2_exp13'] 
labels = ['melting on','melting on','melting off','melting off']
abcd = ['a) ', 'b) ', 'c) ', 'd) ']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
lonh = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC.nc').variables['lonh'][:] # km
im = len(lonh)/2
lath = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC.nc').variables['lath'][:] # km
yh = netCDF4.Dataset(path+'M2_exp4/prog.nc').variables['yh'][:] # km
if os.path.isfile(path+'M2_exp4/IDEAL_IS_IC_1.nc'):
  print path+'M2_exp4/IDEAL_IS_1.nc'
  e   = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC_1.nc').variables['eta'][0,:,:,im]
else:
  e   = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC.nc').variables['eta'][0,:,:,im]
ssh = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im] #.mean(axis=1)
# This is the coordinate of the cell corners (u-points in 1D)
yq = netCDF4.Dataset(path+'M2_exp4/IDEAL_IS_IC.nc').variables['latq'][:]
yq = np.concatenate(([0],yq)) # Inserts left most edge of domain in to coordinate

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(12,10))

for j in range(len(exps)):
  path_to_file = path+exps[j]+'/ocean_month.nc'
  print 'path:',path_to_file
  print '--------------------------------- \n'
  # plot
  if j == 0:
    ax = ax1; var1='temp'; var2='ePBL_h_ML'; axm=1; hidex=True; hidey=False
  elif j == 1:
    ax = ax2; var1='salt'; var2='rhopot2'; axm=2; hidex=True; hidey=True
  elif j == 2:
    ax = ax3; var1='temp'; var2='ePBL_h_ML'; axm=1; hidex=False; hidey=False
  else:
    ax = ax4; var1='salt'; var2='rhopot2'; axm=2; hidex=False; hidey=True

  file_handle = netCDF4.Dataset(path_to_file)
  data1 = file_handle.variables[var1][-2::,:].mean(axis=0).mean(axis=2)  # ave last 2 years
  data2 = file_handle.variables[var2][-2::,:].mean(axis=0).mean(axis=axm) # ave last 2 years
  if j == 0 or j == 2 :
    cs1 = plot_section(file_handle, yq, ax, data1, -data2+ssh,im=im, hidex=hidex, hidey=hidey, m1=-2, m2=1.0)
  else:
    cs2 = plot_section(file_handle, yq, ax, data1, data2, im=im, cmap=cmaps.viridis,hidex=hidex, hidey=hidey, m1=33.85, m2=34.75)

  if j == 3: plt.tight_layout(w_pad=0.5, h_pad=0.5) ;f.subplots_adjust(bottom=0.15) 
  
  # aspect ratio
  ratio = 0.7
  xmin, xmax = ax.get_xlim()
  ymin, ymax = ax.get_ylim()
  ax.set_aspect(np.abs((xmax-xmin)/(ymax-ymin))*ratio, adjustable='box-forced')

  # text
  ax.text(12,-115,abcd[j], color='k', fontsize=18)
  ax.text(180,-1300,labels[j], color='k', fontsize=18)
  # fill IS and bottom
  yh[0]=0.
  val=-900
  e[-1,0:40]= val; e[0,0:40] = val; ssh[0:40] = val # 1km
  #e[-1,0:3]= val; e[0,0:3] = val; ssh[0:3] = val # 10km
  ax.fill_between(yh, e[0,:], 0.0, facecolor='white', interpolate=True)
  ax.fill_between(yh, e[-1,:], -4000, facecolor='lightgray', interpolate=True)
  ax.plot(yh, ssh,'k-',lw=1.5)

f.text(0.35,0.955,r'Downwelling (U$_{asf}$ = -5.0 m s$^{-1}$)', color='k', fontsize=24)
cbar_ax1 = f.add_axes([0.12, 0.07, 0.4, 0.03])
cbar1=f.colorbar(cs1, orientation='horizontal', cax=cbar_ax1)
cbar1.set_label(r'Potential temperature [$^o$C]', fontsize=18)

cbar_ax2 = f.add_axes([0.57, 0.07, 0.4, 0.03])
cbar2=f.colorbar(cs2, orientation='horizontal', cax=cbar_ax2, ticks=[33.8,34.0,34.2,34.4,34.6])
cbar2.set_label(r'Salinity', fontsize=18)
plt.savefig('sections_downwelling.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
