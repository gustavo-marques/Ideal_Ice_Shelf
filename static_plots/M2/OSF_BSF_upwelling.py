from matplotlib import pyplot as plt
import netCDF4
import numpy as np
from misc import *
import wright_eos as eos
from computeOSF import computeOSF
from remapping import mom_remapping
import sys,os
sys.path.append('/work/gmm/pyGVtools/')
sys.path.append('colormap/')
import colormaps as cmaps
import m6toolbox
#plt.style.use('ggplot')
import matplotlib

matplotlib.rcParams.update({'font.size': 18})

def remap(cs,h0,e0,h,data0,z,depth):
   """
   remaps data0 into a regular grid
   """
   km,im = h0.shape
   data1 = numpy.ma.zeros((len(h),im))
   for i in range(im):
     h1 = numpy.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
        h1[h<-e0[0,i]] = 0.0 # top
        h1[h>depth[i]] = 0.0 # bottom
        # need to account for SSH and make sure total thicknesses are
        # the same
        dh = h0[:,i].sum() - h1.sum()
        tmp1 = numpy.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
          if dh > 0.0:
             # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = numpy.ma.masked

     # for debugging
        #if h0[:,i].sum() != h1.sum():
        #   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
        data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
        # mask
        data1[h1==0.0,i] = numpy.ma.masked;
     else: # land/iceshelf
        data1[:,i] = numpy.ma.masked

   return data1


# Define a function to plot a section
def plot_OSF(file_handle, h1, ax, psi, data2,im=250, eta='e',yvar='yh', rep='pcm', xlim=(0,650), ylim=(-2000,0), cmap=plt.cm.bwr, hidex=False, hidey=False, m1=-0.15, m2=0.15):
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
    z = 0.5*(e[0:-1,:]+e[1:,:])
    [Y,TMP] = np.meshgrid(y,z[:,0])

    km, jm = e.shape

    # interpolate psi (can't remap it, it has been integrated already)
    psi = np.zeros((len(h1),jm))
    for j in range(jm):
      if h0[:,j].sum() > 0.001:
         tmp = -0.5*(e0[0:-1,j]+e0[1::,j])
         psi[:,j] = np.interp(h1, tmp, osf[:,j])
      else:
         psi[:,j] = -1.0e32

    # mask bad values
    psi = np.ma.masked_where(np.ma.getmask(rho1), psi)
    psi = np.ma.masked_where(psi==-1.0e32, psi)

    cs = ax.pcolormesh(1e-3*YOut, ZOut, 1e-6*osf,vmin=m1,vmax=m2,cmap=cmap)
    s=ax.contour(Y,z,data2-1000,[37.12],colors='k',lw=10)
    ax.clabel(s, inline=1, fontsize=16,fmt='%4.2f', manual=[(500,-500)])
    
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
exps = ['M2_exp4','M2_exp4','M2_exp14','M2_exp14'] 
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

# remapping params
cs = mom_remapping.Remapping_Cs()
cs.remapping_scheme = 2
cs.degree = 2
z = numpy.linspace(0,4000.,500)
h1 = 0.5* (z[0:-1]+z[1::])
# grid
Y,Z = numpy.meshgrid(y,h1)

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
    # OSF
    cs1 = plot_OSF(file_handle, yq, ax, data1, -data2+ssh,im=im, hidex=hidex, hidey=hidey, m1=-2, m2=1.0)
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



f.text(0.35,0.955,r'Upwelling (U$_{asf}$ = 5.0 m s$^{-1}$)', color='k', fontsize=24)
cbar_ax1 = f.add_axes([0.12, 0.07, 0.4, 0.03])
cbar1=f.colorbar(cs1, orientation='horizontal', cax=cbar_ax1)
cbar1.set_label(r'Potential temperature [$^o$C]', fontsize=18)

cbar_ax2 = f.add_axes([0.57, 0.07, 0.4, 0.03])
cbar2=f.colorbar(cs2, orientation='horizontal', cax=cbar_ax2, ticks=[33.8,34.0,34.2,34.4,34.6])
cbar2.set_label(r'Salinity', fontsize=18)
plt.savefig('OSF_BSF_upwelling.png',format='png',dpi=300,bbox_inches='tight')
plt.show()
