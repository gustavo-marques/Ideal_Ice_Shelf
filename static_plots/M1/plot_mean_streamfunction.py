from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import sys
sys.path.append('../../')
from remapping import mom_remapping
from misc import *
import wright_eos as eos

#plt.style.use('ggplot')
color1 = '#8B4513'
color2 = '#ff6347'
color3 = '#8470ff'
color4 = '#3cb371'

# plot some metrics for runs with varing wind forcing

path='/archive/gmm/Ideal_ice_shelf/Mode1/dx1km/Sigma_zstar/M1_exp5/'
files = ['out1/prog.nc','out1/prog.nc','out2/prog.nc','out3/prog.nc']
indices = [12,24,-1,-1]
titles = ['a)','b)','c)','d)']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
x = Dataset(path+'out1/ocean_geometry.nc').variables['geolon'][:]
y = Dataset(path+'out1/ocean_geometry.nc').variables['geolat'][:]
depth = Dataset(path+'out1/ocean_geometry.nc').variables['D'][:,250]
# remapping params
cs = mom_remapping.Remapping_Cs()
cs.remapping_scheme = 2
cs.degree = 2
z = np.linspace(0,depth.max(),2000)
h = 0.5* (z[0:-1]+z[1::])
# grid
Y,Z = np.meshgrid(y[:,0],h)

def get_data(cs,h,depth,z,exp,t):
    s=Dataset(exp)
    time=s.variables['time'][t]/365. # in yrs
    print 'exp and time (hrs)',exp,time
    rho=s.variables['rhopot2'][t,:,:,:].mean(axis = 2)
    v=s.variables['vh'][t,:,:,:].mean(axis = 2)
    h0=s.variables['h'][t,:,:,:].mean(axis = 2) 
    e0=s.variables['e'][t,:,:,:].mean(axis = 2) 
    km, jm = h0.shape
    # new arrays
    vh = np.zeros((km,jm)); psi = np.zeros((km,jm))
    # v at h points
    vtmp = 0.5 * (v[:,0:-1] + v[:,1::]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    vh[:,1::] = vtmp; vh[:,0] = 0.5*v[:,0] #v_j=1 = 0.5*v_(j=3/2)
    # mask and streamfunction
    
    vh = np.ma.masked_where(h0<0.001, vh)
    rho = np.ma.masked_where(h0<0.001, rho)
    psi[:,:] = (vh.cumsum(axis=0))*0.5
    #remap
    psi1=remap(cs,h0,e0,h,psi,z,depth)
    rho1=remap(cs,h0,e0,h,rho,z,depth)

    return time,rho1,psi1

### Call the function make_cmap which returns your colormap
#colors = [(255,255,255), (0,0,255), (255,0,0)]
colors = [(105,105,105),(119,136,153), (0,0,0)]
my_cmap = make_cmap(colors, bit=True)

i=0
fig, axes = plt.subplots(nrows=2, ncols=2,sharex='col', sharey='row',figsize=(12,12))
for ax in axes.flat:
    time,rho,psi=get_data(cs,h,depth,z,path+files[i],indices[i])
    cs=ax.contourf(Y,Z,psi,cmap=plt.cm.PuBu)
    ax.contour(x,y,psi,15,colors='gray')
    #ax.set_aspect('equal',adjustable = 'box-forced')
    ax.set_xlim([0,1000])
    ax.set_ylim([0,2000])
    ax.set_axis_bgcolor(color1)
    #cs.set_clim(dyelev.min(),dyelev.max())
    #psi_vals = np.arange(psi.min(),psi.max(),2)
    ax.contour(Y,Z,rho,[1037.05,1037.2],colors='k')
    ss = str("ax.set_title('%s time = %2.1f years')"% (titles[i],time))
    eval(ss)
    if i==2:
       ax.set_xlabel('y [km]', fontsize=18)
       ax.set_ylabel('depth [m]', fontsize=18)

    i=i+1

fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.16, 0.1, 0.7, 0.02])
cbar=fig.colorbar(cs, orientation='horizontal', cax=cbar_ax, extend='min',ticks=[1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,1.0,1.0e1,1.0e2])
cbar.set_label(r'Depth-averaged melt water tracer', fontsize=16)
plt.savefig('fig5.png',format='png',dpi=300,bbox_inches='tight')
plt.close()
print 'Done!'


def remap(cs,h0,e0,h,data0,z,depth):
   """
   remaps data...
   """
   km,im = h0.shape
   data1 = np.ma.zeros((len(h),im))
   for i in range(im):
     h1 = np.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
        h1[h<-e0[0,i]] = 0.0 # top
        h1[h>depth[i]] = 0.0 # bottom
        # need to account for SSH and make sure total thicknesses are
        # the same
        dh = h0[:,i].sum() - h1.sum()
        tmp1 = np.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
          if dh > 0.0:
             # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = np.ma.masked

     # for debugging
        #if h0[:,i].sum() != h1.sum():
        #   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
        data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
        # mask
        data1[h1==0.0,i] = np.ma.masked;
     else: # land/iceshelf
        data1[:,i] = np.ma.masked

   return data1
