from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
from matplotlib.colors import Normalize
matplotlib.rcParams.update({'font.size': 16})
import sys
sys.path.append('../../')
from remapping import mom_remapping
from misc import *
import wright_eos as eos

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

#plt.style.use('ggplot')
color1 = '#8B4513'
color2 = '#ff6347'
color3 = '#8470ff'
color4 = '#3cb371'

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx1km/Sigma_zstar/M1_exp11/'
files = ['out3/ocean_month.nc','out1/prog.nc','out2/prog.nc','out3/prog.nc']
indices = [12,24,-12,-12]
titles = ['a)','b)','c)','d)']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
x = Dataset(path+'out1/ocean_geometry.nc').variables['geolon'][:]
y = Dataset(path+'out1/ocean_geometry.nc').variables['geolat'][:]
depth = Dataset(path+'out1/ocean_geometry.nc').variables['D'][:,250]#.mean(axis=1)
# remapping params
cs = mom_remapping.Remapping_Cs()
cs.remapping_scheme = 2
cs.degree = 2
z = np.linspace(0,depth.max(),2000)
h = 0.5* (z[0:-1]+z[1::])
# grid
Y,Z = np.meshgrid(y[:,0],h)

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def get_data(cs,h,depth,z,exp,t):
    s=Dataset(exp)
    time=s.variables['time'][t]/365. # in yrs
    print 'exp and time (hrs)',exp,time
    temp=s.variables['temp'][:,:,:,250].mean(axis=3).mean(axis=0)
    salt=s.variables['salt'][:,:,:,250].mean(axis=3).mean(axis=0)
    rho = eos.wright_eos(temp,salt,2.0e7)
    v=s.variables['vh'][:,:,:,:].sum(axis=3).mean(axis=0)
    h0=s.variables['h'][:,:,:,250].mean(axis=0)
    e0=s.variables['e'][:,:,:,250].mean(axis=0)
    s.close()

    km, jm = h0.shape
    # new arrays
    vh = np.zeros((km,jm)); psi = np.zeros((km,jm))
    # v at h points
    vtmp = 0.5 * (v[:,0:-1] + v[:,1::]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    vh[:,1::] = vtmp; vh[:,0] = 0.5*v[:,0] #v_j=1 = 0.5*v_(j=3/2)
    # mask and streamfunction

    #vh = np.ma.masked_where(h0<0.001, vh)
    rho = np.ma.masked_where(h0<0.001, rho)
    # integrate from bottom to top
    psi[:,:] = -np.cumsum(vh[::-1,:],axis=0)[::-1,:]*0.5
    # interpolate psi (can't remap it, it has been integrated already)
    psi1 = np.zeros((len(h),jm))
    for j in range(jm):
       if h0[:,j].sum() > 0.001:
          tmp = -0.5*(e0[0:-1,j]+e0[1::,j])
          psi1[:,j] = np.interp(h, tmp, psi[:,j])

       else:
          psi1[:,j] = -1.0e32

    #remap rho
    rho1=remap(cs,h0,e0,h,rho,z,depth)
    # mask bad values
    psi1 = np.ma.masked_where(np.ma.getmask(rho1), psi1)

    return time,rho1,psi1/1.0e6 # in sv

### Call the function make_cmap which returns your colormap
colors = [(0,0,255), (255,255,255), (255,0,0)]
my_cmap = make_cmap(colors, bit=True)
psi_vals = np.linspace(-0.2,0.025,10)
norm = MidpointNormalize(midpoint=0)
time,rho,psi=get_data(cs,h,depth,z,path+files[i],indices[i])
fig = plt.figure(facecolor='white')
ax = fig.add_subplot(111,axisbg='black')
ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.coolwarm, norm=norm)
plt.colorbar(ct, orientation='horizontal',ticks=[-0.2,-0.15,-0.1,-0.05,0.,0.025])
#cbar.set_label(r'Depth-averaged melt water tracer', fontsize=16)
ax.contour(Y,-Z,np.abs(psi),5,colors='gray',linewidth=0.5)
#ax.set_aspect('equal',adjustable = 'box-forced')
ax.set_xlim([0,1000])
ax.set_ylim([-4000,0])
#ct.set_clim(dyelev.min(),dyelev.max())
#psi_vals = np.arange(psi.min(),psi.max(),2)
ax.contour(Y,-Z,rho,[1037.05,1037.2],colors='k')
ss = str("ax.set_title('%s time = %2.1f years')"% (titles[0],time))
eval(ss)
ax.set_xlabel('y [km]', fontsize=18)
ax.set_ylabel('depth [m]', fontsize=18)
plt.savefig('fig5.png',format='png',dpi=300,bbox_inches='tight')
plt.close()
print 'Done!'

