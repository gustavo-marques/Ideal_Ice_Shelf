from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import sys
sys.path.append('../../')
from misc import *

#plt.style.use('ggplot')
color1 = '#8B4513'
color2 = '#ff6347'
color3 = '#8470ff'
color4 = '#3cb371'

# plot some metrics for runs with varing wind forcing

path='/archive/gmm/Ideal_ice_shelf/Mode1/dx1km/Sigma_zstar/M1_exp13/'
files = ['out1/prog.nc','out1/prog.nc','out1/prog.nc','out1/prog.nc','out3/prog.nc','out3/prog.nc']
indices = [1,4,12,-1,365,-1]
titles = ['a)','b)','c)','d)', 'e)', 'f)']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
x = Dataset(path+'out1/ocean_geometry.nc').variables['geolon'][:]
y = Dataset(path+'out1/ocean_geometry.nc').variables['geolat'][:]
depth = Dataset(path+'out1/ocean_geometry.nc').variables['D'][:]

def get_data(exp,t):
    s=Dataset(exp)
    time=s.variables['time'][t] # in days
    print 'exp and time (hrs)',exp,time
    tr2=s.variables['tr1'][t,:,:,:]
    u=s.variables['uh'][t,:,:,:].sum(axis = 0)
    v=s.variables['vh'][t,:,:,:].sum(axis = 0)
    h=s.variables['h'][t,:,:,:]
    km, jm, im = h.shape
    # new arrays
    uh = np.zeros((jm,im)); vh = np.zeros((jm,im)); psi = np.zeros((jm,im))
    # u and v at h points
    utmp = 0.5 * (u[:,0:-1] + u[:,1::]) #u_i = 0.5(u_(i+0.5) + u_(i-0.5))
    vtmp = 0.5 * (v[0:-1,:] + v[1::,:]) #v_j = 0.5(v_(j+0.5) + v_(j-0.5))
    uh[:,1::] = utmp; uh[:,0] = 0.5*u[:,0] #u_i=1 = 0.5*u_(i=3/2)
    vh[1::,:] = vtmp; vh[0,:] = 0.5*v[0,:] #v_j=1 = 0.5*v_(j=3/2)
    # depth ave transport 
    uh = uh/(h.sum(axis=0)); vh = vh/(h.sum(axis=0)) # units m^2/s
    # mask and streamfunction
    uh = np.ma.masked_where(np.abs(uh) > 1.0e30, uh)
    vh = np.ma.masked_where(np.abs(vh) > 1.0e30, vh)
    psi[:,:] = (-uh.cumsum(axis=0) + vh.cumsum(axis=1))*0.5

    # depth ave tracer
    tr2_ave = (tr2 * h).sum(axis = 0)/h.sum(axis=0)
    #tr2_ave[tr2_ave<0.00001] = 0.00001
    return time,tr2_ave,uh,vh,psi

### Call the function make_cmap which returns your colormap
#colors = [(255,255,255), (0,0,255), (51,255,51), (255,255,51), (255,0,0), (255,0,255)]
colors = [(255,255,255), (0,0,255), (255,0,0)]
#colors = [(105,105,105),(119,136,153), (0,0,0)]
my_cmap = make_cmap(colors, bit=True)

i=0
fig, axes = plt.subplots(nrows=3, ncols=2,sharex='col', sharey='row',figsize=(12,12))
for ax in axes.flat:
    time,tracer,u,v,psi=get_data(path+files[i],indices[i])
    #cs=ax.pcolor(x,y,tracer,norm=matplotlib.colors.LogNorm(vmin=0.00001, vmax=5.0),cmap=plt.cm.gist_ncar_r)
    cs=ax.pcolor(x,y,tracer,norm=matplotlib.colors.LogNorm(vmin=0.000001, vmax=5.0),cmap=my_cmap)
    ax.contour(x,y,depth,[550,650,750,850,3900],colors='green')
    #ax.set_aspect('equal',adjustable = 'box-forced')
    ax.set_xlim([0,500])
    ax.set_ylim([0,1000])
    #ax.set_axis_bgcolor(color1)
    ax.set_axis_bgcolor('white')
    #cs.set_clim(dyelev.min(),dyelev.max())
    #psi_vals = np.arange(psi.min(),psi.max(),2)
    #strm = ax.streamplot(x, y, u, v, color=np.abs(psi),norm=plt.Normalize(0, 4.0e3),linewidth=1., density=[0.6, 1],arrowsize=7,cmap=plt.cm.copper)
    strm = ax.streamplot(x, y, u, v, color=np.abs(psi),norm=matplotlib.colors.LogNorm(vmin=1.0e-1, vmax=3.0e3),linewidth=1., density=[0.6, 1],arrowsize=7,cmap=plt.cm.gist_yarg)
    #ax.contour(x,y,psi,psi_vals,colors='k')
    ss = str("ax.set_title('%s time = %4.1f days')"% (titles[i],time))
    eval(ss)
    if i==4:
       ax.set_xlabel('x [km]', fontsize=18)
       ax.set_ylabel('y [km]', fontsize=18)

    i=i+1

fig.subplots_adjust(bottom=0.3)
cbar_ax = fig.add_axes([0.16, 0.2, 0.7, 0.02])
cbar=fig.colorbar(cs, orientation='horizontal', cax=cbar_ax, extend='min',ticks=[1.0e-6,1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1,1.0])
cbar.set_label(r'Depth-averaged passive tracer', fontsize=16)
cbar_ax = fig.add_axes([0.16, 0.11, 0.7, 0.02])
cbar=fig.colorbar(strm.lines, orientation='horizontal', cax=cbar_ax,ticks=[1.0e-1,1,10,1.0e2,1.0e3],extend='min')
cbar.set_label(r'Streamfunction magnitude [m$^2$ s$^{-1}$]', fontsize=16)
#fig.tight_layout()
plt.savefig('fig6.png',format='png',dpi=300,bbox_inches='tight')
plt.close()
print 'Done!'
