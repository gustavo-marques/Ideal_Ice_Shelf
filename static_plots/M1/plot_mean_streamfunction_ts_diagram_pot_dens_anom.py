#!/usr/bin/env /home/gmm/Enthought/Canopy_64bit/User/bin/python2.7

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


exp=sys.argv[1]
print 'Processing experiment',exp + '...'

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
violet = '#EE82EE'
orange = '#FF8C00'
red = '#F08080'
grey = '#D3D3D3'
green = '#98FB98'
yellow = '#FFFFE0'
blue = '#6495ED'
navy = '#000080'
salmom = '#FA8072'
royalblue = '#4169E1'
stateblue = '#6A5ACD'
# plot some metrics for runs with varing wind forcing

path='/archive/gmm/Ideal_ice_shelf/Mode1/dx1km/Sigma_zstar/'+str(exp)+'/'
files = ['out3/ocean_month.nc','out1/prog.nc','out2/prog.nc','out3/prog.nc']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
x = Dataset(path+'out1/ocean_geometry.nc').variables['geolon'][:]
y = Dataset(path+'out1/ocean_geometry.nc').variables['geolat'][:]
depth = Dataset(path+'out1/ocean_geometry.nc').variables['D'][:,250]
ssh = Dataset(path+'out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,250]
e0 = Dataset(path+'out4/ocean_month.nc').variables['e'][:,:,:,250].mean(axis=0)
area = np.tile(Dataset(path+'out1/ocean_geometry.nc').variables['Ah'][:],(e0.shape[0]-1,1,1))
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

def get_data(cs,h,depth,z,exp,area):
    s=Dataset(exp)
    print 'Reading file', exp
    temp=s.variables['temp'][:,:,:,:].mean(axis=3).mean(axis=0)
    salt=s.variables['salt'][:,:,:,:].mean(axis=3).mean(axis=0)
    rho = eos.wright_eos(temp,salt,2.0e7)
    v=s.variables['vh'][:,:,:,:].sum(axis=3).mean(axis=0)
    h0=s.variables['h'][:,:,:,250].mean(axis=0)
    e0=s.variables['e'][:,:,:,250].mean(axis=0)
    volume = (s.variables['h'][:,:,:,:].mean(axis=0) * area).mean(axis=2)
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
    #psi[:,:] = -np.cumsum(vh[::-1,:],axis=0)[::-1,:]*0.5
    # top to bottom
    psi[:,:] = np.cumsum(vh,axis=0)*0.5
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
    salt1=remap(cs,h0,e0,h,salt,z,depth)
    temp1=remap(cs,h0,e0,h,temp,z,depth)
    # mask bad values
    psi1 = np.ma.masked_where(np.ma.getmask(rho1), psi1)

    return rho1,psi1/1.0e6,temp1,salt1,volume # in sv

### Call the function make_cmap which returns your colormap
colors = [(255,255,255), (255,255,153), (255,255,0)]
my_cmap = make_cmap(colors, bit=True)
# read data
psi_vals = np.linspace(-0.1,0.025,20)


norm = MidpointNormalize(midpoint=0)
rho,psi,temp,salt,volume=get_data(cs,h,depth,z,path+files[0],area)

### TS diagram ###
n1 = 40
n2 = 10
# Figure out boudaries (mins and maxs)
smin = 34.05  #np.min(salts) - (0.01 * np.min(salts))
smax = 34.85 #np.max(salts) + (0.01 * np.max(salts))
tmin = -2.3 #np.min(temps) - (0.1 * np.min(temps))
tmax = 1.5 #np.max(temps) + (0.1 * np.max(temps))
temps=np.linspace(tmin,tmax,10)
salts=np.linspace(smin,smax,10)
tf = eos.tfreeze(salts,0)
S,T = np.meshgrid(salts,temps)
sig2=eos.wright_eos(T,S,2e7)-1000

# separeta water masses
w1 = np.ma.masked_where(rho-1000<37.25,rho-1000)
w1_levs = np.linspace(37.25,37.3,10)

w2 = np.ma.masked_where(rho-1000<37.2,rho-1000)
w2 = np.ma.masked_where(w2>=37.25,w2)
w2_levs = np.linspace(37.2,37.245,10)

w3 = np.ma.masked_where(rho-1000<=37.05,rho-1000)
w3 = np.ma.masked_where(w3>=37.2,w3)
w3_levs = np.linspace(37.05,37.1999,10)

w4 = np.ma.masked_where(rho-1000>37.05,rho-1000)
w4_levs = np.linspace(36.68,37.05,10)

w5 = np.ma.masked_where(temp>tf.min(),rho-1000)
w5 = np.ma.masked_where(Y>900.,w5)

# new arrays
salt1 = salt[::n1,::n2]
temp1 = temp[::n1,::n2]
w11 = w1[::n1,::n2]
w22 = w2[::n1,::n2]
w33 = w3[::n1,::n2]
w44 = w4[::n1,::n2]
w55 = w5[::n1/5.,::n2/5.]
salt2 = salt[::n1/5.,::n2/5.]
temp2 = temp[::n1/5.,::n2/5.]

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
ax1.set_aspect('auto')
rho_vals = [36.4,36.6,36.8,37,37.2,37.4,37.6]
CS = plt.contour(S,T,sig2,rho_vals, colors='k',lw=0.25)
plt.clabel(CS, fontsize=18, inline=1, fmt='%4.2f') # Label every second level
sc1=ax1.scatter(salt1[w11.mask==False].flatten(),temp1[w11.mask==False].flatten(),s=80, alpha = 1, marker='.',c=green)
sc2=ax1.scatter(salt1[w22.mask==False].flatten(),temp1[w22.mask==False].flatten(),s=80, alpha = 1, marker='.',c=blue)
sc3=ax1.scatter(salt1[w33.mask==False].flatten(),temp1[w33.mask==False].flatten(),s=80, alpha = 1, marker='.',c=salmom)
sc4=ax1.scatter(salt1[w44.mask==False].flatten(),temp1[w44.mask==False].flatten(),s=80, alpha = 1, marker='.',c=stateblue)
sc4=ax1.scatter(salt2[w55.mask==False].flatten(),temp2[w55.mask==False].flatten(),s=80, alpha = 1, marker='.',c='cyan')
ax1.set_ylim(tmin,tmax)
ax1.set_xlim(smin,smax)
ax1.set_xlabel('Salinity', fontsize = 22)
ax1.set_ylabel(r'Temperature [$^o$C]', fontsize = 22)
ax1.plot(salts,tf,'gray')
ax1.text(34.25,-2.2,'ISW',color='cyan',fontsize = 22)
ax1.text(34.7,1.3,'CDW',color=salmom,fontsize = 22)
ax1.text(34.15,-1.55,'AASW',color=stateblue,fontsize = 22)
ax1.text(34.71,-1.,'DSW',color=green,fontsize = 22)
ax1.text(34.7,-0.6,'AABW',color=blue,fontsize = 22)
ax1.text(34.6,-1.8,r'$T_f$',color='k')
plt.tight_layout()
s = str("plt.savefig('%s-ts-diagram.png',format='png',dpi=300,bbox_inches='tight')"% (exp))
eval(s)
plt.show()



fig = plt.figure(facecolor='white',figsize=(10,8))
ax = fig.add_subplot(111,axisbg='black')
ct=ax.contourf(Y,-Z,psi,psi_vals,cmap=plt.cm.bwr, norm=norm)
plt.colorbar(ct, orientation='horizontal',ticks=[-0.1,-0.075,-0.05,-0.025,0.,0.025])
#cbar.set_label(r'Depth-averaged melt water tracer', fontsize=16)
ax.contour(Y,-Z,np.abs(psi),psi_vals,colors='gray',linewidth=0.2)
#ax.set_aspect('equal',adjustable = 'box-forced')
ax.set_xlim([0,1000])
ax.set_ylim([-4000,0])
#ct.set_clim(dyelev.min(),dyelev.max())
#psi_vals = np.arange(psi.min(),psi.max(),2)
s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',lw=0.5)
ax.clabel(s, inline=1, fontsize=14,fmt='%4.2f',)
# fill IS and bottom
ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
ax.plot(Y[0,:], ssh,'k-',lw=1.5)
ax.set_title(r'Mean overturning streamfunction $\psi_{mean}$ [sv]')
ax.set_xlabel('y [km]', fontsize=18)
ax.set_ylabel('depth [m]', fontsize=18)
s = str("plt.savefig('%s-psi.png',format='png',dpi=300,bbox_inches='tight')"% (exp))
eval(s)

# temp
fig = plt.figure(facecolor='white',figsize=(10,8))
ax = fig.add_subplot(111,axisbg='black')
ct=ax.contourf(Y,-Z,temp,np.linspace(-2.2,1.3,50))
plt.colorbar(ct, orientation='horizontal',ticks=[-2.,-1,0.,1])
ax.set_xlim([0,1000])
ax.set_ylim([-4000,0])
s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',linewidths=0.5,zorder=5)
ax.clabel(s, inline=2, fontsize=14,fmt='%4.2f')
# fill IS and bottom
ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
ax.plot(Y[0,:], ssh,'k-',lw=2)
ax.set_title(r'Temperature [$^o$C]')
ax.set_xlabel('y [km]', fontsize=18)
ax.set_ylabel('depth [m]', fontsize=18)
s = str("plt.savefig('%s-temp.png',format='png',dpi=300,bbox_inches='tight')"% (exp))
eval(s)

# salt
fig = plt.figure(facecolor='white',figsize=(10,8))
ax = fig.add_subplot(111,axisbg='black')
ct=ax.contourf(Y,-Z,salt,np.linspace(33.7,34.75,50))
plt.colorbar(ct, orientation='horizontal',ticks=[33.7,34,34.3,34.6])
ax.set_xlim([0,1000])
ax.set_ylim([-4000,0])
s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',linewidths=0.5,zorder=5)
ax.clabel(s, inline=2, fontsize=14,fmt='%4.2f')
# fill IS and bottom
ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
ax.plot(Y[0,:], ssh,'k-',lw=2)
ax.set_title(r'Salinity')
ax.set_xlabel('y [km]', fontsize=18)
ax.set_ylabel('depth [m]', fontsize=18)
s = str("plt.savefig('%s-salt.png',format='png',dpi=300,bbox_inches='tight')"% (exp))
eval(s)

# rho

fig = plt.figure(facecolor='white',figsize=(10,8))
ax = fig.add_subplot(111,axisbg='black')
s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',linewidths=0.5)
ax.clabel(s, inline=2, fontsize=14,fmt='%4.2f')
#ct=ax.contourf(Y,-Z,w1,np.linspace(36.7,37.3,50),)
#ct1=ax.contourf(Y,-Z,w1,w1_levs,cmap=plt.cm.Blues, zorder=1)
ct1=ax.contourf(Y,-Z,w1,w1_levs,colors=green, zorder=1)
ct2=ax.contourf(Y,-Z,w2,w2_levs,colors=blue, zorder=2)
ct3=ax.contourf(Y,-Z,w3,w3_levs,colors=salmom, zorder=3)
#ct4=ax.contourf(Y,-Z,w4,w4_levs, zorder=4)
ct4=ax.contourf(Y,-Z,w4,colors=stateblue, zorder=0)
ct5=ax.contourf(Y,-Z,w5, colors='cyan', zorder=5)
#plt.colorbar(ct, orientation='horizontal',ticks=[36.7,36.8,36.9,37,37.1,37.2,37.3])
#ax.contour(Y,-Z,np.abs(psi),8,colors='gray',linewidth=0.2)
#ax.set_aspect('equal',adjustable = 'box-forced')
ax.set_xlim([0,1000])
ax.set_ylim([-4000,0])
s=ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',linewidths=1)
ax.clabel(s, inline=2, fontsize=14,fmt='%4.2f',colors='k')
#ax.contour(Y,-Z,rho-1000,[37.05,37.20,37.25],colors='k',linewidths=0.5)
# fill IS and bottom
ax.fill_between(Y[0,:], e0[0,:], 0.0, facecolor='white', interpolate=True)
ax.fill_between(Y[0,:], e0[-1,:], -4000, facecolor='#A9A9A9', interpolate=True)
ax.plot(Y[0,:], ssh,'k-',lw=2)
ax.set_title(r'Potential density anomaly $\sigma_{\theta}$ [kg m$^{-3}$]')
ax.set_xlabel('y [km]', fontsize=18)
ax.set_ylabel('depth [m]', fontsize=18)
ax.text(100,-200,'ISW',color='cyan',fontsize = 22)
ax.text(850,-1800,'CDW',color='k',fontsize = 22)
ax.text(600,-180,'AASW',color='w',fontsize = 22)
ax.text(100,-700,'DSW',color='k',fontsize = 22)
ax.text(750,-3800,'AABW',color='k',fontsize = 22)

# Colorbars.
#dy = 0.08
#left, bottom, height, width = 0.18, 0.5, 0.25, 0.02
#rect1 = [left, bottom, height, width]  # Top.
#rect2 = [left, bottom - dy, height, width]  # Center.
#rect3 = [left, bottom - 2*dy, height, width]  # Bottom.
#rect4 = [left, bottom - 3*dy, height, width]  # Bottom.
#cax1 = fig.add_axes(rect1, axisbg='w')
#cax2 = fig.add_axes(rect2, axisbg='w')
#cax3 = plt.axes(rect3, axisbg='w')
#cax4 = plt.axes(rect4, axisbg='w')
#kw = dict(orientation='horizontal')
#cb1 = fig.colorbar(ct1, cax=cax1, ticks=[37.25,37.3],**kw); cb1.set_ticklabels([37.25,37.3])
#cb2 = fig.colorbar(ct2, cax=cax2, ticks=[37.2,37.24],**kw); cb2.set_ticklabels([37.2,37.24])
#cb3 = fig.colorbar(ct3, cax=cax3, ticks=[37.05,37.19],**kw); cb3.set_ticklabels([37.05,37.19])
#cb4 = fig.colorbar(ct4, cax=cax4, ticks=[36.68,37.04],**kw); cb4.set_ticklabels([36.7,37.04])
s = str("plt.savefig('%s-sigma-theta.png',format='png',dpi=300,bbox_inches='tight')"% (exp))
eval(s)

plt.show()
print 'Done!'

