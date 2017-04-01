from matplotlib import pyplot as plt
import netCDF4
import numpy as np
#plt.style.use('ggplot')

# plot some metrics for runs with varing wind forcing

path='/lustre/f1/unswept/Gustavo.Marques/MOM6-examples/ice_ocean_SIS2/IDEAL_IS/dx2km/Sigma_zstar/'
fnames = ['M2_exp14_2km','M2_exp13_2km','M2_exp12_2km','M2_exp2_2km']
exps = ['M2_exp14/','M2_exp13/','M2_exp12/','M2_exp2/']
param = 'cavity_size'
labels = ['1.0','1.4','1.8','2.2']
cc = ['green','red','blue','black']
vars = ['OHT_ice_shelf','OHT_shelf','seaiceVolume','Melt']
units = ['OHT(y=200 km) [TW]','OHT(y=460 km) [TW]','Total Sea ice volume [m$^3$]','<$\dot{m}$>  [m yr$^{-1}$]']
c = ['g','r','b','k']
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

for n in range(len(vars)):
   var = vars[n]
   plt.figure()
   mean_all = np.zeros(len(fnames))
   std_all = np.zeros(len(fnames))
   for i in range(len(fnames)):
       print 'path',path+exps[i]+fnames[i]+'.nc'
       if n <=1:
          data = netCDF4.Dataset(path+exps[i]+fnames[i]+'.nc').variables[var][:]/1.0e12 # TW
       else:
          data = netCDF4.Dataset(path+exps[i]+fnames[i]+'.nc').variables[var][:]
       #units = netCDF4.Dataset(path+fnames[i]+'/'+fnames[i]+'.nc').variables[var].units
       # plot
       mean_all[i]=data.mean()
       std_all[i]=data.std()
       plt.errorbar(i, mean_all[i], yerr=std_all[i], fmt='o',color=cc[i],ms=8)

   plt.plot(range(len(fnames)), mean_all,'k-')

   plt.xticks(range(len(fnames)), labels, rotation='horizontal')
   # Pad margins so that markers don't get clipped by the axes
   plt.margins(0.1)
   # Tweak spacing to prevent clipping of tick-labels
   plt.subplots_adjust(bottom=0.15)
   plt.grid()
   #plt.title(var, fontsize=18)
   plt.xlabel(r'Cavity volume [$\times$ $10^{4}$ km$^3$]', fontsize=16)
   plt.ylabel(r''+units[n], fontsize=16)
   plt.savefig(param+'_sensitivity_'+vars[n]+'.png',format='png',dpi=300,bbox_inches='tight')
   plt.show()

# get cavity volume and plot shape
icname = 'IDEAL_IS_IC'
print path+exps[0]+icname+'.nc'
x = netCDF4.Dataset(path+exps[0]+icname+'.nc').variables['lonh'][:]
y = netCDF4.Dataset(path+exps[0]+icname+'.nc').variables['lath'][:]
[XX,YY] = np.meshgrid(x,y)

cell_area = (x[1]-x[0])*1.0e3 * (y[1]-y[0])*1.0e3

fig, ax = plt.subplots()
for i in range(len(exps)):
    h = netCDF4.Dataset(path+exps[i]+icname+'.nc').variables['h'][0,:,0:100,:]
    print 'Volume is:',((h*cell_area).sum()) * 1.0e-9
    depth = netCDF4.Dataset(path+exps[i]+'ocean_geometry.nc').variables['D'][:]
    if i ==len(exps)-1:
       ax.contour(XX,YY,depth,[600,1000,2000,3000],colors = c[i],lw=1.5)

    tmp = np.zeros(len(x))
    for l in range(len(x)):
        tmp[l] = y[np.nonzero(depth[:,l] > 0.5)[0][0]]

    ax.plot(x,tmp,c[i])

ax.set_xlabel('x [km]',fontsize=15)
ax.set_ylabel('y [km]',fontsize=15)
plt.savefig('cavities_shapes.png')
plt.show()

