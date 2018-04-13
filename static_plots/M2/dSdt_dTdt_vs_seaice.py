from matplotlib import pyplot as plt
from netCDF4 import Dataset, MFDataset
import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
import sys, os


def get_data(path):
    if os.path.isfile(path+'/ocean.stats.nc'):
        if os.path.isfile(path+'/out1/ocean.stats.nc'):
            s=MFDataset(path+'/out?/ocean.stats.nc')
        else:
            s=Dataset(path+'/ocean.stats.nc')
    else:
        s=MFDataset(path+'/out?/ocean.stats.nc')

    print 'Reading exp', path
    tmp = 365*2 # 2 years
    dS_dt=(s.variables['Salt'][-1]-s.variables['Salt'][-tmp])*0.5 # kg/year
    dT_dt=(s.variables['Heat'][-1]-s.variables['Heat'][-tmp])*0.5 # J/year
    s.close()
    return dS_dt, dT_dt

def driver():
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/'
   exps1 = ['M2_exp0','M2_exp1','M2_exp2','M2_exp3','M2_exp4'] # melt on
   exps2 = ['M2_exp13','M2_exp15','M2_exp16','M2_exp17','M2_exp14'] # melt off
   dx = ['1km','2km','5km','10km']
   symbs = ['o','x','^','+','*']
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #

   f, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(12,8))

   for j in range(len(dx)):
     for i in range(len(exps1)):
       path_to_file1 = path+'dx'+dx[j]+'/'+exps1[i]+'/'+exps1[i]+'_'+dx[j]+'_buoyancy.nc'
       path_to_file2 = path+'dx'+dx[j]+'/'+exps2[i]+'/'+exps2[i]+'_'+dx[j]+'_buoyancy.nc'
       print '--------------------------------- \n'
       seaice1 = Dataset(path_to_file1).variables['seaiceVolume'][:].mean()
       seaice2 = Dataset(path_to_file2).variables['seaiceVolume'][:].mean()
 
       # plot some metrics for runs with varing wind forcing
       path1=path+'dx'+dx[j]+'/'+exps1[i]+'/'
       path2=path+'dx'+dx[j]+'/'+exps2[i]+'/'
       # read *stas file data
       dS_dt1,dT_dt1=get_data(path1)
       dS_dt2,dT_dt2=get_data(path2)

       print 'Plotting data...'
       ax1.plot(dT_dt1, dS_dt1,symbs[j],color='blue')
       ax1.plot(dT_dt2, dS_dt2,symbs[j],color='red')
       ax1.set_xlabel('dT/dt [J/yr]'); ax1.set_ylabel('dS/dt [kg/yr]')
       ax2.plot(seaice1,dS_dt1, symbs[j],color='blue')
       ax2.plot(seaice2,dS_dt2, symbs[j],color='red')
       ax2.set_ylabel('dS/dt [kg/yr]'); ax2.set_xlabel('Sea ice volume [km^3]')

   plt.show()
   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': driver()
