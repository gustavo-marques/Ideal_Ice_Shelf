from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#plt.style.use('ggplot')
import matplotlib
from matplotlib.colors import Normalize
matplotlib.rcParams.update({'font.size': 16})
import sys, os
sys.path.append('../../')
from remapping import mom_remapping
from misc import *
import gsw
import wright_eos as eos
import argparse

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Generate files for Idealized Ice Shelf problem.
      ''',
  epilog='Written by Gustavo Marques, Oct. 2016.')

  parser.add_argument('-dx', type=str, default='1',
      help='''Horizontal grid resolution (default = 1).''')

  parser.add_argument('-exp', type=str, default='M2_exp0',
      help='''Experiment name (default = M2_exp0).''')

  parser.add_argument('-save', type=str, default='T',
                help='''Save data into a .txt file (default = T).''')

  parser.add_argument('-out', type=str, default='out6',
      help='''Name of output file (default = out6).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp,area,y_loc):
    s=Dataset(exp+'/prog.nc')
    s1=Dataset(exp+'/ocean_sfc.nc')
    print 'Reading file', exp
    y = s.variables['yh'][:]
    y_tmp = np.nonzero(y<=y_loc)[0][-1]
    tm = len(s.variables['time'][:])
    tmp = 144
    if tm<tmp:
        tmp = tm
    Hout = np.zeros(tmp)
    Hin = np.zeros(tmp)
    Heat = np.zeros(tmp)
    delta_temp = np.zeros(tmp)
    for t in range(tmp):
      time0 = s.variables['time'][t-tmp-1]
      time = s.variables['time'][t-tmp]
      d_time = (time- time0)*3600*24 # in sec
      print 'delta_time',d_time
      print 'Time is (years):',time/365.
      heat_net=s1.variables['net_heat_surface'][t-tmp,:]
      vh = s.variables['vh'][t-tmp,:]
      temp0 = s.variables['temp'][t-tmp-1,:]
      temp = s.variables['temp'][t-tmp,:]
      h0 = s.variables['h'][t-tmp-1,:]
      h = s.variables['h'][t-tmp,:]
      salt0 = s.variables['salt'][t-tmp-1,:]
      salt = s.variables['salt'][t-tmp,:]
      e = s.variables['e'][t-tmp,:]
      depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
      Hout[t], Hin[t] = get_ht(temp,salt,depth,vh,y,y_loc)
      delta_temp[t] = get_temp_tendency(temp0,temp,salt0,salt,h0,h,depth,y_tmp,d_time)

      Heat[t] = (heat_net*area).sum()
    s.close()
    s1.close()
    return Heat.mean()/1.0e12, Hout/1.0e12, Hin/1.0e12 # in TW

def get_temp_tendency(t0,t1,s0,s1,h0,h1,d,yloc,dt):
         dxdy = 5000.0*5000.0
         cp = 3974.0 # heat capacity
         rho0 = 1028.0
         p = gsw.p_from_z(d[:,0:yloc,:],-72.0) * 1.0e4 # in Pascal
         t_freeze0 = eos.tfreeze(s0[:,0:yloc,:],p)
         t_freeze1 = eos.tfreeze(s1[:,0:yloc,:],p)
         dtemp0 = (t0[:,0:yloc,:] - t_freeze0)*h0[:,0:yloc,:]*dxdy
         dtemp1 = (t1[:,0:yloc,:] - t_freeze1)*h1[:,0:yloc,:]*dxdy
         temp_tend = ((dtemp1-dtemp0)*cp*rho0/dt).sum()
         print 'Temp. tendency:',temp_tend
         return temp_tend

def get_ht(t,s,d,vh,y,y_loc):
         '''
         Compute offshore heat transport, as defined in St-Laurent et al JPO 2012
         '''
         cp = 3974.0 # heat capacity
         rho0 = 1028.0
         #cp = 1; rho0 = 1
         tmp = np.nonzero(y<=y_loc)[0][-1]
         # transport at h point
         # mask transport. > 0.
         vhh = 0.5 * (vh[:,tmp-1,:] + vh[:,tmp,:])
         vhnew1 = np.ma.masked_where(vhh<0.0, vhh)
         vhnew2 = np.ma.masked_where(vhh>0.0, vhh)
         t_new = t[:,tmp,:]
         #t_new1 = np.ma.masked_where(vhh<0.0, t[:,tmp,:])
         #t_new2 = np.ma.masked_where(vhh>0.0, t[:,tmp,:])
         #s_new = np.ma.masked_where(vh[:,tmp,:]<0.0, s[:,tmp,:])
         s_new = s[:,tmp,:] # salinity at ice/ocean interface
         #d_new = np.ma.masked_where(vh[:,tmp,:]<0.0, d[:,tmp,:])
         d_new = d[0,tmp,:] # depth at ice/ocean interface
         p = gsw.p_from_z(d_new,-72.0) * 1.0e4 # in Pascal
         t_freeze = eos.tfreeze(s_new,p)
         dt = (t_new - t_freeze)
         print 'dt.min(),dt.max()', dt.min(),dt.max()
         dt = np.ma.masked_where(dt<0,dt)
         #dt2 = (t_new2 - t_freeze)
         Hout = (vhnew1*cp*rho0*dt) # watts
         Hin =  (vhnew2*cp*rho0*dt) # watts
         print 'Hout, Hin, difference:',Hout.sum(), Hin.sum(),Hout.sum()+Hin.sum()
         return Hout.sum(), Hin.sum()

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """

   #plt.style.use('ggplot')
   color1 = '#8B4513'
   color2 = '#ff6347'
   color3 = '#8470ff'
   color4 = '#3cb371'

   # plot some metrics for runs with varing wind forcing
   path='/work/gmm/Projects/Ideal_ice_shelf/Mode2/dx'+args.dx+'km/'+args.exp
   #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   if os.path.exists(path+'/out1/'):
     x = Dataset(path+'/out1/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/out1/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/out1/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/out1/MOM_Shelf_IC.nc').variables['shelf_area'][:]
     ssh = Dataset(path+'/out1/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/out1/ice_month.nc').variables['HI'][0,:]
   else:
     x = Dataset(path+'/ocean_geometry.nc').variables['geolon'][:]
     y = Dataset(path+'/ocean_geometry.nc').variables['geolat'][:]
     jm, im = x.shape
     depth = Dataset(path+'/ocean_geometry.nc').variables['D'][:]
     area = Dataset(path+'/MOM_Shelf_IC.nc').variables['shelf_area'][:]
     ssh = Dataset(path+'/IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,im/2]
     HI = Dataset(path+'/ice_month.nc').variables['HI'][0,:]

   e0 = Dataset(path+'/'+args.out+'/prog.nc').variables['e'][0:2,:,:,im/2].mean(axis=0)
   ISL = y[HI[:,0].mask == True,0][-1] - 10.0 # -10 km to make sure it is under cavity
   ISL = 180.
   print('Ice shelf lenght is (km):',ISL)

   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   heat,Hout,Hin=get_data(path+'/'+args.out,area,ISL)
   print 'net heat from coupler',heat
   print 'Hout',Hout.mean(),Hout.std()
   print 'Hin',Hin.mean(),Hin.std()
   if args.save=='T':
     print 'Saving data...'
     os.system('mkdir TXT')
     text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_heat_budget.txt', "w")
     text_file.write("%f \n" % (Hout.mean()))
     text_file.write("%f \n" % (Hout.std()))
     text_file.write("%f \n" % (Hin.mean()))
     text_file.write("%f \n" % (Hin.std()))
     text_file.write("%f \n" % (heat))
     text_file.close()

   print 'Done!'

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
