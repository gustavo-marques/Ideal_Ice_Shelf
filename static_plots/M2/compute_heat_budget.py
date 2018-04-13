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

  parser.add_argument('-out', type=str, default='',
      help='''Name of output file (default = '').''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)


def get_data(exp,area,y_loc):
    cp = 3974.0 # heat capacity
    rho0 = 1028.0
    s=Dataset(exp+'/prog.nc')
    s1=Dataset(exp+'/ocean_sfc.nc')
    print 'Reading file', exp
    y = s.variables['yh'][:]
    dy = (y[1]-y[0])*1.0e3
    area_cell = dy*dy
    y_tmp = np.nonzero(y<=y_loc)[0][-1]
    tm = len(s.variables['time'][:])
    tmp = 144
    if tm<tmp:
        tmp = tm
    Hout = np.zeros(tmp-1)
    Hin = np.zeros(tmp-1)
    frazil = np.zeros(tmp-1)
    sensible = np.zeros(tmp-1)
    hnet = np.zeros(tmp-1)
    #heat_PmE = np.zeros(tmp-1)
    heat_melt = np.zeros(tmp-1)
    dtemp_dt_cav = np.zeros(tmp-1)
    dtemp_dt_off = np.zeros(tmp-1)
    tt = np.zeros(tmp-1)
    for t in range(1,tmp):
      time = s.variables['time'][t-tmp]
      time0 = s.variables['time'][t-tmp-1]
      dt = (time-time0)*3600*24. # in sec
      print 'Time is (years):',time/365.
      tt[t-1] = time/365.
      if 'frazil' in s1.variables:
        frz=s1.variables['frazil'][tm-tmp+t,:]
      else:
        frz=np.zeros(sen.shape)

      sen=s1.variables['sensible'][t-tmp,:]
      net=s1.variables['net_heat_surface'][t-tmp,:]
      m1=s1.variables['lprec'][t-tmp,:] #/(3600*24.*365) # in m/s
      #heatPmE=s1.variables['Heat_PmE'][t-tmp,:]
      #heatPmE=s1.variables['heat_content_surfwater'][t-tmp,:]
      vh = s.variables['vh'][t-tmp,:]
      temp0 = s.variables['temp'][t-tmp-1,:] + 273.0 # in K
      temp = s.variables['temp'][t-tmp,:] + 273.0 # in K
      h0 = s.variables['h'][t-tmp-1,:]; vol0=h0*area_cell
      h = s.variables['h'][t-tmp,:]; vol = h*area_cell
      
      dtemp_dt = (temp-temp0)*vol*rho0*cp/dt
      #dtemp_dt = (temp-temp0)/dt
      dtemp_dt_cav[t-1] = (dtemp_dt[:,0:y_tmp,:]).sum() # K/s
      dtemp_dt_off[t-1] = (dtemp_dt[:,y_tmp::,:]).sum() 
      h = s.variables['h'][t-tmp,:]
      salt = s.variables['salt'][t-tmp,:]
      e = s.variables['e'][t-tmp,:]
      depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
      Hout[t-1], Hin[t-1] = get_ht(temp,salt,depth,vh,y,y_loc,cp,rho0)
      sensible[t-1] = (sen*area).sum()
      frazil[t-1] = (frz*area).sum()
      hnet[t-1] = (net*area).sum()
      #heat_PmE[t-1] = (heatPmE*area).sum()
      heat_melt[t-1] = (m1*3.34E+05*area).sum()
    s.close()
    s1.close()
    return hnet/1.0e12, frazil/1.0e12, sensible/1.0e12, Hout/1.0e12, Hin/1.0e12, dtemp_dt_cav/1.0e12, dtemp_dt_off/1.0e12, heat_melt/1.0e12, tt # H* in TW

def get_ht(t,s,d,vh,y,y_loc,cp,rho0):
         '''
         Compute offshore heat transport, as defined in St-Laurent et al JPO 2012
         '''
         #cp = 1; rho0 = 1
         tmp = np.nonzero(y<=y_loc)[0][-1]
         # transport at h point
         # mask transport. > 0.
         vhh = 0.5 * (vh[:,tmp-1,:] + vh[:,tmp,:])
         vhnew1 = np.ma.masked_where(vhh<0.0, vhh)
         vhnew2 = np.ma.masked_where(vhh>0.0, vhh)
         t_new = t[:,tmp,:] # temp at cavity exit
         s_new = s[0,tmp,:] # salinity at cavity exit
         d_new = d[0,tmp,:] # depth at cavity exit
         p = gsw.p_from_z(d_new,-72.0) * 1.0e4 # in Pascal
         t_freeze = eos.tfreeze(s_new,p) + 273.0 # in K
         #t_freeze = eos.tfreeze(34.,p) + 273.0 # in K
         print 't_freeze min/max',t_freeze.min(), t_freeze.max()
         dt = (t_new - t_freeze)
         print 'dt max/min', dt.max(), dt.min()
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
   ISL = y[HI[:,0].mask == True,0][-1] - 20.0 # -20 km to make sure it is under cavity
   ISL = 190.
   print('Ice shelf lenght is (km):',ISL)

   ### Call the function make_cmap which returns your colormap
   colors = [(0,0,255), (255,255,255), (255,0,0)]
   my_cmap = make_cmap(colors, bit=True)
   # read data
   Hnet,frazil,sensible,Hout,Hin,dTdt_cav,dTdt_off,heat_melt,time=get_data(path+'/'+args.out,area,ISL)
   print 'heat_melt',heat_melt.mean(), heat_melt.std()
   print 'frazil',frazil.mean(), frazil.std()
   print 'sensible',sensible.mean(), sensible.std()
   #print 'Heat_PmE',Heat_PmE.mean(), Heat_PmE.std()
   print 'Hout',Hout.mean(),Hout.std()
   print 'Hin',Hin.mean(),Hin.std()
   print 'Hnet',Hnet.mean(),Hnet.std()
   print 'Fraz_Lat',frazil.mean() + sensible.mean() #+ Heat_PmE.mean()
   print 'Mass_in',Hin.mean()+Hout.mean()-sensible.mean()-frazil.mean()
   print 'dTdt cavity',dTdt_cav.mean(), dTdt_cav.std()
   print 'dTdt offshore',dTdt_off.mean(), dTdt_off.std()

   if args.save=='T':
     print 'Saving data...'
     # NC
     os.system('mkdir ncfiles')
     exp_name = args.exp+'_dx'+args.dx
     create_ncfile(exp_name, x[0,:], y[:,0], time)
     ncwrite(exp_name,'Frazil',frazil)
     #ncwrite(exp_name,'Heat_PmE',Heat_PmE)
     ncwrite(exp_name,'Lat',sensible)
     ncwrite(exp_name,'Hin',Hin)
     ncwrite(exp_name,'Hout',Hout)
     ncwrite(exp_name,'HNET',Hnet)
     ncwrite(exp_name,'dTdt',dTdt_cav)
     #ncwrite(exp_name,'Fraz_Lat',frazil+sensible+Heat_PmE)
     ncwrite(exp_name,'Fraz_Lat',frazil+sensible)
     ncwrite(exp_name,'Mass_in',-frazil-sensible+Hin+Hout)
     # TXT 
     os.system('mkdir TXT')
     text_file = open('TXT/'+args.exp+'_dx'+args.dx+'_heat_budget.txt', "w")
     text_file.write("%f \n" % (Hout.mean()))
     text_file.write("%f \n" % (Hout.std()))
     text_file.write("%f \n" % (Hin.mean()))
     text_file.write("%f \n" % (Hin.std()))
     text_file.write("%f \n" % (dTdt_cav.mean()))
     text_file.write("%f \n" % (dTdt_cav.std()))
     text_file.write("%f \n" % (dTdt_off.mean()))
     text_file.write("%f \n" % (dTdt_off.std()))
     text_file.close()

   print 'Done!'

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset('ncfiles/'+name+'_heat_budget_cavity.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def create_ncfile(exp_name, xx, yy, ocean_time): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset('ncfiles/'+exp_name+'_heat_budget_cavity.nc','w',format='NETCDF4')
   # dimensions
   nx = len(xx) ; ny = len(yy)
   # create dimensions.
   ncfile.createDimension('time', len(ocean_time))
   ncfile.createDimension('nx',nx)
   ncfile.createDimension('ny',ny)

   # create variables, assign units and provide decription
   x = ncfile.createVariable('x',np.dtype('float32').char,('nx'))
   x.units = 'km'
   x.description = 'x location of cell centers'
   x.long_name = 'x location of cell centers'
   y = ncfile.createVariable('y',np.dtype('float32').char,('ny'))
   y.units = 'km'
   y.description = 'y location of cell centers'
   y.long_name = 'y location of cell centers'

   time = ncfile.createVariable('time',np.dtype('float32').char,('time'))
   time.units = 'years'
   time.description = 'last two years of simulation'
   time.long_name = 'last two years of simulation'

   Hin = ncfile.createVariable('Hin',np.dtype('float32').char,('time'))
   Hin.units = 'TW'; Hin.description = 'Heat flux entering the cavity'

   dTdt = ncfile.createVariable('dTdt',np.dtype('float32').char,('time'))
   dTdt.units = 'TW'; dTdt.description = 'Temporal change in temperature inside cavity'

   Hout = ncfile.createVariable('Hout',np.dtype('float32').char,('time'))
   Hout.units = 'TW'; Hout.description = 'Heat flux leaving the cavity'

   Lat = ncfile.createVariable('Lat',np.dtype('float32').char,('time'))
   Lat.units = 'TW'; Lat.description = 'Total latent heat flux inside cavity'

   Frazil = ncfile.createVariable('Frazil',np.dtype('float32').char,('time'))
   Frazil.units = 'TW'; Frazil.description = 'Total frazil heat flux inside cavity'

   Fraz_Lat = ncfile.createVariable('Fraz_Lat',np.dtype('float32').char,('time'))
   Fraz_Lat.units = 'TW'; Fraz_Lat.description = 'Net heat flux in the cavity from frazil+latent+Heat_PmE'

   HNET = ncfile.createVariable('HNET',np.dtype('float32').char,('time'))
   HNET.units = 'TW'; HNET.description = 'Net heat flux in the cavity from coupler'

   Heat_PmE = ncfile.createVariable('Heat_PmE',np.dtype('float32').char,('time'))
   Heat_PmE.units = 'TW'; Heat_PmE.description = 'Net heat flux in the cavity from mass flux into the ocean'
    
   Mass_in = ncfile.createVariable('Mass_in',np.dtype('float32').char,('time'))
   Mass_in.units = 'TW'; Mass_in.description = 'Net heat flux in the cavity from water crossing ocean surface (relative to freezing point)'
   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]-ocean_time[0]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'_heat_budget_cavity.nc!')
   return
# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
