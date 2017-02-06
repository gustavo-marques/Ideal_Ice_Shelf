#!/usr/bin/env python

# Compute various diagnostics for the Idealized Ice Shelf problem.
# Gustavo Marques

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import wright_eos as eos
import gsw 
import matplotlib.pyplot as plt
import warnings
import os

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

  parser.add_argument('-exp_name', type=str, default='test', help='''Name of the experiment (default is test).''')
  
  parser.add_argument('-prog_file', type=str, default='prog.nc', help='''Name of the prog. file (default is prog.nc).''')
  
  parser.add_argument('-ice_file', type=str, default='ice_month.nc', help='''Name of the ice file (default is ice_month.nc).''')
  
  parser.add_argument('-sfc_file', type=str, default='ocean_sfc.nc', help='''Name of the ocean surface fluxes file (default is ocean_sfc.nc).''')

  parser.add_argument('-ice_shelf_file', type=str, default='MOM_Shelf_IC.nc', help='''Name of the file that has the initial conditions for the ice shelf (default is MOM_Shelf_IC.nc).''')

  parser.add_argument('-energy', type=str, default='', help='''Extract the total energy (KE + PE) per unit mass (m2 s-2) from the given file and plot it as a function of time.''')
  
  parser.add_argument('-total_tke', type=str, default='', help='''Plot the total domain turbulent kinetic energy per unit mass (m2 s-2) from the given file.''')

  parser.add_argument('-cshelf_lenght', type=float, default=470.,
     help='''Continental shelf lenght in the y direction (km). Default is 470.''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   # load a few variables
   time = Dataset(args.prog_file).variables['time'][:]/365. # in years
   #time = Dataset(args.ice_file).variables['time'][:]/365. # in years
   x = Dataset(args.prog_file).variables['xh'][:] # in km
   y = Dataset(args.prog_file).variables['yh'][:]
   HI = Dataset(args.ice_file).variables['HI'][0,:]
   # ice shelf lenth
   args.ISL = y[HI[:,0].mask == True][-1]
   name = args.exp_name
   # create ncfile and zero fields. 
   create_ncfile(name,x,y,time,args)

   # lists used to save ncdata
   varname = []; var = []

   # create arrays
   salt_flux = np.zeros(len(time)); var.append('salt_flux'); varname.append('totalSaltFlux')
   FW_IS = np.zeros(len(time)); var.append('FW_IS'); varname.append('FwIS')
   FW = np.zeros(len(time)); var.append('FW'); varname.append('totalFw')
   polynya_area = np.zeros(len(time)); var.append('polynya_area'); varname.append('polynyaArea')
   ice_area = np.zeros(len(time)); var.append('ice_area'); varname.append('seaiceArea')
   ice_volume = np.zeros(len(time)); var.append('ice_volume'); varname.append('seaiceVolume')
   HI_max = np.zeros(len(time)); var.append('HI_max'); varname.append('maxSeaiceThick')
   AABW_transp = np.zeros(len(time)); var.append('AABW_transp'); varname.append('AABW')
   CDW_transp = np.zeros(len(time)); var.append('CDW_transp'); varname.append('CDW')
   NHT_shelf = np.zeros(len(time)); var.append('NHT_shelf'); varname.append('NorthwardTranspShelf')
   SHT_shelf = np.zeros(len(time)); var.append('SHT_shelf'); varname.append('SouthwardTranspShelf')
   NHT_ice_shelf = np.zeros(len(time)); var.append('NHT_ice_shelf'); varname.append('NorthwardTranspIceShelf')
   SHT_ice_shelf = np.zeros(len(time)); var.append('SHT_ice_shelf'); varname.append('SouthwardTranspIceShelf')
   oht1 = np.zeros(len(time)); var.append('oht1'); varname.append('OHT_shelf')
   oht2 = np.zeros(len(time)); var.append('oht2'); varname.append('OHT_ice_shelf')
   AABW_transp_x = np.zeros((len(time),len(x)))
   AABW_h = np.zeros((len(time),len(x)))
   MO_lenght = np.zeros((len(time),len(y),len(x)))
   var.append('MO_lenght'); varname.append('Monin_Obukhov_lenght')
   B0 = np.zeros((len(time),len(y),len(x)))
   var.append('B0'); varname.append('B0')
   B0_mean = np.zeros(len(time)); var.append('B0_mean'); varname.append('B0_mean')
   B0_shelf_mean = np.zeros(len(time)); var.append('B0_shelf_mean'); varname.append('B0_shelf_mean')
   B0_IS_mean = np.zeros(len(time))
   var.append('B0_IS_mean'); varname.append('B0_iceshelf_mean')

   # loop in time
   for t in range(len(time)):
           print 'Time (years):', time[t]
	   # load data
	   saltf = mask_bad_values(Dataset(args.ice_file).variables['SALTF'][t,:])
	   CI_tot = mask_bad_values(Dataset(args.ice_file).variables['CI_tot'][t,:])
	   HI = mask_bad_values(Dataset(args.ice_file).variables['HI'][t,:])
	   h = mask_bad_values(Dataset(args.prog_file).variables['h'][t,:])
	   vh = mask_bad_values(Dataset(args.prog_file).variables['vh'][t,:])
	   rhopot2 = mask_bad_values(Dataset(args.prog_file).variables['rhopot2'][t,:])
	   temp = mask_bad_values(Dataset(args.prog_file).variables['temp'][t,:])
	   salt = mask_bad_values(Dataset(args.prog_file).variables['salt'][t,:])
	   e = Dataset(args.prog_file).variables['e'][t,:]
	   mass_flux = mask_bad_values(Dataset(args.prog_file).variables['mass_flux'][t,:])
           # PRCmE: net surface water flux (lprec + melt + lrunoff - evap + calcing)
	   PRCmE = mask_bad_values(Dataset(args.sfc_file).variables['PRCmE'][t,:])
           depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
	   #tr2 = mask_bad_values(Dataset(args.prog_file).variables['tr2'][t,:])

           # diags functions
	   salt_flux[t] = get_saltf(x,y,saltf,CI_tot,args)
	   polynya_area[t] = get_polynya_area(x,y,CI_tot,args)
	   ice_area[t],ice_volume[t] = get_ice_diags(x,y,CI_tot,HI)
	   HI_max[t] = HI.max()
           AABW_transp[t],AABW_transp_x[t,:], AABW_h[t,:] = get_transport(x,y,vh,h,rhopot2,args)
           CDW_transp[t] = get_CDW(x,y,vh,salt,temp,args) 
           NHT_shelf[t] = get_total_transp(y,vh,args.cshelf_lenght,0) # northward
           SHT_shelf[t] = get_total_transp(y,vh,args.cshelf_lenght,1) # southward
           NHT_ice_shelf[t] = get_total_transp(y,vh,args.ISL,0) # northward
           SHT_ice_shelf[t] = get_total_transp(y,vh,args.ISL,1) # southward
           oht1[t] = get_oht(temp,salt,depth,vh,y,y_loc = 460.)
           oht2[t] = get_oht(temp,salt,depth,vh,y,y_loc = args.ISL)
           FW_IS[t], FW[t] = get_fw(PRCmE,mass_flux,x,y)

           #return l, B0, B0.mean(), B0_shelf, B0_IS
           MO_lenght[t,:], B0[t,:], B0_mean[t], B0_shelf_mean[t], B0_IS_mean[t] = compute_B0_MO_lenght(temp[0,:],salt[0,:],PRCmE,depth[0,:],t,y,args)


   if args.energy:
      print("Computing total energy...")
      compute_energy(args)

   if args.total_tke:
      print("Computing turbulent kinetic energy...")
      compute_total_tke(args)

   print 'Saving netcdf data...'

   for i in range(len(varname)):
       s = str("ncwrite(name,'%s',%s)"% (varname[i],var[i]))
       print s
       eval(s)
       #ncwrite(name,varname,var)

   print 'Done!'

def get_fw(PRCmE,mass_flux,x,y):
    '''
    Compute diags. related to freshwater/melt fluxes.
    '''
    area = np.ones(PRCmE.shape) * (x[1]-x[0]) * (y[1]-y[0]) * 1.0e3 # in m^2
    net = (area * PRCmE).sum()
    fw_is = mass_flux.sum()
    return fw_is, net

def compute_B0_MO_lenght(temp,salt,PRCmE,depth,t,y,args):
    '''
    Compute net surface buoyancy flux and Monin-Obukhov lenght scale
    '''
    tmp1 = np.nonzero(y<=args.ISL)[0][-1]
    tmp2 = np.nonzero(y<=args.cshelf_lenght)[0][-1]
    # constants
    rho_0 = 1028.0
    vonKar = 0.41
    g = 9.8
    Cp = 3974.0
    # get EoS coeffs
    p = gsw.p_from_z(depth,-70.0) * 1.0e4 # in Pascal
    beta = eos.beta_wright_eos(temp,salt,p)/rho_0 
    alpha = eos.alpha_wright_eos(temp,salt,p)/rho_0
    #print 'depth, beta, alpha',depth.min(), depth.max(), beta.min(), beta.max(), alpha.min(), alpha.max()
    # load local data
    ustar = mask_bad_values(Dataset(args.sfc_file).variables['ustar'][t,:])    
    sensible = mask_bad_values(Dataset(args.sfc_file).variables['sensible'][t,:])    
    latent = mask_bad_values(Dataset(args.sfc_file).variables['latent'][t,:])    
    shelf_area = Dataset(args.ice_shelf_file).variables['shelf_area'][0,:]
    # buoyancy flux
    B0 = -g * (alpha*(-(sensible + latent)/(rho_0 * Cp)) - beta*(PRCmE*salt/rho_0))
    B0_shelf = B0[tmp1:tmp2,:].mean()
    B0_IS = B0[0:tmp1,:].mean()
    # Monin-Obukhov Length
    l = ustar**3/(vonKar * B0)
    # mask values outside cavity
    l[l==0.0] = -1e+34
#    l = np.ma.masked_where(depth == depth.max(), l)
#    print 'Monin-Obukhov Length min/max',l.min(), l.max()

    return l, B0, B0.mean(), B0_shelf, B0_IS

def ncwrite(name,var_name,data):
    '''
    Write data (with var_name) into netCDF file (name.nc).
    '''
    file = Dataset(name+'.nc','r+')
    file.variables[var_name][:]=data[:]
    file.close()
    print('*** Variable '+var_name+' was written successfully. \n')
    return

def mask_bad_values(data,value = -1.30856803e+26):
         return np.ma.masked_where(data == value, data)

def get_oht(t,s,d,vh,y,y_loc):
         '''
         Compute the onshore heat transport, as defined in St-Laurent et al JPO 2012
         '''
         cp = 3974.0 # heat capacity
         rho0 = 1028.0
         tmp = np.nonzero(y<=y_loc)[0][-1]
         # mask transport. > 0.
         vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])
         t_new = np.ma.masked_where(vh[:,tmp,:]>0.0, t[:,tmp,:])
         s_new = np.ma.masked_where(vh[:,tmp,:]>0.0, s[:,tmp,:])
         d_new = np.ma.masked_where(vh[:,tmp,:]>0.0, d[:,tmp,:])
         p = gsw.p_from_z(d_new,-70.0) * 1.0e4 # in Pascal
         t_freeze = eos.tfreeze(s_new,p)
         dt = (t_new - t_freeze)
         #print 't-tf min/max',dt.min(),dt.max()
         
         oht = - (vhnew*cp*rho0*dt).sum() # watts
         return oht

def get_total_transp(y,vh,loc_y,opt):
         '''
         Compute the total cross-shelf vol. transport at the shelf break.
         '''
         tmp = np.nonzero(y<=loc_y)[0][-1]
         if opt == 0: # mask southward flow
            vhnew = np.ma.masked_where(vh[:,tmp,:]<0.0, vh[:,tmp,:]) 
         else: # mask northward flow
            vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])

         return vhnew.sum()/1.0e6 # in sv

def get_CDW(x,y,vh,s,t,args):
         '''
         Compute the onshore volume transport of CDW at the shelf break.
         '''
         t_cdw = 0.0
         s_cdw = 34.6
         tmp = np.nonzero(y<=args.cshelf_lenght)[0][-1]
         vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])         
         tnew = np.ma.masked_where(t[:,tmp,:]>=t_cdw, t[:,tmp,:])
         snew = np.ma.masked_where(s[:,tmp,:]>=s_cdw, s[:,tmp,:])
         transp = (vhnew * tnew * snew / (tnew * snew)).sum() # in m^3/s
        
         print 'CDW Transport (sv)', transp/1.0e6
 
         return transp/1.0e6 # in sv

def get_transport(x,y,vh,h,rhopot2,args):
	 '''
         Compute the total volume transport, outflow thickness.
	 '''
         T = -0.75; S = 34.7
         rho_min = eos.wright_eos(T,S,2.0e7)-1000.
         rho_min = 37.18
         rhopot2 = rhopot2 - 1000.0
         print 'rho_AABW, rhopot2 min/max',rho_min, rhopot2.min(), rhopot2.max()
	 #dx = np.ones((h[:,0,:].shape)) * (x[1]-x[0]) * 1.0e3 # in m
         # end of cont. shelf
	 tmp = np.nonzero(y<=args.cshelf_lenght)[0][-1]
	 # mask h below 1.0e-12 m
	 hnew = np.ma.masked_where(h[:,tmp,:]<=1.0e-12, h[:,tmp,:])
         # mask vel. <= 0.
	 vhnew = np.ma.masked_where(vh[:,tmp,:]<=0.0, vh[:,tmp,:])
         # mask rhopot2 < rho_min
         sig2 = np.ma.masked_where(rhopot2[:,tmp,:]<rho_min, rhopot2[:,tmp,:])

         total_transp = (vhnew * sig2 / sig2).sum() # in m^3/s
         transp_x = (vhnew * sig2 / sig2).sum() # in m^3/s
         thickness = (hnew * sig2 *vhnew / (vhnew*sig2)).sum(axis=0) # in m^3/s
         # set masked values to zero
         if total_transp is np.ma.masked: total_transp = 0.0

         print 'AABW Transport (sv)', total_transp/1.0e6, '\n'
 
         return total_transp/1.0e6, transp_x, thickness # in sv

def get_saltf(x,y,saltf,CI_tot,args):
	'''
        Compute the total salt flux into (+) the ocean in the cont. shelf 
	region (cshelf_lenght). Salt flux is just computed when CI_tot < 0.7. 
	'''
	# get indices for region of interest
        tmp = np.nonzero(y<=args.cshelf_lenght)[0]
	# mask points where CI_tot >= 0.7
	CI_tot = np.ma.masked_where(CI_tot >= 0.7, CI_tot)
	saltf = saltf * (CI_tot/CI_tot)
        tot_saltf = saltf[tmp,:].sum()
        # set masked values to zero
        if tot_saltf is np.ma.masked: tot_saltf = 0.0

	return tot_saltf

def get_polynya_area(x,y,CI_tot,args):
	'''
        Compute the total (coastal, i.e., for y <= ISL + 100 km) polynya area.
	'''
	grid_area = np.ones((CI_tot.shape)) * (x[1]-x[0]) * (y[1]-y[0])
        # get indices for region of interest
        tmp = np.nonzero(y<=args.ISL + 100.0)[0]
        # mask points where CI_tot >= 0.7
        CI_tot = np.ma.masked_where(CI_tot >= 0.7, CI_tot)
        grid_area = grid_area * (CI_tot/CI_tot)
        poly_area = grid_area[tmp,:].sum()
        # set masked values to zero
        if poly_area is np.ma.masked: poly_area = 0.0

	return poly_area

def get_ice_diags(x,y,CI_tot,HI):
	'''
        Compute total area covered by sea ice. Sea ice is present when CI_tot >= 0.7.
	'''
        grid_area = np.ones((CI_tot.shape)) * (x[1]-x[0]) * (y[1]-y[0])
	CI_tot = np.ma.masked_where(CI_tot < 0.85, CI_tot)
        ice_area = (grid_area * (CI_tot/CI_tot)).sum()
        ice_volume = (grid_area * 1.0e6 * HI).sum()
        # set masked values to zero
        if ice_area is np.ma.masked: ice_area = 0.0

	return ice_area,ice_volume

def compute_total_tke(args):
    file = args.total_tke
    time = Dataset(file).variables['time'][:]
    TKEu = np.zeros(len(time))
    TKEv = np.zeros(len(time))
    ubar = Dataset(file).variables['u'][0] * 0.; ubar2=0*ubar;
    vbar = Dataset(file).variables['v'][0] * 0.; vbar2=0*vbar;
    # using the following identity
    # <u'^2> = <(u - ubar)^2>, where <> = domain ave; u' = u-ubar
    #        = <u^2> - <ubar^2>
    for t in range(len(time)):
        print 'Time is (days):', time[t]
        u = Dataset(file).variables['u'][t,:]
        v = Dataset(file).variables['v'][t,:]
        ubar = ubar + u
        ubar2 = ubar2 + u*u
        vbar = vbar + v
        vbar2 = vbar2 + v*v
   
    ubar, vbar = ubar/len(time), vbar/len(time)
    print 'Finish computing mean, now tke...'
    for t in range(len(time)):
        print 'Time is (days):', time[t]
        u = Dataset(file).variables['u'][t,:]
        v = Dataset(file).variables['v'][t,:]
        up = (u - ubar)**2; vp = (v - vbar)**2
        TKEu[t] = up.sum() 
        TKEv[t] = vp.sum() 

    if args.savefig:
       # plot
       plt.figure()
       plt.plot(time/365.,TKEu,'r',lw=2.5)
       plt.plot(time/365.,TKEv,'b',lw=2.5)
       plt.xlabel('Time [years]')
       plt.ylabel('Total TKE (red = u, blue = v) [m2 s-2]')
       plt.grid()
       plt.savefig(args.exp_name+'_total_tke.png')

def compute_energy(args):
    file = args.energy
    os.system("awk '/MOM Day/{print $3}' " + file + " > tmp0" ) 
    os.system("awk '/MOM Day/{print $6}' " + file + " > tmp1" ) 
    os.system("awk -F',' '{print $1}' tmp1 > tmp2")
    time = np.loadtxt('tmp0')/365. # in yr.
    energy = np.loadtxt('tmp2')
    os.system('rm tmp?')
    if args.savefig:
       plt.figure()
       plt.plot(time,energy,lw=2.5)
       plt.xlabel('Time [years]')
       plt.ylabel('Total energy [m2 s-2]')
       plt.grid()
       plt.show()

def create_ncfile(exp_name, xx, yy, ocean_time, args): # may add exp_type
   """
   Creates a netcdf file with the fields required for the isomip+ experiments. Different fields are generated based on the type of experiment that is being analyzed (Ocean0, Ocean1 etc).
   """
   # open a new netCDF file for writing.
   ncfile = Dataset(exp_name+'.nc','w',format='NETCDF4')
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
   time.description = 'time since start of simulation'
   time.long_name = 'time since start of simulation'

   AABW = ncfile.createVariable('AABW',np.dtype('float32').char,('time'))
   AABW.units = 'sv'; AABW.description = 'Offshore transport of AABW computed at the shelf break'

   CDW = ncfile.createVariable('CDW',np.dtype('float32').char,('time'))
   CDW.units = 'sv'; CDW.description = 'Onshore transport of CDW computed at the shelf break'

   FwIS = ncfile.createVariable('FwIS',np.dtype('float32').char,('time'))
   FwIS.units = 'kg/s'; FwIS.description = 'Total mass flux of freshwater across the ice-ocean interface.'
   totalFw = ncfile.createVariable('totalFw',np.dtype('float32').char,('time'))
   totalFw.units = 'kg/s'; totalFw.description = 'Total net mass flux of freshwater across the entire domain (lprec+melt).'

   NorthwardTranspShelf = ncfile.createVariable('NorthwardTranspShelf',np.dtype('float32').char,('time'))
   NorthwardTranspShelf.units = 'sv' 
   NorthwardTranspShelf.description = 'Northward cross-shelf volume transport computed at the shelf break'

   SouthwardTranspShelf = ncfile.createVariable('SouthwardTranspShelf',np.dtype('float32').char,('time'))
   SouthwardTranspShelf.units = 'sv'
   SouthwardTranspShelf.description = 'Southward cross-shelf volume transport computed at the shelf break'

   NorthwardTranspIceShelf = ncfile.createVariable('NorthwardTranspIceShelf',np.dtype('float32').char,('time'))
   NorthwardTranspIceShelf.units = 'sv'
   NorthwardTranspIceShelf.description = 'Northward cross-shelf volume transport computed at the ice shelf front'

   SouthwardTranspIceShelf = ncfile.createVariable('SouthwardTranspIceShelf',np.dtype('float32').char,('time'))
   SouthwardTranspIceShelf.units = 'sv'
   SouthwardTranspIceShelf.description = 'Southward cross-shelf volume transport computed at the ice shelf front '

   OHT_shelf = ncfile.createVariable('OHT_shelf',np.dtype('float32').char,('time'))
   OHT_shelf.units = 'Watts'; OHT_shelf.description = 'Onshore heat transport at shelf break'

   OHT_ice_shelf = ncfile.createVariable('OHT_ice_shelf',np.dtype('float32').char,('time'))
   OHT_ice_shelf.units = 'Watts'; OHT_ice_shelf.description = 'Onshore heat transport at the ice shelf edge'

   totalSaltFlux = ncfile.createVariable('totalSaltFlux',np.dtype('float32').char,('time'))
   totalSaltFlux.units = 'kg/(m^2 s)'; totalSaltFlux.description = 'total net salt flux into the ocean from coupler'

   polynyaArea = ncfile.createVariable('polynyaArea',np.dtype('float32').char,('time'))
   polynyaArea.units = 'km^2'; polynyaArea.description = 'total polynya (open water) area'

   seaiceArea = ncfile.createVariable('seaiceArea',np.dtype('float32').char,('time'))
   seaiceArea.units = 'km^2'; seaiceArea.description = 'total area covered by sea ice'

   seaiceVolume = ncfile.createVariable('seaiceVolume',np.dtype('float32').char,('time'))
   seaiceVolume.units = 'm^3'; seaiceArea.description = 'total volume of sea ice'

   maxSeaiceThick = ncfile.createVariable('maxSeaiceThick',np.dtype('float32').char,('time'))
   maxSeaiceThick.units = 'm'; maxSeaiceThick.description = 'maximum sea ice thickness'

   Monin_Obukhov_lenght = ncfile.createVariable('Monin_Obukhov_lenght',np.dtype('float32').char,('time','ny','nx'))
   Monin_Obukhov_lenght.units = 'm'
   Monin_Obukhov_lenght.description = 'Monin-Obukhov lenght scale'

   B0 = ncfile.createVariable('B0',np.dtype('float32').char,('time','ny','nx'))
   B0.units = 'm2/s3'
   B0.description = 'Net surface buoyancy flux'

   B0_mean = ncfile.createVariable('B0_mean',np.dtype('float32').char,('time'))
   B0_mean.units = 'm2/s3'
   B0_mean.description = 'Domain ave. net surface buoyancy flux'

   B0_shelf_mean = ncfile.createVariable('B0_shelf_mean',np.dtype('float32').char,('time'))
   B0_shelf_mean.units = 'm2/s3'
   B0_shelf_mean.description = 'Net surface buoyancy flux ave. within the shelf region'

   B0_iceshelf_mean = ncfile.createVariable('B0_iceshelf_mean',np.dtype('float32').char,('time'))
   B0_iceshelf_mean.units = 'm2/s3'
   B0_iceshelf_mean.description = 'Net surface buoyancy flux ave. within the ice shelf region'

      
   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:] 
   time[:] = ocean_time[:]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
