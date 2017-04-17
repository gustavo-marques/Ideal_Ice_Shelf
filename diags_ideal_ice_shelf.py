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

  parser.add_argument('-n', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-prog_file', type=str, default='prog.nc', help='''Name of the prog. file (default is prog.nc).''')

  parser.add_argument('-ice_file', type=str, default='ice_month.nc', help='''Name of the ice file (default is ice_month.nc).''')

  parser.add_argument('-sfc_file', type=str, default='ocean_sfc.nc', help='''Name of the ocean surface fluxes file (default is ocean_sfc.nc).''')

  parser.add_argument('-ice_shelf_file', type=str, default='MOM_Shelf_IC.nc', help='''Name of the file that has the initial conditions for the ice shelf (default is MOM_Shelf_IC.nc).''')

  parser.add_argument('-cshelf_lenght', type=float, default=470.,
     help='''Continental shelf lenght in the y direction (km). Default is 470.''')

  parser.add_argument('-t0', type=int, default=0,
     help='''Initial time indice to start computations. Default is 0.''')

  parser.add_argument('-tf', type=int, default=-1,
     help='''Final time indice for computations. Default is -1 (all starting from t0).''')

  parser.add_argument('-AABW_rho', type=float, default=1037.15,
     help='''Minimum density (sigma2) to compute AABW transport (kg/m3). Default is 1037.15''')

  parser.add_argument('-CDW_rho', type=float, default=1037.05,
     help='''Minimum density (sigma2) to compute CDW transport (kg/m3). Default is 1037.05''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
   """
   This is where all the action happends, i.e., calls to different functions that generate
   files that are needed.
   """
   # load a few variables
   if args.tf == -1:
      time = Dataset(args.prog_file).variables['time'][args.t0::]/365. # in years
   else:
      time = Dataset(args.prog_file).variables['time'][args.t0:args.tf+1]/365.

   x = Dataset(args.prog_file).variables['xh'][:] # in km
   y = Dataset(args.prog_file).variables['yh'][:]
   HI = Dataset(args.ice_file).variables['HI'][0,:]
   shelf_area = Dataset(args.ice_shelf_file).variables['shelf_area'][0,:,:]
   # ice shelf lenth
   args.ISL = y[HI[:,0].mask == True][-1] - 10.0 # -10 km to make sure it is under cavity
   print('Ice shelf lenght is (km):',args.ISL)
   name = args.n
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
   CDW1_transp = np.zeros(len(time)); var.append('CDW1_transp'); varname.append('CDW1')
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
   melt = np.zeros(len(time))
   var.append('melt'); varname.append('Melt')
   total_mass_flux = np.zeros(len(time))
   var.append('total_mass_flux'); varname.append('TotalMassFlux')
   dyn_pump = np.zeros(len(time))
   #var.append('dyn_pump'); varname.append('DynamicalPump')

   # loop in time
   for t in range(args.t0,len(time)+args.t0):
           tt = t - args.t0 # time indice used in the arrays
           print 'Time (years):', time[tt]
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
	   melt_all = mask_bad_values(Dataset(args.prog_file).variables['melt'][t,:])
           # PRCmE: net surface water flux (lprec + melt + lrunoff - evap + calcing)
	   PRCmE = mask_bad_values(Dataset(args.sfc_file).variables['PRCmE'][t,:])
           depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
	   #tr1 = mask_bad_values(Dataset(args.prog_file).variables['tr1'][t,:])
	   #tr2 = mask_bad_values(Dataset(args.prog_file).variables['tr2'][t,:])
           pressure = gsw.p_from_z(depth,-75.) * 1.0e4 # in Pa [1 db = 10e4 Pa]
           rho = eos.wright_eos(temp,salt,pressure)

           # diags functions
	   salt_flux[tt] = get_saltf(x,y,saltf,CI_tot,args)
	   polynya_area[tt] = get_polynya_area(x,y,CI_tot,args)
	   ice_area[tt],ice_volume[tt] = get_ice_diags(x,y,CI_tot,HI)
	   HI_max[tt] = HI.max()
           AABW_transp[tt],AABW_transp_x[tt,:], AABW_h[tt,:] = get_transport(x,y,vh,h,rhopot2,args)
           CDW_transp[tt] = get_CDW(x,y,vh,rhopot2,args.cshelf_lenght)
           CDW1_transp[tt] = get_CDW(x,y,vh,rhopot2,300.0)
           NHT_shelf[tt] = get_total_transp(y,vh,args.cshelf_lenght,0) # northward
           SHT_shelf[tt] = get_total_transp(y,vh,args.cshelf_lenght,1) # southward
           NHT_ice_shelf[tt] = get_total_transp(y,vh,args.ISL,0) # northward
           SHT_ice_shelf[tt] = get_total_transp(y,vh,args.ISL,1) # southward
           oht1[tt] = get_oht(temp,salt,depth,vh,y,y_loc = 460.)
           oht2[tt] = get_oht(temp,salt,depth,vh,y,y_loc = args.ISL)
           FW_IS[tt], FW[tt] = get_fw(PRCmE,mass_flux,x,y)
           melt[tt] = get_melt(melt_all,shelf_area)
           total_mass_flux[tt] = get_total_mass_flux(mass_flux,shelf_area)
           #return l, B0, B0.mean(), B0_shelf, B0_IS
           #MO_lenght[tt,:], B0[tt,:], B0_mean[tt], B0_shelf_mean[tt], B0_IS_mean[tt] = compute_B0_MO_lenght(temp[0,:],salt[0,:],PRCmE,depth[0,:],t,y,args)
           #dyn_pump[tt] = get_dyn_pump(y,vh,tr1,tr2,rho,total_mass_flux[tt],np.abs(SHT_ice_shelf[tt]),args.ISL)
           print '\n'

   print 'Saving netcdf data...'

   for i in range(len(varname)):
       s = str("ncwrite(name,'%s',%s)"% (varname[i],var[i]))
       print s
       eval(s)
       #ncwrite(name,varname,var)

   print 'Done!'

def get_dyn_pump(y,vh,dye1,dye2,rho,mass_flux,transp_in,loc_y):
   '''
   Compute the dynamical pump efficiency following Jourdan et al 2017 JGR
   '''
   tmp = np.nonzero(y<=loc_y)[0][-1]
   #vh_in = np.ma.masked_where(dye1[:,tmp,:]<0.5, vh[:,tmp,:]) # CAUTION, they are in diff points
                                                              # this needs to be changed
   vh_in = np.ma.masked_where(dye2[:,tmp,:]>1.0, vh[:,tmp,:])
   vh_in = np.ma.masked_where(vh[:,tmp,:]>0.0, vh_in[:,:])

   vh_out = np.ma.masked_where(vh[:,tmp,:]<0.0, vh[:,tmp,:])
   vh_out = np.ma.masked_where(dye2[:,tmp,:]<5., vh_out[:,:])
   vh_out = np.ma.masked_where(dye1[:,tmp,:]>0.5, vh_out[:,:])

   rho_melt = mass_flux/(vh_out.sum()-np.abs(vh_in.sum()))
   print 'rho_melt,vh_in.sum(),vh_out.sum(),mass_flux',rho_melt,vh_in.sum(),vh_out.sum(),mass_flux
   dyn_pump = np.abs(vh_in.sum())*rho_melt/mass_flux
   print 'dyn_pump',dyn_pump

   return dyn_pump

def get_melt(melt_all,shelf_area):
    '''
    Compute the spatial mean ice shelf melting/freezing.
    '''
    melt_all = np.ma.masked_where(shelf_area == 0, melt_all)

    return melt_all.mean()

def get_total_mass_flux(mass_flux,shelf_area):
    '''
    Compute the total mass flux of melt water.
    '''
    mass_flux = np.ma.masked_where(shelf_area == 0, mass_flux)
    mass_flux = np.ma.masked_where(mass_flux == 0, mass_flux)

    return mass_flux.sum()

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

def get_CDW(x,y,vh,rhopot2,yloc):
         '''
         Compute the onshore volume transport of CDW at the shelf break.
         '''
         rhopot2 = rhopot2 - 1000.0
         rho_max = args.AABW_rho - 1000.
         rho_min = args.CDW_rho - 1000.
         tmp = np.nonzero(y<=yloc)[0][-1]
         vhnew = np.ma.masked_where(vh[:,tmp,:]>0.0, vh[:,tmp,:])
         sig2_min = np.ma.masked_where(rhopot2[:,tmp,:]<rho_min, rhopot2[:,tmp,:])
         sig2_max = np.ma.masked_where(rhopot2[:,tmp,:]>rho_max, rhopot2[:,tmp,:])
         transp = (vhnew * sig2_min * sig2_max / (sig2_min * sig2_max)).sum() # in m^3/s

         print 'CDW Transport (sv)', transp/1.0e6

         return transp/1.0e6 # in sv

def get_transport(x,y,vh,h,rhopot2,args):
	 '''
         Compute the total volume transport, outflow thickness.
	 '''
         T = -0.75; S = 34.7
         #rho_min = eos.wright_eos(T,S,2.0e7)-1000.
         rho_min = args.AABW_rho - 1000.
         rhopot2 = rhopot2 - 1000.0
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

         print 'AABW Transport (sv)', total_transp/1.0e6

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

   CDW1 = ncfile.createVariable('CDW1',np.dtype('float32').char,('time'))
   CDW1.units = 'sv'; CDW1.description = 'Onshore transport of CDW computed at y = 300 km'

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

   Melt = ncfile.createVariable('Melt',np.dtype('float32').char,('time'))
   Melt.units = 'm/year'
   Melt.description = 'Domain ave. melt rate'

   TotalMassFlux = ncfile.createVariable('TotalMassFlux',np.dtype('float32').char,('time'))
   TotalMassFlux.units = 'kg/s'
   TotalMassFlux.description = 'Total mass flux of melt water'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
