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
from get_ML_state import ml_average

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

  parser.add_argument('-AABW_rho', type=float, default=1037.2,
                help='''Minimum density (sigma2) to compute AABW mass flux (kg/m^3). Default is 1037.2''')

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
   depth = Dataset('ocean_geometry.nc').variables['D'][:]
   area = Dataset('ocean_geometry.nc').variables['Ah'][:]
   HI = Dataset(args.ice_file).variables['HI'][0,:]
   area_ice_shelf = Dataset(args.ice_shelf_file).variables['shelf_area'][0,:,:]
   area_cshelf = np.ma.masked_where(area_ice_shelf>0,area)
   tmp1 = np.nonzero(y<=args.cshelf_lenght)[0][-1]
   area_cshelf[tmp1::,:] = 0.0
   area_cshelf = np.ma.masked_where(area_cshelf==0.,area_cshelf)
   # ice shelf lenth
   args.ISL = y[HI[:,0].mask == True][-1] - 10.0 # -10 km to make sure it is under cavity
   print('Ice shelf lenght is (km):',args.ISL)
   name = args.n
   # create ncfile and zero fields.
   create_ncfile(name,x,y,time,args)

   # lists used to save ncdata
   varname = []; var = []

   # create arrays
   Bseaice = np.zeros((len(time)))
   var.append('Bseaice'); varname.append('Bseaice')
   Biceshelf = np.zeros((len(time)))
   var.append('Biceshelf'); varname.append('Biceshelf')
   Nshelf = np.zeros((len(time)))
   var.append('Nshelf'); varname.append('Nshelf')
   Tshelf = np.zeros((len(time)))
   var.append('Tshelf'); varname.append('Tshelf')
   Sshelf = np.zeros((len(time)))
   var.append('Sshelf'); varname.append('Sshelf')
   Ustar_shelf = np.zeros((len(time)))
   var.append('Ustar_shelf'); varname.append('Ustar_shelf')
   qflux = np.zeros((len(time)))
   var.append('qflux'); varname.append('qflux')
   AABW_transp = np.zeros(len(time)); var.append('AABW_transp'); varname.append('AABW')
   sflux = np.zeros(len(time)); var.append('sflux'); varname.append('sflux')
   hml_shelf = np.zeros(len(time)); var.append('hml_shelf'); varname.append('hml_shelf')
   ice_volume = np.zeros(len(time)); var.append('ice_volume'); varname.append('seaiceVolume')
   # loop in time
   for t in range(args.t0,len(time)+args.t0):
           tt = t - args.t0 # time indice used in the arrays
           print 'Time (years):', time[tt]
	   # load data
	   saltf = -mask_bad_values(Dataset(args.ice_file).variables['SALTF'][t,:])
           Qall = mask_bad_values(Dataset(args.sfc_file).variables['net_heat_surface'][t,:])
           lprec = mask_bad_values(Dataset(args.sfc_file).variables['lprec'][t,:])
           ustar = mask_bad_values(Dataset(args.sfc_file).variables['ustar'][t,:])
           vh = mask_bad_values(Dataset(args.prog_file).variables['vh'][t,:])
	   CI_tot = mask_bad_values(Dataset(args.ice_file).variables['CI_tot'][t,:])
	   HI = mask_bad_values(Dataset(args.ice_file).variables['HI'][t,:])
	   h = mask_bad_values(Dataset(args.prog_file).variables['h'][t,:])
	   hml = mask_bad_values(Dataset(args.sfc_file).variables['ePBL_h_ML'][t,:])
	   temp = mask_bad_values(Dataset(args.prog_file).variables['temp'][t,:])
	   salt = mask_bad_values(Dataset(args.prog_file).variables['salt'][t,:])
	   e = Dataset(args.prog_file).variables['e'][t,:]
	   #mass_flux = mask_bad_values(Dataset(args.prog_file).variables['mass_flux'][t,:])
	   #melt_all = mask_bad_values(Dataset(args.prog_file).variables['melt'][t,:])
           # PRCmE: net surface water flux (lprec + melt + lrunoff - evap + calcing)
           depth = 0.5*(e[0:-1,:,:]+e[1::,:,:]) # vertical pos. of cell
	   #tr1 = mask_bad_values(Dataset(args.prog_file).variables['tr1'][t,:])
	   #tr2 = mask_bad_values(Dataset(args.prog_file).variables['tr2'][t,:])
           pressure = gsw.p_from_z(depth,-75.) * 1.0e4 # in Pa [1 db = 10e4 Pa]
           rho = eos.wright_eos(temp,salt,pressure)
           # get ML properties
           SST,SSS=ml_average(temp,salt,hml,h)
           SST = np.ma.masked_where(np.abs(SST)>1.0e15,SST)
           SST = np.ma.masked_where(hml<1.0e-2,SST)
           SSS = np.ma.masked_where(np.abs(SSS)>1.0e15,SSS)
           SSS = np.ma.masked_where(hml<1.0e-2,SSS)
           # diags functions
	   sflux[tt] = (area_cshelf *saltf).sum()
           print 'Total salt flux on the shelf (kg/(s m^2)):',sflux[tt]
           hml_shelf[tt] = (np.ma.masked_where(area_cshelf.mask == True, hml)).mean()
           print 'Mixed layer depth ave. over the cont. shelf ():',hml_shelf[tt]

           Bseaice[tt] = buoyancy_flux_cshelf(Qall,saltf,lprec,area_cshelf,hml,temp,salt,h,SST,SSS)
           print 'B sea ice (m4/s3):',Bseaice[tt]

           Biceshelf[tt] = buoyancy_flux_ishelf(Qall,lprec,area_ice_shelf,hml,temp,salt,h,SST,SSS)
           print 'B shelf (m4/s3):',Biceshelf[tt]

           qflux[tt] = (np.ma.masked_where(area_cshelf.mask == True, Qall)).sum()
           print 'Q flux shelf (W/m2):',qflux[tt]
           Ustar_shelf[tt] = (np.ma.masked_where(area_cshelf.mask == True, ustar)).mean()
           print 'Ustar shelf (m/s):',Ustar_shelf[tt]
           h1 = -np.diff(depth,axis=0)
           # ambient strat.
           N = np.sqrt((9.8/1028.0)*np.diff(rho,axis=0)/h1)
           Ntmp = ((N*h1).sum(axis=0)/h.sum(axis=0))
           Nshelf[tt] = (np.ma.masked_where(area_cshelf.mask == True,Ntmp)).mean()
           print 'N shelf (1/s):',Nshelf[tt]
           # salt
           Stmp = ((salt*h).sum(axis=0)/h.sum(axis=0))
           Sshelf[tt] = (np.ma.masked_where(area_cshelf.mask == True,Stmp)).mean()
           print 'Salt shelf (psu):',Sshelf[tt]
           # temp
           Ttmp = ((temp*h).sum(axis=0)/h.sum(axis=0))
           Tshelf[tt] = (np.ma.masked_where(area_cshelf.mask == True,Ttmp)).mean()
           print 'Temp shelf (C):',Tshelf[tt]
           # AABW mass flux
           AABW_transp[tt] = get_mass_flux(x,y,vh,h,rho,args)
           # sea ice vol
           tmp,ice_volume[tt] = get_ice_diags(x,y,CI_tot,HI)
           print 'Total sea ice vol. (m^3):',ice_volume[tt]
           print '\n'

   print 'Saving netcdf data...'

   for i in range(len(varname)):
       s = str("ncwrite(name,'%s',%s)"% (varname[i],var[i]))
       print s
       eval(s)
       #ncwrite(name,varname,var)

   print 'Done!'

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


def get_mass_flux(x,y,vh,h,rhopot2,args):
         '''
         Compute the total volume transport, outflow thickness.
         '''
         rho_min = args.AABW_rho - 1000.
         rho = rhopot2 - 1000.0
         # end of cont. shelf
         tmp = np.nonzero(y<=args.cshelf_lenght)[0][-1]
         # mask h below 1.0e-12 m
         hnew = np.ma.masked_where(h[:,tmp,:]<=1.0e-12, h[:,tmp,:])
         # mask vel. <= 0.
         vhnew = np.ma.masked_where(vh[:,tmp,:]<=0.0, vh[:,tmp,:])
         # mask rhopot2 < rho_min
         sig2 = np.ma.masked_where(rho[:,tmp,:]<rho_min, rho[:,tmp,:])

         total_transp = (vhnew * rhopot2[:,tmp,:] * sig2 / sig2).sum() # in kg/s
         # set masked values to zero
         if total_transp is np.ma.masked: total_transp = 0.0

         print 'AABW mass flux (kg/s)', total_transp

         return total_transp # in kg/s

def buoyancy_flux_ishelf(Qall,lprec,area_ice_shelf,hml,temp,salt,h,SST,SSS):
   '''
   Surface buoyancy flux due to sub-ice-shelf melting - Eq. 9 in Chapman et al, 1998 JPO
   '''
   RHO_ICE = 918.0
   RHO0 = 1028.0
   L = 3.34e5
   g = 9.8
   Q = (Qall)/(RHO_ICE*L)# sea ice growth rate
   # change in rho due to lprec
   volume = hml * area_ice_shelf
   # fresh water mass
   mass_fw = (lprec * area_ice_shelf * 3600.0*24) # in kg and over one day
   rho = eos.wright_eos(SST,SSS,2.0e7)
   mass_ocean = volume*rho # in kg
   salt2 = mass_ocean * SSS/(mass_ocean - mass_fw) # salinity due to lprec flux
   print 'ice shelf melting, mass_fw,mass_ocean',mass_fw.sum(), mass_ocean.sum()
   salt2 = salt2 * area_ice_shelf/area_ice_shelf
   print 'salt2 min/max',salt2.min(), salt2.max()
   delta_rho = rho - eos.wright_eos(SST,salt2,2.0e7)
   S = (delta_rho * g)/RHO0
   B = (S*Q*area_ice_shelf).sum() # in m4/s3
   print 'ice shelf melting, delta_rho min/max',delta_rho.min(), delta_rho.max()

   return B


def buoyancy_flux_cshelf(Qall,saltf,lprec,area_cshelf,hml,temp,salt,h,SST,SSS):
   '''
   Surface buoyancy flux due to sea ice formation - Eq. 9 in Chapman et al, 1998 JPO
   '''
   RHO_ICE = 905.0
   RHO0 = 1028.0
   L = 3.34e5
   g = 9.8
   Q = (Qall)/(RHO_ICE*L)# sea ice growth rate
   # change in rho due to salt and lprec
   saltf = saltf * area_cshelf * 1.0e3 * 3600.0*24 # in g and over one day
   volume = hml * area_cshelf
   # fresh water mass
   mass_fw = (lprec * area_cshelf * 3600.0*24) # in kg and over one day
   rho = eos.wright_eos(SST,SSS,2.0e7)
   mass_ocean = volume*rho # in kg
   total_salt = mass_ocean * SSS # in g
   # old_mass should be > mass_ocean if sea ice is produced
   old_mass = mass_ocean - mass_fw # in kg
   old_total_salt = total_salt + saltf # in g
   old_SSS = old_total_salt/old_mass
   old_rho = eos.wright_eos(SST,old_SSS,2.0e7)
   print 'old_rho min/max',old_rho.min(),old_rho.max()
   print 'old_SSS min/max',old_SSS.min(),old_SSS.max()
   print 'old_mass',old_mass.sum()
   print 'mass_fw,mass_ocean',mass_fw.sum(), mass_ocean.sum()
   #salt2 = saltf/(mass_ocean + mass_fw) # salinity due to salt flux
   #delta_rho = rho - eos.wright_eos(SST,SSS-salt2,2.0e7)
   delta_rho = rho - old_rho
   S = (delta_rho * g)/RHO0
   B = (S*Q*area_cshelf).sum() # in m4/s3
   print 'delta_rho min/max',delta_rho.min(), delta_rho.max()

   return B

#def get_ML_state(temp,salt,hml,h):
#    '''
#    Compute mean tracer properties in the mixed layer.
#    '''
#    km,jm,im = h.shape
#    depth = np.zeros((jm,im)); SST = np.zeros((jm,im)); SSS = np.zeros((jm,im))
#    for j in range(jm):
#        for i in range(im):
#            for k in range(km):
#                if (depth[j,i] + h[k,j,i] < hml[j,i]):
#                    dh = h[k,j,i]
#                elif (depth[j,i] < hml[j,i]):
#                    dh = hml[j,i] - depth[j,i]
#                else:
#                    dh = 0.0
#
#                SST[j,i] = SST[j,i] + dh * temp[k,j,i]
#                SSS[j,i] = SSS[j,i] + dh * salt[k,j,i]
#
#            depth[j,i] = depth[j,i] + dh
#
#    return SST/depth, SSS/depth

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
   AABW.units = 'kg/s'; AABW.description = 'Offshore mass flux of AABW computed at the shelf break'
   hml_shelf = ncfile.createVariable('hml_shelf',np.dtype('float32').char,('time'))
   hml_shelf.units = 'm'; hml_shelf.description = 'ML depth ave. over continental shelf region'

   Bseaice = ncfile.createVariable('Bseaice',np.dtype('float32').char,('time'))
   Bseaice.units = 'm^4/s^3'; Bseaice.description = 'Surface buoyancy flux due to sea ice formation (integrated over the shelf area).'

   Biceshelf = ncfile.createVariable('Biceshelf',np.dtype('float32').char,('time'))
   Biceshelf.units = 'm^4/s^3'; Biceshelf.description = 'Surface buoyancy flux due to sub-ice-shelf melting (integrated over the ice shelf area).'

   Nshelf = ncfile.createVariable('Nshelf',np.dtype('float32').char,('time'))
   Nshelf.units = '1/s'; Nshelf.description = 'Mean ambient strat. freq. over the cont. shelf'

   Tshelf = ncfile.createVariable('Tshelf',np.dtype('float32').char,('time'))
   Tshelf.units = 'C'; Tshelf.description = 'Mean pot. temp. over the cont. shelf'

   Sshelf = ncfile.createVariable('Sshelf',np.dtype('float32').char,('time'))
   Sshelf.units = 'ppt'; Sshelf.description = 'Mean salinity over the cont. shelf'

   Ustar_shelf = ncfile.createVariable('Ustar_shelf',np.dtype('float32').char,('time'))
   Ustar_shelf.units = 'm/s'; Ustar_shelf.description = 'Mean ustar over the cont. shelf'

   sflux = ncfile.createVariable('sflux',np.dtype('float32').char,('time'))
   sflux.units = 'kg/(m^2 s)'; sflux.description = 'Total salt flux over the cont. shelf'

   qflux = ncfile.createVariable('qflux',np.dtype('float32').char,('time'))
   qflux.units = 'W/m2'; qflux.description = 'Total heat flux over the cont. shelf'

   seaiceVolume = ncfile.createVariable('seaiceVolume',np.dtype('float32').char,('time'))
   seaiceVolume.units = 'm^3'; seaiceVolume.description = 'total volume of sea ice'

   # write data to coordinate vars.
   x[:] = xx[:]
   y[:] = yy[:]
   time[:] = ocean_time[:]  # in years

   # close the file.
   ncfile.close()
   print ('*** SUCCESS creating '+exp_name+'.nc!')

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
