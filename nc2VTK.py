#!/usr/bin/env python

# generate 3D diagnostics
# Gustavo Marques, Sep. 2016

import argparse
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from remapping import mom_remapping
import remapping # for pp to work
import warnings
import pyvtk
import os,sys
import pp
import time

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Generated 3D VTK diagnostics for the ISOMIP+ test cases. This script
      assumes that files/variables names are consistent with the default ISOMIP/GFDL
      options (see: https://github.com/NOAA-GFDL/MOM6-examples/tree/dev/master/ocean_only/ISOMIP).
      ''',
  epilog='Written by Gustavo Marques, Sep. 2016.')

  parser.add_argument('-n', type=str, default='Ocean0_COM_MOM6',
      help='''The name of the experiment following the ISOMIP+ definition (expt_COM_model). This name is used to save the VTK files. Default is Ocean0_COM_MOM6.''')

  parser.add_argument('-time', type=int, default=0,
      help='''The time indice to save the VTK files. Default value is 0, which saves the entire dataset. If --time > 0, only one time instance will be saved.''')

  parser.add_argument('-dt', type=int, default=1,
      help='''The time indice interval to save the VTK files. Default value is 1, which saves the entire dataset.''')

  parser.add_argument('-ncpus', type=int, default=1,
      help='''The number of cpus to use when time = 0 (i.e., the entire dataset will be saved). Default value is 1.''')

  parser.add_argument('--oceanfile', type=str, default='prog.nc',
      help='''Name of the netCDF file with the ocean variables. Default is prog.nc.''')

  parser.add_argument('--icefile', type=str, default='ice_month.nc',
      help='''Name of the netCDF file with the sea ice variables. Default is ice_month.nc.''')

  parser.add_argument('--tracer', type=str, default='salt',
      help='''Name of the tracer in netCDF file. Default is salt.''')

  parser.add_argument('--bergs', help='''Generates 3D VTK data using icebergs.''', action="store_true")

  parser.add_argument('--seaice', help='''Write sea ice thickness.''', action="store_true")

  parser.add_argument('--save_vel', help='''If true, save u and v.''', action="store_true")

  parser.add_argument('-time_stats', help='''Prints the approx. time spent at different blocks of code.''', action="store_true")

  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver3D(optCmdLineArgs)

def driver3D(args):
    t0 = time.time()
    if args.time > 0 and args.ncpus > 2:
       print 'Error! Parameter -ncpus cannot be > 1 when time =! 0.'
       quit()

    # tuple of all parallel python servers to connect with
    ppservers = ()

    if args.ncpus > 1: # run in parallel
       # Creates jobserver with ncpus workers
       job_server = pp.Server(args.ncpus, ppservers=ppservers)
    else:
       # Creates jobserver with automatically detected number of workers
       job_server = pp.Server(ppservers=ppservers)

    print "Starting pp with", job_server.get_ncpus(), "workers"

    args.oceanfile = args.oceanfile
    # ocean grid
    D=netCDF4.Dataset('ocean_geometry.nc').variables['D'][:]
    lonh=netCDF4.Dataset('ocean_geometry.nc').variables['lonh'][:]
    lonq=netCDF4.Dataset('ocean_geometry.nc').variables['lonq'][:]
    lath=netCDF4.Dataset('ocean_geometry.nc').variables['lath'][:]
    latq=netCDF4.Dataset('ocean_geometry.nc').variables['latq'][:]
    dx = (lonh[1] - lonh[0]) * 1.0e3 # in m
    dy = (lath[1] - lath[0]) * 1.0e3 # in m
    lonqs, latqs = np.meshgrid(lonq,latq)
    lons, lats = np.meshgrid(lonh,lath)
    D=np.ma.masked_where(D <= 1, D)
    D.mask = np.ma.array(D); D.mask[:,:]=False

    # ice shelf base
    ssh = netCDF4.Dataset('IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,:]

    t1 = time.time()
    if args.time_stats:
       print 'Time in initialization:',(t1-t0)

    t0 = time.time()
    # load time
    if args.time > 0 :
        times = np.array(netCDF4.MFDataset(args.oceanfile).variables['time'][args.time])
        tind = [args.time]
    else:
        times = netCDF4.MFDataset(args.oceanfile).variables['time'][:]
        tind = range(0,len(times),args.dt)

    t1 = time.time()
    if args.time_stats:
       print 'Time loading time:',(t1-t0)

    t0 = time.time()
    if args.bergs:
       rho_berg = 918.0
       rho = 1030.
       # mass of bergs can time dependent
       mass_berg = netCDF4.MFDataset(args.oceanfile).variables['mass_berg'][tind,:,:]
       IS=(mass_berg/rho_berg)

    else:
       # ice shelf thickness, static for now
       IS = netCDF4.MFDataset('MOM_Shelf_IC.nc').variables['h_shelf'][:]

    t1 = time.time()
    if args.time_stats:
       print 'Time loading ice shelf data:',(t1-t0)

    t0 = time.time()
    # interface and layer thickness
    e=netCDF4.MFDataset(args.oceanfile).variables['e'][0,:,:,:]
    h=netCDF4.MFDataset(args.oceanfile).variables['h'][0,:,:,:]
    # correct top and bottom, for vis pourposes
    h[0,:,:]=e[0,:,:]; h[-1,:,:]=e[-1,:,:]
    NZ,NY,NX=h.shape


    # create VTK bathymetry
    VTKgen(lats,lons,D.mask,depth=D,h=h,fname=name)

    t1 = time.time()
    if args.time_stats:
       print 'Time writing topography:',(t1-t0)

    t0 = time.time()

    #if not tind>1:
        # create VTK ice-shelf
    print 'Saving ice shelf...'
    VTKgen(lats,lons,D.mask,h=h,shelf_base=ssh,shelf_thick=IS,fname=name)

    t1 = time.time()
    if args.time_stats:
       print 'Time writing ice shelf data:',(t1-t0)

    t0 = time.time()

    # list to tuple
    time_ind = tuple(tind)
    # Execution starts as soon as one of the workers will become available
    jobs = [(t, job_server.submit(get_data, (t, lats, lons, D, dx, dy, NX, NY, args,),
           (VTKgen,nan_helper,f1,f3,),
           ("os","pyvtk","netCDF4","numpy","time",))) for t in time_ind]

#    jobs = [(t, job_server.submit(get_data1, (t,args, ), (VTKgen,),
#        ("netCDF4",))) for t in tind]

    for input, job in jobs:
       job()

    job_server.print_stats()

    print ' \n' + '==> ' + ' Done saving vtk files!\n' + ''
    # the following is for displaying the time in the PNG figures generated by Visit
    if args.time==0:
       print ' \n' + '==> ' + ' Writting time to asn ascii file... \n' + ''
       f = open('VTK/time.txt', 'wb')
       for i in range(len(times)):
            f.write('%8.4f \n' % (times[i]))
       f.close()


    print ' \n' + '==> ' + '  DONE!\n' + ''

def get_data(t, lats, lons, D, dx, dy, NX, NY, args):
    """
    Reads data from netcdf at time indice t and writes a VTK file.
    """
    # check if data has been saved
    if args.save_vel:
       path_to_file = str('VTK/%s-uv-%05d.vtk' % (args.n,t))
    else:
       path_to_file = str('VTK/%s-%s-%05d.vtk' % (args.n,args.tracer,t))

    if os.path.isfile(path_to_file):
       print ' \n' + '==> ' + 'FILE EXISTS, MOVING TO THE NEXT ONE ...\n' + ''
    else:
       print 'Saving time indice:', t
       # ocean
       # structure
       t0 = time.time()
       #file = netCDF4.MFDataset(args.oceanfile)
       e=netCDF4.MFDataset(args.oceanfile).variables['e'][t,:,:,:]
       h=netCDF4.MFDataset(args.oceanfile).variables['h'][t,:,:,:]
       h_dum=numpy.abs(numpy.diff(e,axis=0))
       hz=0.5*(e[0:-1,:,:]+e[1::,:,:])# cell depth
       if args.save_vel:
         data=0
       else:
	 data=netCDF4.MFDataset(args.oceanfile).variables[args.tracer][t,:,:,:]

       t1 = time.time()
       if args.time_stats:
	  print 'Main loop, time loading data:',(t1-t0)

       t0 = time.time()

       #data1 =remap(cs,h0,e0,h,data0,z,D)
       # mark values where h_dum<10 cm with Nan
       data[h_dum<0.01]=numpy.nan
       # interpolate nan values
       for j in range(NY):
	     for i in range(NX):
		 nans, tmp = nan_helper(data[:,j,i])
		 if nans.mask.any() == False:
			if h_dum[:,j,i].sum() > 0.0:
			   data[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], data[~nans,j,i])

       t1 = time.time()
       if args.time_stats:
	  print 'Main loop, time processing data:',(t1-t0)

       t0 = time.time()

       # write files
       print 'Saving ocean data... \n'
       VTKgen(lats,lons,D.mask,h=hz,tracer=data,fname=args.n,varname=args.tracer,t=t)
       t1 = time.time()

       if args.seaice:
	  print 'Saving sea ice data... \n'
          HI=netCDF4.MFDataset(args.icefile).variables['HI'][t,:,:]
          HI.fill_value=numpy.nan
          VTKgen(lats,lons,D.mask,h=hz,seaice=HI,fname=args.n,varname='seaice',t=t)

   #    if args.bergs:
   #       print 'Saving bergs data... \n'
   #       # save ice shelf made of icebergs
#	  VTKgen(lats,lons,D.mask,h=hz,shelf_base=e[0,:,:],shelf_thick=IS[t,:,:],fname=args.n,t=t)

       if args.time_stats:
          print 'Main loop, time writing VTK:',(t1-t0)

    return

def VTKgen(lat,lon,mask,depth=None,h=None,tracer=None,u=None,v=None,w=None,seaice=None,shelf_base=None,shelf_thick=None,writebottom=False,fname='test',dirname='VTK',varname='salt',date=None, t=0):
    """ Writes ocean and ice shelf geometry and data (e.g., tracer, vel., sea-ice) into VTK files."""

    import numpy as np

    NY,NX=lat.shape
    if not os.path.isdir(dirname):
        os.system('mkdir ' + dirname)

    NY,NX=lat.shape
    if depth is not None:
       newlat=np.resize(lat,(2,NY,NX))
       newlon=np.resize(lon,(2,NY,NX))
       newdepth,bottom=get_depth(h,depth,mask)
       pp = f3(newlon,newlat,newdepth)
       structure=pyvtk.StructuredGrid([2,NY,NX],pp)
       path_to_file = str('%s/%s-bathymetry.vtk' % (dirname,fname))
       if os.path.isfile(path_to_file):
           print ' \n' + '==> ' + 'Bathymetry has already been written, moving on ...\n' + ''
       else:
          # create bottom/shape and depths
          newdepth=f1(newdepth)
          bottom=f1(bottom)
          pointdata = pyvtk.PointData(pyvtk.Scalars(newdepth,name='Depth'), pyvtk.Scalars(bottom,name='Bottom9999'))
          # saving the data
          vtk = pyvtk.VtkData(structure,pointdata)
          vtk.tofile(dirname+'/'+fname+'-bathymetry','binary')


    if shelf_base is not None and shelf_thick is not None:
       NZ_IS = 2
       newlat=np.resize(lat,(NZ_IS,NY,NX))
       newlon=np.resize(lon,(NZ_IS,NY,NX))
       dum,z = get_iceshelf(shelf_base,shelf_thick,NZ_IS)
       iceshelf = f1(dum)
       pp = f3(newlon,newlat,z)
       structure=pyvtk.StructuredGrid([NZ_IS,NY,NX],pp)
       pointdata = pyvtk.PointData(pyvtk.Scalars(iceshelf,name='IceShelf9999'))
       vtk = pyvtk.VtkData(structure,pointdata)
       if t>0:
           s = str("vtk.tofile('%s/%s-ice-shelf-%05d','binary')" % (dirname,fname,t))
           eval(s)

       else:
           vtk.tofile(dirname+'/'+fname+'-ice-shelf','binary')

    if writebottom==False:
       data=[]
       if tracer is not None:
         tracer=f1(tracer)
         #s = str('data.append("pyvtk.Scalars(tracer,name=''%s'')")' % (varname))
         #eval(s)
         if varname == 'salt':
            data.append("pyvtk.Scalars(tracer,name='salt')")
         elif varname == 'temp':
            data.append("pyvtk.Scalars(tracer,name='temp')")
         elif varname == 'rho':
            data.append("pyvtk.Scalars(tracer,name='rho')")
         elif varname == 'tr1':
            data.append("pyvtk.Scalars(tracer,name='tr1')")
         elif varname == 'tr2':
            data.append("pyvtk.Scalars(tracer,name='tr2')")
         else:
            data.append("pyvtk.Scalars(tracer,name='tracer')")

       if u is not None and v is not None:
         if w is not None:
            vel=f3(u,v,w)
         else:
            w=np.zeros(u.shape)
            vel=f3(u,v,w)

         data.append("pyvtk.Vectors(vel,name='Velocity')")

       if seaice is not None:
         NZ,NY,NX=h.shape
         sice1=np.zeros((NZ,NY,NX))
         sice1[0,:,:]=seaice[:,:]
         #sice2=np.zeros((NZ,NY,NX))
         #seaice[seaice>=0.15]=1.0 # all values >=15% are unit
         #sice2[0,:,:]=seaice[:,:]
         seaice1=f1(sice1)
         #seaice2=f1(sice2)
         data.append("pyvtk.Scalars(seaice1,name='SeaIceThickness')")
         #data.append("pyvtk.Scalars(seaice2,name='Sea-ice-binary')")


       if tracer is not None or u is not None or seaice is not None:
         NZ,NY,NX=h.shape
         # resize lon lat for real mesh
         newlat=np.resize(lat,(NZ,NY,NX))
         newlon=np.resize(lon,(NZ,NY,NX))
         pp = f3(newlon,newlat,h)
         structure=pyvtk.StructuredGrid([NZ,NY,NX],pp)

         for d in range(len(data)):
            if d==0:
               tmp=data[d]
            else:
               tmp=tmp+','+data[d]

         s = str("pyvtk.PointData(%s)" % (tmp))
         pointdata = eval(s)
         # saving the data
         vtk = pyvtk.VtkData(structure,pointdata)
         if date is not None:
             s = str("vtk.tofile('%s/%s-%s-%s-%05d','binary')" % (dirname,fname,varname,date,t))
             eval(s)

         else:
             s = str("vtk.tofile('%s/%s-%s-%05d','binary')" % (dirname,fname,varname,t))
             eval(s)


def get_iceshelf(bottom,thick,NZ):
    ''' Return an array with the vertical position of the ice shelf.'''
    NY,NX=bottom.shape
    tmp=np.resize(thick,(NZ,NY,NX))
    ice_shelf=np.zeros((NZ,NY,NX))
    z_shelf=np.zeros((NZ,NY,NX))
    z_shelf[0,:,:]=bottom[:,:]+0.1
    for i in range(NX):
	    z_shelf[1,:,i]=(bottom[:,0] + thick[0,:,0]) + 0.1

    z_shelf = np.ma.masked_where(tmp==0,z_shelf)
    z_shelf.fill_value=np.nan
    ice_shelf[:]=9999
    ice_shelf = np.ma.masked_where(tmp==0,ice_shelf)
    ice_shelf.fill_value=np.nan
    #return ice_shelf,z_shelf
    return z_shelf,z_shelf

def get_depth(h,D,mask):
    NZ,NY,NX=h.shape
    dep=np.zeros((2,NY,NX))
    bot=np.zeros((2,NY,NX))
    bot[0,:,:]=9999
    # start at min(D)
    dep[0,:,:]=-D.max()
    for i in range(NX):
      for j in range(NY):
         if mask[j,i]==False:
                 dep[1,j,i]=-D[j,i] # ocean/ice-shelf
         else:
                 dep[1,j,i]=0 # land

    return dep,bot

def f1(X):
    """function to sort one 3D array"""

    a,b,c=X.shape
    return [(X[k,i,j]) for j in range(c) for i in range(b) for k in range(a)]

def f3(X,Y,Z):
    """ function to sort three 3D matrices"""
    a,b,c=X.shape
    return [(X[k,i,j],Y[k,i,j],Z[k,i,j]) for j in range(c) for i in range(b) for k in range(a)]


def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return numpy.isnan(y), lambda z: z.nonzero()[0]
# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
