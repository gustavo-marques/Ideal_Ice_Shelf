#!/usr/bin/env python

# generate 3D diagnostics 
# Gustavo Marques, Sep. 2016

import argparse
import netCDF4 
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pyvtk
import os,sys
import pp

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

  parser.add_argument('--bergs', help='''Generates 3D VTK data using icebergs.''', action="store_true")

  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver3D(optCmdLineArgs)

def driver3D(args):
    if args.time != 0 and args.ncpus > 1:
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

    # load time
    if args.time > 0 :
        time = np.array(netCDF4.Dataset(args.oceanfile).variables['time'][args.time])
        tind = [args.time]
    else:
        time = netCDF4.Dataset(args.oceanfile).variables['time'][:]
        tind = range(0,len(time),args.dt)

    if args.bergs:
       rho_berg = 918.0
       rho = 1030.
       # mass of bergs can time dependent
       mass_berg = netCDF4.Dataset(args.oceanfile).variables['mass_berg'][tind,:,:]
       IS=(mass_berg/rho_berg)

    else:
       # ice shelf thickness, static for now
       IS = netCDF4.Dataset('MOM_Shelf_IC.nc').variables['h_shelf'][:]
    
    # interface and layer thickness
    e=netCDF4.Dataset(args.oceanfile).variables['e'][0,:,:,:]
    h=netCDF4.Dataset(args.oceanfile).variables['h'][0,:,:,:]
    # correct top and bottom, for vis pourposes
    h[0,:,:]=e[0,:,:]; h[-1,:,:]=e[-1,:,:]
    NZ,NY,NX=h.shape
    # create VTK bathymetry
    VTKgen(lats,lons,D.mask,depth=D,h=h,fname=name)
    
    #if not tind>1:
        # create VTK ice-shelf
    print 'Saving ice shelf...'
    VTKgen(lats,lons,D.mask,h=h,shelf_base=ssh,shelf_thick=IS,fname=name)

    # list to tuple
    time_ind = tuple(tind)
    # Execution starts as soon as one of the workers will become available
    jobs = [(t, job_server.submit(get_data, (t, lats, lons, D, dx, dy, NX, NY, args, ), 
           (VTKgen,nan_helper,f1,f3,),
           ("os","pyvtk","netCDF4","numpy",))) for t in time_ind]

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
       for i in range(len(time)):
            f.write('%8.4f \n' % (time[i]))
       f.close()


    print ' \n' + '==> ' + '  DONE!\n' + ''

def get_data(t, lats, lons, D, dx, dy, NX, NY, args):
    """
    Reads data from netcdf at time indice t and writes a VTK file.
    """
    # check if data has been saved
    path_to_file = str('VTK/%s-%05d.vtk' % (args.n,t))
    if os.path.isfile(path_to_file):
       print ' \n' + '==> ' + 'FILE EXISTS, MOVING TO THE NEXT ONE ...\n' + ''
    else:
       print 'Saving time indice:', t
       # ocean
       # structure
       # layer thickness
       file = netCDF4.Dataset(args.oceanfile)
       e=netCDF4.Dataset(args.oceanfile).variables['e'][t,:,:,:]
       h=netCDF4.Dataset(args.oceanfile).variables['h'][t,:,:,:]
       h_dum=numpy.abs(numpy.diff(e,axis=0))
       hz=0.5*(e[0:-1,:,:]+e[1::,:,:])# cell depth
       temp=netCDF4.Dataset(args.oceanfile).variables['temp'][t,:,:,:]
       salt=netCDF4.Dataset(args.oceanfile).variables['salt'][t,:,:,:]
       sig2=netCDF4.Dataset(args.oceanfile).variables['rhopot2'][t,:,:,:] - 1000.
       if 'tr2' in file.variables:
          tr2=netCDF4.Dataset(args.oceanfile).variables['tr2'][t,:,:,:]
       else:
          tr2=numpy.zeros(salt.shape)

       if 'tr3' in file.variables:
          tr3=netCDF4.Dataset(args.oceanfile).variables['tr3'][t,:,:,:]
       else:
          tr3=numpy.zeros(salt.shape)

       # for isopycnal models
       # temp, mark values where h_dum<10 cm with Nan
       temp[h_dum<0.01]=numpy.nan; salt[h_dum<0.01]=numpy.nan
       sig2[h_dum<0.01] = numpy.nan; tr2[h_dum<0.01] = numpy.nan; tr3[h_dum<0.01] = numpy.nan
       # same for salt, u and v
       v=netCDF4.Dataset(args.oceanfile).variables['vh'][t,:,:,:]/(h*dx)
       u=netCDF4.Dataset(args.oceanfile).variables['uh'][t,:,:,:]/(h*dy)
       u=numpy.ma.masked_where(numpy.abs(u)>5,u)
       v=numpy.ma.masked_where(numpy.abs(v)>5,v)
       v.fill_value=numpy.nan
       u.fill_value=numpy.nan
       u[h_dum<0.01]=numpy.nan
       v[h_dum<0.01]=numpy.nan
       v=v.filled(); u=u.filled()
       # interpolate nan values
       for j in range(NY):
	     for i in range(NX):
		 nans, tmp = nan_helper(temp[:,j,i])
		 if nans.mask.any() == False:
			if h_dum[:,j,i].sum() > 0.0:
			   temp[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], temp[~nans,j,i])
			   salt[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], salt[~nans,j,i])
			   sig2[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], sig2[~nans,j,i])
			   tr2[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], tr2[~nans,j,i])
			   tr3[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], tr3[~nans,j,i])
			   u[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], u[~nans,j,i])
			   v[nans,j,i]= numpy.interp(-hz[nans,j,i], -hz[~nans,j,i], v[~nans,j,i])

       # write just bottom values
       #VTKgen(lats,lons,D.mask,depth=D,h=h,temp=tt,salt=ss,rho=gamma,u=uu,v=vv,writebottom=True,fname=reg,t=ind)

       print 'Saving ocean data... \n'
       VTKgen(lats,lons,D.mask,h=hz,temp=temp,salt=salt,rho=sig2,dye2=tr2,dye3=tr3,u=u,v=v,fname=args.n,t=t)
       if args.bergs:
           print 'Saving bergs data... \n'
           # save ice shelf made of icebergs
           VTKgen(lats,lons,D.mask,h=hz,shelf_base=e[0,:,:],shelf_thick=IS[t,:,:],fname=args.n,t=t)

    return

def VTKgen(lat,lon,mask,depth=None,h=None,temp=None,salt=None,rho=None,dye1=None,dye2=None,dye3=None,u=None,v=None,w=None,seaice=None,shelf_base=None,shelf_thick=None,writebottom=False,fname='test',dirname='VTK',date=None, t=0):
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

       if writebottom == True:
          print ' \n' + '==> ' + 'Writing tracers/vel. just at the bottom layer ...\n' + ''
          data=[]
          if temp is not None:
             tmp=np.zeros((2,NY,NX)) 
             if len(temp.shape)==2: # in case the user provides 2D array with bottom data
                 tmp[:,:,:]=temp[:,:]
             else:
                 tmp[:,:,:]=temp[-1,:,:]
                 
             temp=f1(tmp)
             data.append("pyvtk.Scalars(temp,name='Temp')")

          if salt is not None:
             tmp=np.zeros((2,NY,NX))
             if len(salt.shape)==2:
                 tmp[:,:,:]=salt[:,:] 
             else:
                tmp[:,:,:]=salt[-1,:,:]

             salt=f1(tmp)
             data.append("pyvtk.Scalars(salt,name='Salt')")

          if rho is not None:
             tmp=np.zeros((2,NY,NX))
             if len(rho.shape)==2:
                 tmp[:,:,:]=rho[:,:]
             else:
                tmp[:,:,:]=rho[-1,:,:]

             rho=f1(tmp)
             data.append("pyvtk.Scalars(rho,name='Rho')")

          if dye1 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye1.shape)==2:
                 tmp[:,:,:]=dye1[:,:]
             else:
                tmp[:,:,:]=dye1[-1,:,:]

             dye1=f1(tmp)
             data.append("pyvtk.Scalars(dye1,name='Dye1')") 

          if dye2 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye2.shape)==2:
                 tmp[:,:,:]=dye2[:,:]
             else:
                tmp[:,:,:]=dye2[-1,:,:]

             dye2=f1(tmp)
             data.append("pyvtk.Scalars(dye2,name='Dye2')")

          if dye3 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye3.shape)==2:
                 tmp[:,:,:]=dye3[:,:]
             else:
                tmp[:,:,:]=dye3[-1,:,:]

             dye3=f1(tmp)
             data.append("pyvtk.Scalars(dye3,name='Dye3')")

          if u is not None and v is not None:
                w=np.zeros((2,NY,NX)) # no w vel for now
                tmpu=np.zeros((2,NY,NX))
                tmpv=np.zeros((2,NY,NX))
                if len(u.shape)==2:
                    tmpu[:,:,:]=u[:,:]
                    tmpv[:,:,:]=v[:,:]
                else:
                    tmpu[:,:,:]=u[-1,:,:]
                    tmpv[:,:,:]=v[-1,:,:]

                vel=f3(tmpu,tmpv,w)
                data.append("pyvtk.Vectors(vel,name='Velocity')")

          if temp is not None or salt is not None or rho is not None or u is not None:

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
                 s = str("vtk.tofile('%s/%s-%s-bottom-%05d','binary')" % (dirname,fname,date,t))
                 eval(s)

              else:
                 s = str("vtk.tofile('%s/%s-bottom-%05d','binary')" % (dirname,fname,t))
                 eval(s)
   
            
 
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
       if temp is not None:
         temp=f1(temp)
         data.append("pyvtk.Scalars(temp,name='Temp')")
         
       if salt is not None:
         salt=f1(salt)
         data.append("pyvtk.Scalars(salt,name='Salt')")

       if rho is not None:
         rho=f1(rho)
         data.append("pyvtk.Scalars(rho,name='Rho')")

       if dye1 is not None:
         dye1=f1(dye1)
         data.append("pyvtk.Scalars(dye1,name='Dye1')")

       if dye2 is not None:
         dye2=f1(dye2)
         data.append("pyvtk.Scalars(dye2,name='Dye2')")

       if dye3 is not None:
         dye3=f1(dye3)
         data.append("pyvtk.Scalars(dye3,name='Dye3')")

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
         sice2=np.zeros((NZ,NY,NX))
         seaice[seaice>=0.15]=1.0 # all values >=15% are unit
         sice2[0,:,:]=seaice[:,:]
         seaice1=f1(sice1) 
         seaice2=f1(sice2) 
         data.append("pyvtk.Scalars(seaice1,name='Sea-ice')")
         data.append("pyvtk.Scalars(seaice2,name='Sea-ice-binary')")
       

       if temp is not None or salt is not None or rho is not None or u is not None or seaice is not None:
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
             s = str("vtk.tofile('%s/%s-%s-%05d','binary')" % (dirname,fname,date,t))
             eval(s)

         else:
             s = str("vtk.tofile('%s/%s-%05d','binary')" % (dirname,fname,t))
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
