#!/usr/bin/env python

# generate 3D diagnostics 
# Gustavo Marques, Sep. 2016

import argparse
from netCDF4 import MFDataset, Dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
from pyvtk import *
import os,sys
#import pymp

#pymp.config.thread_limit = 3
#pymp.config.nested = True

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

  parser.add_argument('--oceanfile', type=str, default='prog.nc',
      help='''Name of the netCDF file with the ocean variables. Default is prog.nc.''')

  parser.add_argument('--bergs', help='''Generates 3D VTK data using icebergs.''', action="store_true")

  optCmdLineArgs = parser.parse_args()
  global name
  name = optCmdLineArgs.n
  driver3D(optCmdLineArgs)

def driver3D(args):
    ocean_file = args.oceanfile 
    # ocean grid
    D=Dataset('ocean_geometry.nc').variables['D'][:]
    lonh=Dataset('ocean_geometry.nc').variables['lonh'][:]
    lonq=Dataset('ocean_geometry.nc').variables['lonq'][:]
    lath=Dataset('ocean_geometry.nc').variables['lath'][:]
    latq=Dataset('ocean_geometry.nc').variables['latq'][:]
    dx = (lonh[1] - lonh[0]) * 1.0e3 # in m
    dy = (lath[1] - lath[0]) * 1.0e3 # in m
    lonqs, latqs = np.meshgrid(lonq,latq)
    lons, lats = np.meshgrid(lonh,lath)
    D=np.ma.masked_where(D <= 1, D)
    D.mask = np.ma.array(D); D.mask[:,:]=False

    # ice shelf base
    ssh = Dataset('IDEAL_IS_IC.nc').variables['ave_ssh'][0,:,:]

    # load time
    print('Time indice is:' + str(args.time))
    if args.time > 0 :
        time = np.array(Dataset(ocean_file).variables['time'][args.time])
        tind = [args.time]
    else:
        time = Dataset(ocean_file).variables['time'][:]
        tind = range(len(time))

    if args.bergs:
       rho_berg = 918.0
       rho = 1030.
       # mass of bergs can time dependent
       mass_berg = Dataset(ocean_file).variables['mass_berg'][tind,:,:]
       IS=(mass_berg/rho_berg)

    else:
       # ice shelf thickness, static for now
       IS = Dataset('MOM_Shelf_IC.nc').variables['h_shelf'][:]
    
    # interface and layer thickness
    e=Dataset(ocean_file).variables['e'][0,:,:,:]
    h=Dataset(ocean_file).variables['h'][0,:,:,:]
    # correct top and bottom, for vis pourposes
    h[0,:,:]=e[0,:,:]; h[-1,:,:]=e[-1,:,:]
    NZ,NY,NX=h.shape
    # create VTK bathymetry
    VTKgen(lats,lons,D.mask,depth=D,h=h,fname=name)
    
    #if not tind>1:
        # create VTK ice-shelf
    print 'Saving ice shelf...'
    VTKgen(lats,lons,D.mask,h=h,shelf_base=ssh,shelf_thick=IS,fname=name)

    time_list=[] # list for time 

    print 'len(tind)',len(tind)
    # loop through time and plot surface fields
    for t in range(0,len(tind),args.dt):
        print 'Time is:',t
        #print 'Time is:',time[t]
        # save data in the lists
        #time_list.append(date)
        # check if data has been saved
        path_to_file = str('VTK/%s-%05d.vtk' % (name,tind[t]))
        if os.path.isfile(path_to_file):
           print ' \n' + '==> ' + 'FILE EXISTS, MOVING TO THE NEXT ONE ...\n' + ''
        else:
           print 'Saving time indice:', tind[t] 
           # ocean
           # structure
           # layer thickness
           e=Dataset(ocean_file).variables['e'][tind[t],:,:,:]
           h=Dataset(ocean_file).variables['h'][tind[t],:,:,:]
           h_dum=np.abs(np.diff(e,axis=0))
           hz=0.5*(e[0:-1,:,:]+e[1::,:,:])# cell depth
           temp=Dataset(ocean_file).variables['temp'][tind[t],:,:,:]
           salt=Dataset(ocean_file).variables['salt'][tind[t],:,:,:]
	   sig2=Dataset(ocean_file).variables['rhopot2'][tind[t],:,:,:] - 1000.
	   tr2=Dataset(ocean_file).variables['tr2'][tind[t],:,:,:]
	   tr3=Dataset(ocean_file).variables['tr3'][tind[t],:,:,:]
           # for isopycnal models
           # temp, mark values where h_dum<10 cm with Nan
           temp[h_dum<0.01]=np.nan; salt[h_dum<0.01]=np.nan
	   sig2[h_dum<0.01] = np.nan; tr2[h_dum<0.01] = np.nan; tr3[h_dum<0.01] = np.nan
           # same for salt, u and v
           v=Dataset(ocean_file).variables['vh'][tind[t],:,:,:]/(h*dx)
           u=Dataset(ocean_file).variables['uh'][tind[t],:,:,:]/(h*dy)
	   u=np.ma.masked_where(np.abs(u)>5,u)
	   v=np.ma.masked_where(np.abs(v)>5,v)
	   v.fill_value=np.nan
	   u.fill_value=np.nan
           u[h_dum<0.01]=np.nan
	   v[h_dum<0.01]=np.nan
	   v=v.filled(); u=u.filled()
	   # interpolate nan values
	   for j in range(NY):
		 for i in range(NX):
			 nans, tmp = nan_helper(temp[:,j,i])
			 if nans.mask.any() == False:
		            if h_dum[:,j,i].sum() > 0.0:
			       temp[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], temp[~nans,j,i])
			       salt[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], salt[~nans,j,i])
			       sig2[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], sig2[~nans,j,i])
			       tr2[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], tr2[~nans,j,i])
			       tr3[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], tr3[~nans,j,i])
			       u[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], u[~nans,j,i])
			       v[nans,j,i]= np.interp(-hz[nans,j,i], -hz[~nans,j,i], v[~nans,j,i])

           # write just bottom values
           #VTKgen(lats,lons,D.mask,depth=D,h=h,temp=tt,salt=ss,rho=gamma,u=uu,v=vv,writebottom=True,fname=reg,t=ind)

           print 'Saving ocean data... \n'
           VTKgen(lats,lons,D.mask,h=hz,temp=temp,salt=salt,rho=sig2,dye2=tr2,dye3=tr3,u=u,v=v,fname=name,t=tind[t])
           if args.bergs:
              print 'Saving bergs data... \n'
              # save ice shelf made of icebergs
              VTKgen(lats,lons,D.mask,h=hz,shelf_base=e[0,:,:],shelf_thick=IS[t,:,:],fname=name,t=tind[t])

    print ' \n' + '==> ' + ' Done saving vtk files!\n' + ''
    # the following is for displaying the time in the PNG figures generated by Visit
    if args.time==0:
       print ' \n' + '==> ' + ' Writting time to asn ascii file... \n' + ''
       f = open('VTK/time.txt', 'wb')
       for i in range(len(time)):
            f.write('%8.4f \n'.format(time[i])) 
       f.close()


    print ' \n' + '==> ' + '  DONE!\n' + ''

def VTKgen(lat,lon,mask,depth=None,h=None,temp=None,salt=None,rho=None,dye1=None,dye2=None,dye3=None,u=None,v=None,w=None,seaice=None,shelf_base=None,shelf_thick=None,writebottom=False,fname='test',dirname='VTK',date=None, t=0):
    """ Writes ocean and ice shelf geometry and data (e.g., tracer, vel., sea-ice) into VTK files."""

    NY,NX=lat.shape
    if not os.path.isdir(dirname):
        os.system('mkdir ' + dirname)

    NY,NX=lat.shape

    if depth is not None:
       newlat=np.resize(lat,(2,NY,NX))
       newlon=np.resize(lon,(2,NY,NX))
       newdepth,bottom=get_depth(h,depth,mask)
       pp = f3(newlon,newlat,newdepth)
       structure=StructuredGrid([2,NY,NX],pp)
       path_to_file = str('%s/%s-bathymetry.vtk' % (dirname,fname))
       if os.path.isfile(path_to_file):
           print ' \n' + '==> ' + 'Bathymetry has already been written, moving on ...\n' + ''
       else:
          # create bottom/shape and depths
          newdepth=f1(newdepth)
          bottom=f1(bottom)
          pointdata = PointData(Scalars(newdepth,name='Depth'), Scalars(bottom,name='Bottom9999'))
          # saving the data
          vtk = VtkData(structure,pointdata)
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
             data.append("Scalars(temp,name='Temp')")

          if salt is not None:
             tmp=np.zeros((2,NY,NX))
             if len(salt.shape)==2:
                 tmp[:,:,:]=salt[:,:] 
             else:
                tmp[:,:,:]=salt[-1,:,:]

             salt=f1(tmp)
             data.append("Scalars(salt,name='Salt')")

          if rho is not None:
             tmp=np.zeros((2,NY,NX))
             if len(rho.shape)==2:
                 tmp[:,:,:]=rho[:,:]
             else:
                tmp[:,:,:]=rho[-1,:,:]

             rho=f1(tmp)
             data.append("Scalars(rho,name='Rho')")

          if dye1 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye1.shape)==2:
                 tmp[:,:,:]=dye1[:,:]
             else:
                tmp[:,:,:]=dye1[-1,:,:]

             dye1=f1(tmp)
             data.append("Scalars(dye1,name='Dye1')") 

          if dye2 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye2.shape)==2:
                 tmp[:,:,:]=dye2[:,:]
             else:
                tmp[:,:,:]=dye2[-1,:,:]

             dye2=f1(tmp)
             data.append("Scalars(dye2,name='Dye2')")

          if dye3 is not None:
             tmp=np.zeros((2,NY,NX))
             if len(dye3.shape)==2:
                 tmp[:,:,:]=dye3[:,:]
             else:
                tmp[:,:,:]=dye3[-1,:,:]

             dye3=f1(tmp)
             data.append("Scalars(dye3,name='Dye3')")

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
                data.append("Vectors(vel,name='Velocity')")

          if temp is not None or salt is not None or rho is not None or u is not None:

              for d in range(len(data)):
                 if d==0:
                    tmp=data[d]
                 else:
                    tmp=tmp+','+data[d]

              s = str("PointData(%s)" % (tmp))
              pointdata = eval(s)
              # saving the data
              vtk = VtkData(structure,pointdata)
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
       structure=StructuredGrid([NZ_IS,NY,NX],pp)
       pointdata = PointData(Scalars(iceshelf,name='IceShelf9999'))
       vtk = VtkData(structure,pointdata)
       if t>0:
           s = str("vtk.tofile('%s/%s-ice-shelf-%05d','binary')" % (dirname,fname,t))
           eval(s)

       else:
           vtk.tofile(dirname+'/'+fname+'-ice-shelf','binary')

    if writebottom==False:
       data=[]
       if temp is not None:
         temp=f1(temp)
         data.append("Scalars(temp,name='Temp')")
         
       if salt is not None:
         salt=f1(salt)
         data.append("Scalars(salt,name='Salt')")

       if rho is not None:
         rho=f1(rho)
         data.append("Scalars(rho,name='Rho')")

       if dye1 is not None:
         dye1=f1(dye1)
         data.append("Scalars(dye1,name='Dye1')")

       if dye2 is not None:
         dye2=f1(dye2)
         data.append("Scalars(dye2,name='Dye2')")

       if dye3 is not None:
         dye3=f1(dye3)
         data.append("Scalars(dye3,name='Dye3')")

       if u is not None and v is not None:
         if w is not None:
            vel=f3(u,v,w)
         else:
            w=np.zeros(u.shape)
            vel=f3(u,v,w)

         data.append("Vectors(vel,name='Velocity')")
    
       if seaice is not None:
         NZ,NY,NX=h.shape
         sice1=np.zeros((NZ,NY,NX))
         sice1[0,:,:]=seaice[:,:]
         sice2=np.zeros((NZ,NY,NX))
         seaice[seaice>=0.15]=1.0 # all values >=15% are unit
         sice2[0,:,:]=seaice[:,:]
         seaice1=f1(sice1) 
         seaice2=f1(sice2) 
         data.append("Scalars(seaice1,name='Sea-ice')")
         data.append("Scalars(seaice2,name='Sea-ice-binary')")
       

       if temp is not None or salt is not None or rho is not None or u is not None or seaice is not None:
         NZ,NY,NX=h.shape
         # resize lon lat for real mesh
         newlat=np.resize(lat,(NZ,NY,NX))
         newlon=np.resize(lon,(NZ,NY,NX))
         pp = f3(newlon,newlat,h)
         structure=StructuredGrid([NZ,NY,NX],pp)

         for d in range(len(data)):
            if d==0:
               tmp=data[d]
            else:
               tmp=tmp+','+data[d]

         s = str("PointData(%s)" % (tmp))
         pointdata = eval(s)
         # saving the data
         vtk = VtkData(structure,pointdata)
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

    return np.isnan(y), lambda z: z.nonzero()[0]
# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()
