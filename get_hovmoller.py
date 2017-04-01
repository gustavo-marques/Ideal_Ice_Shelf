#!/usr/bin/env python

# Gustavo Marques

from remapping import mom_remapping
from scipy import signal
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import os
import argparse
from numpy.fft import fft
import numpy
from scipy.stats import pearsonr

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """
  parser = argparse.ArgumentParser(description=
      '''
      Compute... .
      ''',
  epilog='Written by Gustavo Marques, Mar. 2017.')

  parser.add_argument('-name', type=str, default='test', help='''Name of the experiment (default is test).''')

  parser.add_argument('-j', type=int, default='450', help='''j indice (default is 450).''')

  parser.add_argument('-n', type=int, default='1', help='''Compute every n values (default is 1).''')

  parser.add_argument('-ti', type=int, default='0', help='''Initial time indice (default is 0).''')
  parser.add_argument('-tf', type=int, default='-1', help='''Final time indice (default is len(time)).''')

  optCmdLineArgs = parser.parse_args()
  driver(optCmdLineArgs)

def driver(args):
  j = args.j # section location in y
  n = args.n # plot every n
  ti = args.ti
  if args.tf == -1:
     tf = len(netCDF4.Dataset('prog.nc').variables['time'][:])
  else:
     tf = args.tf

  # create dir locally
  os.system('mkdir -p PNG')

  # static variables
  depth = netCDF4.Dataset('ocean_geometry.nc').variables['D'][j,:]
  x = netCDF4.Dataset('ocean_geometry.nc').variables['geolon'][j,:]
  ave_ssh = netCDF4.Dataset('IDEAL_IS_IC.nc').variables['ave_ssh'][0,j,:]
  z = np.linspace(0,depth.max(),100)
  h = 0.5* (z[0:-1]+z[1::])

  # time
  tind = np.arange(ti,tf,n)
  tm = len(tind)
  # remapping params
  cs = mom_remapping.Remapping_Cs()
  # PCM
  #cs.remapping_scheme = 0 # PCM
  # PPM-H4
  cs.remapping_scheme = 2
  cs.degree = 2

  # points/sections
  time = np.zeros((tm))
  temp_p1 = np.zeros((tm))
  v_bar_p1 = np.zeros((tm))
  v_p1 = np.zeros((tm))
  v_p2 = np.zeros((tm))
  v_p3 = np.zeros((tm))
  v_s1 = np.zeros((tm,len(x))); v_s1_f = np.zeros((tm,len(x)))
  v_s2 = np.zeros((tm,len(x))); v_s2_f = np.zeros((tm,len(x)))
  v_s3 = np.zeros((tm,len(x))); v_s3_f = np.zeros((tm,len(x)))
  v_bar = np.zeros((tm,len(x))); v_bar_f = np.zeros((tm,len(x)))
  # main loop
  for t in range(tm):
      time[t] = netCDF4.Dataset('prog.nc').variables['time'][tind[t]] # in days
      print 'Time is:',time[t]
      data0 = netCDF4.Dataset('prog.nc').variables['u'][tind[t],:,j,:]
      h0 = netCDF4.Dataset('prog.nc').variables['h'][tind[t],:,j,:]
      e0 = netCDF4.Dataset('prog.nc').variables['e'][tind[t],:,j,:]
      temp0 = netCDF4.Dataset('prog.nc').variables['temp'][tind[t],:,j,:]

      # remap
      data1 = remap(cs,h0,e0,h,data0,z,depth)
      temp1 = remap(cs,h0,e0,h,temp0,z,depth)
      # vbar
      v_bar[t,:] = (data0 * h0).sum(axis = 0)/h0.sum(axis=0)

      # sections
      v_s1[t,:] = data1[0,:] # surface
      v_s2[t,:] = data1[44,:] # mid-depth
      v_s3[t,:] = data1[-3,:] # bottom
      # points
      v_bar_p1[t] = v_bar[t,140]
      v_p1[t] = data1[0,140]
      v_p2[t] = data1[44,140]
      v_p3[t] = data1[-3,140]
      # temp
      temp_p1[t] = temp1[25,125]

  print 'Done with loop, now plotting...'

  # auto-corr temp.
  #temp_p1 = temp_p1 - temp_p1.mean()
  #temp_p1_corr = lagcorr(temp_p1,temp_p1,lag=25,verbose=True)
  #print 'temp_p1_corr.shape', temp_p1_corr.shape, temp_p1_corr[:,0]
  #plt.figure()
  #plt.subplot(211)
  #plt.plot(time,temp_p1,'k')
  #plt.title('Temperature anomaly')
  #plt.grid()
  #plt.subplot(212)
  #plt.plot(temp_p1_corr[0,:])
  #plt.title('Lagged-correlation')
  #plt.grid()
  #filename = str('PNG/%s_temp_time_series' % (args.name)) + '.png'
  #plt.savefig(filename, dpi=100)

  # filtering
  fps = 2
  lowcut = 0.075; highcut = 0.25

  filtered_p1 = butter_bandpass_filter(v_p1, lowcut, highcut, fps, order=6)
  filtered_p2 = butter_bandpass_filter(v_p2, lowcut, highcut, fps, order=6)
  filtered_p3 = butter_bandpass_filter(v_p3, lowcut, highcut, fps, order=6)
  filtered_v_bar_p1 = butter_bandpass_filter(v_bar_p1, lowcut, highcut, fps, order=6)
  #filtered_p1 = butter_highpass_filter(v_p1,0.1,fps)
  #filtered_p2 = butter_highpass_filter(v_p2,0.1,fps)
  #filtered_p3 = butter_highpass_filter(v_p3,0.1,fps)

  for i in range(len(x)):
     v_s1_f[:,i] = butter_bandpass_filter(v_s1[:,i],lowcut, highcut, fps, order=6)
     v_s2_f[:,i] = butter_bandpass_filter(v_s2[:,i],lowcut, highcut, fps, order=6)
     v_s3_f[:,i] = butter_bandpass_filter(v_s3[:,i],lowcut, highcut, fps, order=6)
     v_bar_f[:,i] = butter_bandpass_filter(v_bar[:,i],lowcut, highcut, fps, order=6)
     #v_s1_f[:,i] = butter_highpass_filter(v_s1[:,i],0.1,fps)
     #v_s2_f[:,i] = butter_highpass_filter(v_s2[:,i],0.1,fps)
     #v_s3_f[:,i] = butter_highpass_filter(v_s3[:,i],0.1,fps)

  # compute spectrum
  computeSpectrum(args.name,time,v_p1,'v_p1')
  computeSpectrum(args.name,time,v_p2,'v_p2')
  computeSpectrum(args.name,time,v_p3,'v_p3')
  computeSpectrum(args.name,time,v_bar_p1,'v_bar_p1')

  # plots
  plt.figure(figsize=(20,10))
  plt.subplot(211)
  plt.plot(time,v_p1,'k',time,v_p2,'b',time,v_p3,'r',time,v_bar_p1,'g')
  plt.title('raw signal')
  plt.grid()
  plt.subplot(212)
  plt.plot(time,filtered_p1,'k',time,filtered_p2,'b',time,filtered_p3,'r',time,filtered_v_bar_p1,'g')
  plt.title('filtered signal')
  plt.grid()
  #ax.set_ylabel(r'Across-slope vel. [m s$^{-1}$] ',fontsize=15)
  #ax.set_xlabel('time [days]',fontsize=15)
  filename = str('PNG/%s_time_series' % (args.name)) + '.png'
  plt.savefig(filename, dpi=100)

  X,Z = np.meshgrid(x,time)

  tmp = np.abs(v_s1_f.mean())
  tmp = 0.05
  levs = np.linspace(-tmp,tmp,100)
  fig, ax = plt.subplots()
  ct = ax.contourf(X,Z,v_s1_f,levs, cmap=plt.cm.bwr)
  #cbar = fig.colorbar(ct, ticks=[-0.1, -0.05, 0, 0.05, 0.1], orientation = 'horizontal')
  cbar = fig.colorbar(ct, orientation = 'horizontal')
  cbar.set_label(r'Across-slope vel. [m s$^{-1}$]',fontsize=15)
  ax.set_xlim(x.min(),x.max())
  #ax.set_ylim(-720,0)
  #ss = str("ax.set_title('Time: %6.4f years')"% (time))
  ax.set_xlabel('x [km]',fontsize=15)
  ax.set_ylabel('time [days]',fontsize=15)
  filename = str('PNG/%s_v_s1' % (args.name)) + '.png'
  plt.savefig(filename, dpi=100)

  #tmp = np.abs(v_s2_f.mean())
  #levs = np.linspace(-tmp,tmp,25)
  fig, ax = plt.subplots()
  ct = ax.contourf(X,Z,v_s2_f,levs, cmap=plt.cm.bwr)
  cbar = fig.colorbar(ct, orientation = 'horizontal')
  cbar.set_label(r'Across-slope vel. [m s$^{-1}$]',fontsize=15)
  ax.set_xlim(x.min(),x.max())
  #ax.set_ylim(-720,0)
  #ss = str("ax.set_title('Time: %6.4f years')"% (time))
  ax.set_xlabel('x [km]',fontsize=15)
  ax.set_ylabel('time [days]',fontsize=15)
  filename = str('PNG/%s_v_s2' % (args.name)) + '.png'
  plt.savefig(filename, dpi=100)

  #tmp = np.abs(v_s3_f.mean())
  #levs = np.linspace(-tmp,tmp,25)
  fig, ax = plt.subplots()
  ct = ax.contourf(X,Z,v_s3_f,levs, cmap=plt.cm.bwr)
  cbar = fig.colorbar(ct, orientation = 'horizontal')
  cbar.set_label(r'Across-slope vel. [m s$^{-1}$]',fontsize=15)
  ax.set_xlim(x.min(),x.max())
  #ax.set_ylim(-720,0)
  #ss = str("ax.set_title('Time: %6.4f years')"% (time))
  ax.set_xlabel('x [km]',fontsize=15)
  filename = str('PNG/%s_v_s3' % (args.name)) + '.png'
  plt.savefig(filename, dpi=100)

  fig, ax = plt.subplots()
  ct = ax.contourf(X,Z,v_bar_f,levs, cmap=plt.cm.bwr)
  cbar = fig.colorbar(ct, orientation = 'horizontal')
  cbar.set_label(r'Across-slope depth ave. vel. [m s$^{-1}$]',fontsize=15)
  ax.set_xlim(x.min(),x.max())
  #ax.set_ylim(-720,0)
  #ss = str("ax.set_title('Time: %6.4f years')"% (time))
  ax.set_xlabel('x [km]',fontsize=15)
  ax.set_ylabel('time [days]',fontsize=15)
  filename = str('PNG/%s_vbar' % (args.name)) + '.png'
  plt.savefig(filename, dpi=100)


  return

def remap(cs,h0,e0,h,data0,z,depth):
   km,im = h0.shape
   data1 = np.ma.zeros((len(h),im))
   for i in range(im):
     h1 = np.diff(z)
     if h0[:,i].sum() > 0.01: # ocean
	h1[h<-e0[0,i]] = 0.0 # top
	h1[h>depth[i]] = 0.0 # bottom
	# need to account for SSH and make sure total thicknesses are
	# the same
	dh = h0[:,i].sum() - h1.sum()
        tmp1 = np.nonzero(h1!=0.0)[0]
        if len(tmp1)>0:
	  if dh > 0.0:
	     # correct thickness in the top non-vanished layer
             h1[tmp1[0]] = h1[tmp1[0]] + dh # add
          elif dh < 0.0:
             h1[tmp1[0]] = h1[tmp1[0]] - dh # remove
        else:
           data1[:,i] = np.ma.masked

        # for debugging
        #if h0[:,i].sum() != h1.sum():
	#   print 'WARNING: dh, h0[:,i].sum(), h1.sum()',dh, h0[:,i].sum(), h1.sum()

        # remap
	data1[:,i] = mom_remapping.remapping_core_h(h0[:,i], data0[:,i], h1, cs)
	# mask
	data1[h1==0.0,i] = np.ma.masked;
     else: # land/iceshelf
	data1[:,i] = np.ma.masked

   return data1

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def computeSpectrum(id,time,y,var):
    """
    Compute a Single-Sided Amplitude Spectrum of y(t)
    id : id of experiment, as a string
    time : time vector, in days
    y : signal
    var : name of variable to be saved

    """
    time=np.asarray(time)
    y=np.asarray(y)
    tmp=np.nonzero(time>=50)[0]
    time=time[tmp]*3600.*24. # in sec
    y=y[tmp]

    n = len(y)
    k = np.arange(n)

    Fs = 1 / ((time[1]-time[0])) * 86400.# in cpd
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    tmp = fft(y) # fft computing and normalization
    Y = tmp*tmp.conjugate()
    Y = Y[range(n/2)]

    # save data

    #file1 = str("'TXT/%s-spectra-%s.txt'" % (id,var))
    #s1 = csv.writer(open(eval(file1), 'wb'), delimiter=',')
    #for l in range(1,len(frq)):
    #   #s1.writerow([str(frq[l])]  + [str(abs(Y[l]/n))])
    #   s1.writerow([str(frq[l])]  + [str(abs(Y[l]))])

    #print 'Done saving spectra!'

    plt.figure()
    plt.loglog(frq,abs(Y))
    plt.xlabel('Freq (cpd)')
    plt.ylabel('PSD')
    plt.grid()
    filename = str('PNG/%s-spectra-%s' % (id,var)) + '.png'
    plt.savefig(filename, dpi=100)

    return

def lagcorr(x,y,lag=None,verbose=True):
    '''Compute lead-lag correlations between 2 time series.

    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    '''

    if len(x)!=len(y):
        raise('Input variables of different lengths.')

    #--------Unify types of <lag>-------------
    if numpy.isscalar(lag):
        if abs(lag)>=len(x):
            raise('Maximum lag equal or larger than array.')
        if lag<0:
            lag=-numpy.arange(abs(lag)+1)
        elif lag==0:
            lag=[0,]
        else:
            lag=numpy.arange(lag+1)
    elif lag is None:
        lag=[0,]
    else:
        lag=numpy.asarray(lag)

    #-------Loop over lags---------------------
    result=[]
    if verbose:
        print '\n#<lagcorr>: Computing lagged-correlations at lags:',lag

    for ii in lag:
        if ii<0:
            result.append(pearsonr(x[:ii],y[-ii:]))
        elif ii==0:
            result.append(pearsonr(x,y))
        elif ii>0:
            result.append(pearsonr(x[ii:],y[:-ii]))

    result=numpy.asarray(result)

    return result

# Invoke parseCommandLine(), the top-level prodedure
if __name__ == '__main__': parseCommandLine()

