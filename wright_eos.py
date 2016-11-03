import numpy
import gsw
import matplotlib.pyplot as plt

"""
==============================

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================
"""


def wright_eos(T,S,p):
  """
  
 **********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
    Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns density [kg m-3]
 
 ***********************************************************************
 
 """

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  rho = (p + p0) * I_denom

  return rho 

def alpha_wright_eos(T,S,p):
  """

**********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
   Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns partial derivative of density with respect to temperature [kg m-3 C-1]

 ***********************************************************************
 
 """
  
  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  I_denom2 = I_denom*I_denom

  drho_dT =  I_denom2*(lam*(b1+T*(2*b2 + 3*b3*T) + b5*S) - (p+p0)*((p+p0)*a1 + (c1+T*(2*c2 + 3*c3*T) + c5*S)))


  return drho_dT

def beta_wright_eos(T,S,p):
  """

 **********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
   Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns partial derivative of density with respect to salinity [kg m-3 PSU-1]
 
 ***********************************************************************
 
 """
  
  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  I_denom2 = I_denom*I_denom
  drho_dS =  I_denom2*(lam*(b4+b5*T) - (p+p0)*((p+p0)*a2 + (c4+c5*T)))

  return drho_dS


def tfreeze(S,P):
  TFr_S0_P0 = 0.0788; dTFr_dS = -0.0572; dTFr_dp = -7.77E-08
  T_Fr = (TFr_S0_P0 + dTFr_dS*S) + dTFr_dp*P
  return T_Fr


def linearize_tfreeze():
  '''
   Compute The coefficients for the freezing point.
   Copied from Nico's script
  '''
  (S, p) = numpy.meshgrid(numpy.linspace(20.,40.,21), numpy.linspace(0, 1e3, 21))
  lat = -75.*numpy.pi/180.
  lon = -105.*numpy.pi/180.
  SA = gsw.SA_from_SP(S, p, lon, lat)
  CT_freezing = gsw.CT_freezing(SA, p)
  pt_freezing = gsw.pt_from_CT(SA, CT_freezing)
  t_freezing = gsw.t_from_CT(SA, CT_freezing, p)
  t_lin = -0.0573*S + 0.0832 -7.53e-8*(p*1e4)
  nSamples = len(S.ravel())
  bases = numpy.zeros((nSamples,3))
  bases[:,0] = S.ravel()
  bases[:,1] = 1.
  bases[:,2] = p.ravel()
  A = numpy.matmul(bases.T, bases)
  b = numpy.matmul(bases.T, pt_freezing.ravel())
  coeffs = numpy.linalg.solve(A,b)
  print coeffs[0], coeffs[1], coeffs[2]*1e-4
  pt_lin = coeffs[0]*S + coeffs[1] + coeffs[2]*p
  plt.figure(1)
  plt.pcolor(S, p, pt_freezing)
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.figure(2)
  plt.pcolor(S, p, pt_freezing - t_freezing)
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.figure(3)
  plt.pcolor(S, p, pt_freezing-pt_lin)
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.xlabel('S (PSU)')
  plt.ylabel('p (dbar)')
  plt.figure(4)
  plt.pcolor(S, p, t_freezing-t_lin)
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.figure(5)
  plt.pcolor(S, p, pt_freezing-t_lin)
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.xlabel('S (PSU)')
  plt.ylabel('p (dbar)')
  plt.show()
  return
