# ==========================================
# Signal Fitting Functions and 
#   Initial Parameters
#
# Migle Stankaityte
# 6 May 2019

#Modified for Z->bb calibration
#QT PhD project 2019-2020
# ==========================================

#import
import sys, re
import numpy as np
import ROOT
from ROOT import TF1, gROOT
gROOT.SetBatch(True)

# ==========================================

# Fitted functions
#-------------------------------------------
# Getting the right function for template

def getFitFn(template,fname,mmin,mmax,wbin,conf,cand) :
  
  # --------- Conf Note ----------
  if conf :  
  # Higgs : 3 gaussians
    if 'Higgs' in template : 
      print 'Fitting Higgs template with 3 gaussians (conf. note setup).'
      f = fitHiggs_confnote(fname,re.sub(r'f\_Higgs','',fname),mmin,mmax,wbin)

    # Z boson : 3 gaussians
    elif 'Zboson' in template :  
      print 'Fitting Z boson template with 3 gaussians (conf. note setup).'
      f = fitZboson_confnote(fname,re.sub(r'f\_Zboson','',fname),mmin,mmax,wbin)
    
    # W boson : 3 gaussians
    elif 'Wboson' in template:  
      print 'Fitting W boson template with 3 gaussians (conf. note setup).'
      f = fitWboson_confnote(fname,re.sub(r'f\_Wboson','',fname),mmin,mmax,wbin)
    
    # ttbar : DSCB
    elif 'ttbar' in template : 
      print 'Fitting ttbar template with DSCB (conf. note setup)'
      f = fit_ttbar_confnote(fname,re.sub(r'f\_ttbar','',fname),mmin,mmax,wbin)
  
  # --------- Full Dataset ----------
  elif 'Higgs' in template : 
    print 'Fitting Higgs template with DSCB.'
    f = fitHiggs(fname,re.sub(r'f\_Higgs','',fname),mmin,mmax,wbin,cand)

  # Z boson : 3 gaussians
  elif 'Zboson' in template :  
    print 'Fitting Z boson template with DCSB'
    f = fitZboson(fname,re.sub(r'f\_Zboson','',fname),mmin,mmax,wbin,cand)
    #f = fitZboson_bukin(fname,re.sub(r'f\_Zboson','',fname),mmin,mmax,wbin)
  
  # W boson : 3 gaussians
  elif 'Wboson' in template:  
    print 'Fitting W boson template with DCSB'
    f = fitWboson(fname,re.sub(r'f\_Wboson','',fname),mmin,mmax,wbin,cand)
    #f = fitWboson(fname,re.sub(r'f\_Wboson','',fname),mmin,mmax,wbin,cand)

  # ttbar : DSCB
  #elif 'ttbar' in template :  #inclusive ttbar
    #print 'Fitting ttbar template with DSCB.'
    #f = fit_ttbar(fname,re.sub(r'f\_ttbar','',fname),mmin,mmax,wbin,cand)

  elif 'ttbar1' in template : #ttbar top in large-R jet
    print 'Fitting ttbar Category 1 template with DSCB.'
    f = fit_ttbar(fname,re.sub(r'f\_ttbar1','',fname),mmin,mmax,wbin,cand) 

  elif 'ttbar2' in template : #ttbar W in large-R jet+Others
    print 'Fitting ttbar Category 2 template with DSCB.'
    f = fit_ttbar_cat2(fname,re.sub(r'f\_ttbar2','',fname),mmin,mmax,wbin,cand) 

  elif 'ttbar2W' in template : #ttbar W in large-R jet+Others
    print 'Fitting ttbar Category 2 W template with DSCB.'
    f = fit_ttbar_cat2W(fname,re.sub(r'f\_ttbar2','',fname),mmin,mmax,wbin,cand) 

  elif 'ttbar2O' in template : #ttbar W in large-R jet+Others
    print 'Fitting ttbar Category 2 Other template with DSCB.'
    f = fit_ttbar_cat2O(fname,re.sub(r'f\_ttbar2','',fname),mmin,mmax,wbin,cand)   

  # dijets : expo polynomial
  elif 'dijets' in template :
    print 'Fitting dijets template with expo polynomial.'
    f = fit_dijets(fname,re.sub(r'f\_dijets','',fname),mmin,mmax,wbin,cand)

    # --------- Unknown ----------
  else :
    raise Exception('Unknown template: {}'.format(template))
 
  return f

# End of: getFitFn()
# ==========================================

# Double sided CB definition

def doublecrystal(inx,par) :
  x   = inx[0]
  N   = par[0] 
  mu  = par[1] 
  sig = par[2] 
  a1  = par[3]
  p1  = par[4]
  a2  = par[5]
  p2  = par[6]

  t = (x-mu)/sig
 
  if t < -a1 :
    a = np.exp(-0.5*a1*a1)
    b = p1/a1 - a1
    return N*a/np.power(a1/p1*(b - t), p1)
  
  elif t > a2 : 
    a = np.exp(-0.5*a2*a2)
    b = p2/a2 - a2
    return N*(a/np.power(a2/p2*(b + t), p2))
  
  return N*np.exp(-0.5*t*t)


# End of: doublecrystal()
# ==========================================
# ******************************************
def bukin(x, par):

    Xp = par[0]
    sigp = par[1]
    xi = par[2]
    rho1 = par[3]
    rho2 = par[4]
    ap = par[5]
    consts = 2*np.sqrt(2*np.log(2.0))
    r1=0
    r2=0
    r3=0
    r4=0
    r5=0
    hp=0
    x1 = 0
    x2 = 0
    fit_result = 0
    hp=sigp*consts
    r3=np.log(2.)
    r4=np.sqrt(xi*xi+1)
    r1=xi/r4
    if abs(xi) > np.exp(-6.):
        r5=xi/np.log(r4+xi)
    else:
        r5=1
    x1 = Xp + (hp / 2) * (r1-1)
    x2 = Xp + (hp / 2) * (r1+1)
    #--- Left Side
    if x[0] < x1:
        r2=rho1*(x[0] - x1)*(x[0]-x1)/(Xp-x1)/(Xp-x1)-r3+4*r3*(x[0]-x1)/hp*r5*r4/(r4-xi)/(r4-xi)
    #--- Center
    elif x[0] < x2:
        if abs(xi) > np.exp(-6.):
            r2=np.log(1 + 4 * xi * r4 * (x[0] - Xp)/hp)/np.log(1 + 2*xi*(xi - r4))
            r2=-r3*r2*r2
        else:
            r2=-4*r3*(x[0] - Xp)*(x[0] - Xp)/hp/hp
    #--- Right Side
    else:
         r2=rho2*(x[0] - x2)*(x[0] - x2)/(Xp - x2)/(Xp - x2)-r3 - 4 * r3 * (x[0] - x2)/hp * r5 * r4/(r4 + xi)/(r4 + xi)
    if abs(r2) > 100:
        fit_result = 0
    else:
    #---- Normalize the result
        fit_result = np.exp(r2)
    return fit_result *ap


# STARTING PARAMTERS / FIT FUNCTIONS
# (here be dragons)

# Higgs: DSCB

def fitHiggs(fname,systname,mmin,mmax,wbin,cand) :

  print 'Systematic: ', systname 
  
  f = TF1(fname,doublecrystal,mmin,mmax,7)

  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,6.75089e+01)    # norm
  f.SetParameter(1,1.24151e+02)    # mean
  f.SetParameter(2,9.26263e+00)    # sigma
  f.SetParameter(3,6.98991e-01)    # a1
  f.SetParameter(4,9.99886e+04)    # p1
  f.SetParameter(5,2.06096e+00)    # a2
  f.SetParameter(6,8.00505e-01)    # p2
  f.FixParameter(7,0)
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,4.26329e+01)    # norm
    f.SetParameter(1,1.14136e+02)    # mean
    f.SetParameter(2,1.26103e+01)    # sigma
    f.SetParameter(3,6.70565e-01)    # a1
    f.SetParameter(4,9.99993e+04)    # p1
    f.SetParameter(5,2.05157e+00)    # a2
    f.SetParameter(6,6.82528e-01)    # p2
    f.FixParameter(7,0)

  # ------------------------------------------
  # constraints
  
  M = 125.18 
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
 
  return f

# End of: fitHiggs_confnote()
# ==========================================
#Trying 3 Gaussians initially for fitting Zboson

def fitZboson_3gaus(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 
  f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",mmin,mmax);
  #f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) ) + [6]",mmin,mmax);

  # best fit values with 85-115% constraints on two out of 3 mu
  f.SetParameter(0,3.10995e+03)      # norm
  f.SetParameter(1,8.20689e+01)      # mean1 
  f.SetParameter(2,9.23617e+00)      # sigma1
  f.SetParameter(3,5.44394e+00)      # N2/norm
  f.SetParameter(4,9.43995e+01)      # mean2
  f.SetParameter(5,8.33186e+00)      # sigma2
  f.SetParameter(6,1.63850e+00)      # N3/norm
  f.SetParameter(7,1.00305e+02)      # mean3
  f.SetParameter(8,4.93243e+01)      # sigma3
  f.SetParameter(9,2.04081e+01)      # const.
  

  # ------------------------------------------
  # constraints
  
  M = 91.1876 
  
  f.SetParLimits(0, 0.    , 1e+05  )      # norm
  f.SetParLimits(1, 0.85*M , 1.3*M )      # mean1
  f.SetParLimits(2, wbin  , 50.   )      # sigma1
  f.SetParLimits(3, 0.    , 1e+05  )      # N2/norm
  f.SetParLimits(4, 0.9*M , 1.1*M )      # mean2
  f.SetParLimits(5, wbin  , 50.   )      # sigma2
  f.SetParLimits(6, 0.    , 1e+05  )      # N3/norm
  f.SetParLimits(7, 0.9*M , 1.1*M )      # mean3
  f.SetParLimits(8, wbin  , 50.   )      # sigma3
  f.SetParLimits(9, 0.    , 1e+05  )      # const.

  return f


# W boson: 3 gaussian
def fitWboson_3gaus(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 
  #f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",mmin,mmax);
  f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1)) + [6]",mmin,mmax);
  
  # best fit values with 85-115% constraints on two out of 3 mu
  f.SetParameter(0,1.88833e+03)      # norm
  f.SetParameter(1,7.23412e+01)      # mean1
  f.SetParameter(2,2.40493e+01)      # sigma1
  f.SetParameter(3,7.05691e-01)      # N2/norm
  f.SetParameter(4,8.02729e+01)      # mean2
  f.SetParameter(5,5.00000e+00)      # sigma2
  f.SetParameter(6,4.88712e-01)      # N3/norm
  f.SetParameter(7,8.84169e+01)      # mean3
  f.SetParameter(8,2.06415e+01)      # sigma3
  f.SetParameter(9,2.02217e-01)      # const.

  if 'herwig' in systname :
    # initialise sub-peaks in the center since Herwig
    # does not have the high mass stat-fluctuations
    f.SetParameter(1,8.03790e+01)      # mean1
    f.SetParameter(7,8.03790e+01)      # mean3



  # ------------------------------------------
  # constraints
  
  M = 80.379
  
  f.SetParLimits(0, 0.    , 1e+05  )      # norm
  f.SetParLimits(1, 0.95*M , 1.0*M )      # mean1
  f.SetParLimits(2, wbin  , 50.   )      # sigma1
  f.SetParLimits(3, 0.    , 1e+05  )      # N2/norm
  f.SetParLimits(4, 0.9*M , 1.1*M )      # mean2
  f.SetParLimits(5, wbin  , 50.   )      # sigma2
  f.SetParLimits(6, 0.    , 1e+05  )      # N3/norm
  f.SetParLimits(7, 0.9*M , 1.1*M )      # mean3
  f.SetParLimits(8, wbin  , 50.   )      # sigma3
  f.SetParLimits(9, 0.    , 1e+05  )      # const.
  
  return f
# Z boson: DSCB

def fitZboson_bukin(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 

  f = TF1(fname,bukin,mmin,mmax,6)

  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,9.11876e+01)    # Xp
  f.SetParameter(1,5.00000e+01)    # sigp
  f.SetParameter(2,3.11022e-03)    # xi      
  f.SetParameter(3,4.10995e-01)    # rho1
  f.SetParameter(4,3.87565e-01)    # rho2
  f.SetParameter(5,5.91673e+01)    # ap
  #f.SetParameter(6,7.08804e-01)    # p2
  
  # subleading
  #if cand == 's' :
    #f.SetParameter(0,1.88038e+03)    # norm
    #f.SetParameter(1,8.74564e+01)    # mean
    #f.SetParameter(2,9.79574e+00)    # sigma
    #f.SetParameter(3,5.49219e-01)    # a1
    #f.SetParameter(4,9.99998e+04)    # p1
    #f.SetParameter(5,1.38511e+00)    # a2
    #f.SetParameter(6,9.25181e-01)    # p2
  
  # ------------------------------------------
  # constraints
  
  #M = 90
  #M = 85
  M = 91.1876
  #M=94.0

  f.SetParLimits(0, .85*M, 1.15*M )    # Xp
  f.SetParLimits(1, wbin  , 50 )    # sigp
  f.SetParLimits(2, 0.    , 1e+03 )       # xi
  f.SetParLimits(3, 0.    , 1e+03 )    # rho1
  f.SetParLimits(4, 0.    , 1e+03 )    # rho2
  f.SetParLimits(5, 0.    , 1e+03 )    # ap
  #f.SetParLimits(6, 0.    , 1e+05 )    # p2

  return f


# End of: fitZboson()
# ==========================================
def fitZboson(fname,systname,mmin,mmax,wbin,cand) :

  print 'Systematic: ', systname 

  f = TF1(fname,doublecrystal,mmin,mmax,7)

  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,2.91673e+03)    # norm
  f.SetParameter(1,9.28007e+01)    # mean
  f.SetParameter(2,1.39571e+01)    # sigma 8.39571e+00
  f.SetParameter(3,9.37115e-01)    # a1
  f.SetParameter(4,9.99999e+04)    # p1
  f.SetParameter(5,1.68133e+00)    # a2
  f.SetParameter(6,7.08804e-01)    # p2
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,1.88038e+03)    # norm
    f.SetParameter(1,8.74564e+01)    # mean
    f.SetParameter(2,9.79574e+00)    # sigma
    f.SetParameter(3,5.49219e-01)    # a1
    f.SetParameter(4,9.99998e+04)    # p1
    f.SetParameter(5,1.38511e+00)    # a2
    f.SetParameter(6,9.25181e-01)    # p2
  
  # ------------------------------------------
  # constraints
  
  #M = 90
  #M = 85
  M = 91.1876
  #M=94.0
  
  #Some of these parameters needs change based on shape of the templates, better to look for cross-check plost 
  #once and decides these parameters
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.90*M , 1.1*M )    # mean
  #f.SetParLimits(1, 0.95*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 60 )       # sigma #for 2-b tag, Xbb
  #f.SetParLimits(2, wbin  , 50 )       # sigma  #trials for 1-b tag
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2

  return f
# W boson: DSCB

def fitWboson(fname,systname,mmin,mmax,wbin,cand) :

  print 'Systematic: ', systname 
  
  f = TF1(fname,doublecrystal,mmin,mmax,7)
  
  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,4.56374e+02)    # norm
  f.SetParameter(1,7.98614e+01)    # mean
  f.SetParameter(2,9.82938e+00)    # sigma
  f.SetParameter(3,1.60292e+00)    # a1
  f.SetParameter(4,1.94456e-08)    # p1
  f.SetParameter(5,1.35934e+00)    # a2
  f.SetParameter(6,6.64446e-01)    # p2
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,1.42669e+02)    # norm
    f.SetParameter(1,8.05915e+01)    # mean
    f.SetParameter(2,5.00000e+01)    # sigma
    f.SetParameter(3,9.80596e+02)    # a1
    f.SetParameter(4,2.29606e+04)    # p1
    f.SetParameter(5,6.24730e+00)    # a2
    f.SetParameter(6,1.20403e+03)    # p2
  
  #if 'herwig' in systname :
    # initialise sub-peaks in the center since Herwig
    # does not have the high mass stat-fluctuations
    #f.SetParameter(1,8.03790e+01)      # mean1
    #f.SetParameter(7,8.03790e+01)      # mean3

  # ------------------------------------------
  # constraints
  
  M = 80.379
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
  
  return f


# End of: fitWboson()
# ==========================================

# ttbar: DSCB

def fit_ttbar(fname,systname,mmin,mmax,wbin,cand) :

  print 'Systematic: ', systname 

  f = TF1(fname,doublecrystal,mmin,mmax,7)
  
  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,1.48757e+03)    # norm
  f.SetParameter(1,1.79555e+02)    # mean
  f.SetParameter(2,1.33839e+01)    # sigma
  f.SetParameter(3,3.96075e-01)    # a1
  f.SetParameter(4,3.70981e+00)    # p1
  f.SetParameter(5,1.19178e+00)    # a2
  f.SetParameter(6,6.22289e-01)    # p2
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,1.08837e+03)    # norm
    f.SetParameter(1,1.55700e+02)    # mean
    f.SetParameter(2,2.27623e+01)    # sigma
    f.SetParameter(3,3.34465e-01)    # a1
    f.SetParameter(4,1.00000e+05)    # p1
    f.SetParameter(5,1.08504e+00)    # a2
    f.SetParameter(6,5.00240e-01)    # p2
  
  # ------------------------------------------
  # constraints
  
  M = 173.000
  #M = 80.379
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50  )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
  
  return f

def fit_ttbar_cat2(fname,systname,mmin,mmax,wbin,cand) :
 
  print 'Systematic: ', systname 
  
  f = TF1(fname,doublecrystal,mmin,mmax,7)

  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,6.75089e+01)    # norm
  f.SetParameter(1,1.24151e+02)    # mean
  f.SetParameter(2,9.26263e+00)    # sigma
  f.SetParameter(3,6.98991e-01)    # a1
  f.SetParameter(4,9.99886e+04)    # p1
  f.SetParameter(5,2.06096e+00)    # a2
  f.SetParameter(6,8.00505e-01)    # p2
  f.FixParameter(7,0)
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,4.26329e+01)    # norm
    f.SetParameter(1,1.14136e+02)    # mean
    f.SetParameter(2,1.26103e+01)    # sigma
    f.SetParameter(3,6.70565e-01)    # a1
    f.SetParameter(4,9.99993e+04)    # p1
    f.SetParameter(5,2.05157e+00)    # a2
    f.SetParameter(6,6.82528e-01)    # p2
    f.FixParameter(7,0)

  # ------------------------------------------
  # constraints
  
  M = 125.18 
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
 
  return f  

# End of: fit_ttbar()
def fit_ttbar_cat2W(fname,systname,mmin,mmax,wbin,cand) :
 
  print 'Systematic: ', systname 
  
  f = TF1(fname,doublecrystal,mmin,mmax,7)
  
  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,4.56374e+00)    # norm
  f.SetParameter(1,7.98614e+01)    # mean
  f.SetParameter(2,9.82938e+00)    # sigma
  f.SetParameter(3,1.60292e+00)    # a1
  f.SetParameter(4,1.94456e-08)    # p1
  f.SetParameter(5,1.35934e+00)    # a2
  f.SetParameter(6,6.64446e-01)    # p2
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,1.42669e+02)    # norm
    f.SetParameter(1,8.05915e+01)    # mean
    f.SetParameter(2,5.00000e+01)    # sigma
    f.SetParameter(3,9.80596e+02)    # a1
    f.SetParameter(4,2.29606e+04)    # p1
    f.SetParameter(5,6.24730e+00)    # a2
    f.SetParameter(6,1.20403e+03)    # p2
  
  #if 'herwig' in systname :
    # initialise sub-peaks in the center since Herwig
    # does not have the high mass stat-fluctuations
    #f.SetParameter(1,8.03790e+01)      # mean1
    #f.SetParameter(7,8.03790e+01)      # mean3

  # ------------------------------------------
  # constraints
  
  M = 80.379
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
  
  return f

def fit_ttbar_cat2O(fname,systname,mmin,mmax,wbin,cand) :
 
  #print 'Systematic: ', systname 
  print 'Systematic: ', systname 
  
  f = TF1(fname,doublecrystal,mmin,mmax,7)

  # ------------------------------------------
  
  # leading (default)
  f.SetParameter(0,6.75089e+01)    # norm
  f.SetParameter(1,1.24151e+02)    # mean
  f.SetParameter(2,9.26263e+00)    # sigma
  f.SetParameter(3,6.98991e-01)    # a1
  f.SetParameter(4,9.99886e+04)    # p1
  f.SetParameter(5,2.06096e+00)    # a2
  f.SetParameter(6,8.00505e-01)    # p2
  f.FixParameter(7,0)
  
  # subleading
  if cand == 's' :
    f.SetParameter(0,4.26329e+01)    # norm
    f.SetParameter(1,1.14136e+02)    # mean
    f.SetParameter(2,1.26103e+01)    # sigma
    f.SetParameter(3,6.70565e-01)    # a1
    f.SetParameter(4,9.99993e+04)    # p1
    f.SetParameter(5,2.05157e+00)    # a2
    f.SetParameter(6,6.82528e-01)    # p2
    f.FixParameter(7,0)

  # ------------------------------------------
  # constraints
  
  M = 125.18 
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
 
  return f
    
# fit_dijets
def fit_dijets(fname,systname,mmin,mmax,wbin,cand) :
  
  print 'Systematic: ', systname
  #f = TF1(fname,"TMath::Exp([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x )",mmin,mmax);
  #f = TF1(fname,"TMath::Exp([0] + [1]*x )",mmin,mmax);
  #f = TF1(fname,"TMath::Exp([0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x)",mmin,mmax); 
  f = TF1(fname,"TMath::Exp([0] + [1]*x + [2]*x*x + [3]*x*x*x)",mmin,mmax); 
  f.SetParameter(0,1)
 # f.SetParameter(1,10)
 # f.SetParameter(2,10)
 # f.SetParameter(3,10)
 # f.SetParameter(4,10)
  return f


# ==========================================
# CONFNOTE 
# ==========================================
# Higgs, conf, note : 3 gaussians

def fitHiggs_confnote(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 
  f = TF1(fname,"[0]*(TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",mmin,mmax); 
 
  # best fit values using 90-110% constraints
  f.SetParameter(0,3.34620e+02)      # norm
  f.SetParameter(1,1.13398e+02)      # mean1
  f.SetParameter(2,1.53517e+01)      # sigma1
  f.SetParameter(3,1.63080e+00)      # N2/norm
  f.SetParameter(4,1.24809e+02)      # mean2
  f.SetParameter(5,8.29813e+00)      # sigma2
  f.SetParameter(6,1.46720e-01)      # N3/norm
  f.SetParameter(7,1.26176e+02)      # mean3
  f.SetParameter(8,5.00000e+01)      # sigma3
  f.SetParameter(9,2.55738e-01)      # const.


  # ------------------------------------------
  # constraints
  
  M = 125.18 
  
  f.SetParLimits(0, 0.    , 1e+05  )      # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )      # mean1
  f.SetParLimits(2, wbin  , 50.   )      # sigma1
  f.SetParLimits(3, 0.    , 1e+05  )      # N2/norm
  f.SetParLimits(4, 0.9*M , 1.1*M )      # mean2
  f.SetParLimits(5, wbin  , 50.   )      # sigma2
  f.SetParLimits(6, 0.    , 1e+05  )      # N3/norm
  f.SetParLimits(7, 0.9*M , 1.1*M )      # mean3
  f.SetParLimits(8, wbin  , 50.   )      # sigma3
  f.SetParLimits(9, 0.    , 1e+05  )      # const.

  return f

# End of: fitHiggs_confnote()
# ==========================================

# Z boson, conf. note : 3 gaussians

def fitZboson_confnote(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 
  f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",mmin,mmax);
 
  # best fit values with 85-115% constraints on two out of 3 mu
  f.SetParameter(0,3.10995e+03)      # norm
  f.SetParameter(1,8.20689e+01)      # mean1 
  f.SetParameter(2,9.23617e+00)      # sigma1
  f.SetParameter(3,5.44394e+00)      # N2/norm
  f.SetParameter(4,9.43995e+01)      # mean2
  f.SetParameter(5,8.33186e+00)      # sigma2
  f.SetParameter(6,1.63850e+00)      # N3/norm
  f.SetParameter(7,1.00305e+02)      # mean3
  f.SetParameter(8,4.93243e+01)      # sigma3
  f.SetParameter(9,2.04081e+01)      # const.
  

  # ------------------------------------------
  # constraints
  
  M = 91.1876 
  
  f.SetParLimits(0, 0.    , 1e+05  )      # norm
  f.SetParLimits(1, 0.85*M , 1.3*M )      # mean1
  f.SetParLimits(2, wbin  , 50.   )      # sigma1
  f.SetParLimits(3, 0.    , 1e+05  )      # N2/norm
  f.SetParLimits(4, 0.9*M , 1.1*M )      # mean2
  f.SetParLimits(5, wbin  , 50.   )      # sigma2
  f.SetParLimits(6, 0.    , 1e+05  )      # N3/norm
  f.SetParLimits(7, 0.9*M , 1.1*M )      # mean3
  f.SetParLimits(8, wbin  , 50.   )      # sigma3
  f.SetParLimits(9, 0.    , 1e+05  )      # const.

  return f


# End of: fitZboson_confnote()
# ==========================================

# W boson, conf. note : 3 gaussians

def fitWboson_confnote(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 
  f = TF1(fname,"[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",mmin,mmax);
  
  
  # best fit values with 85-115% constraints on two out of 3 mu
  f.SetParameter(0,1.88833e+03)      # norm
  f.SetParameter(1,7.23412e+01)      # mean1
  f.SetParameter(2,2.40493e+01)      # sigma1
  f.SetParameter(3,7.05691e-01)      # N2/norm
  f.SetParameter(4,8.02729e+01)      # mean2
  f.SetParameter(5,5.00000e+00)      # sigma2
  f.SetParameter(6,4.88712e-01)      # N3/norm
  f.SetParameter(7,8.84169e+01)      # mean3
  f.SetParameter(8,2.06415e+01)      # sigma3
  f.SetParameter(9,2.02217e-01)      # const.

  if 'herwig' in systname :
    # initialise sub-peaks in the center since Herwig
    # does not have the high mass stat-fluctuations
    f.SetParameter(1,8.03790e+01)      # mean1
    f.SetParameter(7,8.03790e+01)      # mean3



  # ------------------------------------------
  # constraints
  
  M = 80.379
  
  f.SetParLimits(0, 0.    , 1e+05  )      # norm
  f.SetParLimits(1, 0.95*M , 1.0*M )      # mean1
  f.SetParLimits(2, wbin  , 50.   )      # sigma1
  f.SetParLimits(3, 0.    , 1e+05  )      # N2/norm
  f.SetParLimits(4, 0.9*M , 1.1*M )      # mean2
  f.SetParLimits(5, wbin  , 50.   )      # sigma2
  f.SetParLimits(6, 0.    , 1e+05  )      # N3/norm
  f.SetParLimits(7, 0.9*M , 1.1*M )      # mean3
  f.SetParLimits(8, wbin  , 50.   )      # sigma3
  f.SetParLimits(9, 0.    , 1e+05  )      # const.
  
  return f


# End of: fitWboson_confnote()
# ==========================================

# ttbar, conf. note : DSCB

def fit_ttbar_confnote(fname,systname,mmin,mmax,wbin) :

  print 'Systematic: ', systname 

  f = TF1(fname,doublecrystal,mmin,mmax,7)
  
  # best fit values with 85-115% constraints on two out of 3 mu
  f.SetParameter(0,6.34917e+02)    # norm
  f.SetParameter(1,1.71817e+02)    # mean
  f.SetParameter(2,1.46715e+01)    # sigma
  f.SetParameter(3,3.91665e-01)    # a1
  f.SetParameter(4,5.21391e+00)    # p1
  f.SetParameter(5,1.20031e+00)    # a2
  f.SetParameter(6,3.56933e-01)    # p2
  f.SetParameter(7,1.0)
  # ------------------------------------------
  # constraints
  
  M = 173.0
  
  f.SetParLimits(0, 0.    , 1e+05 )    # norm
  f.SetParLimits(1, 0.9*M , 1.1*M )    # mean
  f.SetParLimits(2, wbin  , 50.   )    # sigma
  f.SetParLimits(3, 0.    , 1e+05 )    # a1
  f.SetParLimits(4, 0.    , 1e+05 )    # p1
  f.SetParLimits(5, 0.    , 1e+05 )    # a2
  f.SetParLimits(6, 0.    , 1e+05 )    # p2
  
  return f

# End of: fit_ttbar_confnote()
# ==========================================
