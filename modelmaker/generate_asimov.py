# Original Author Francesco Maria Follega, 19th April 2019
# Functions to generate asimov data for H,V,ttbar and QCD
# Refurbished for QT PhD 2019-2020

import sys
from ROOT import TF1,TFile,TH1F,TMath,TRandom3

class fnc_dscb:
    def __call__(self, xx, pp):
      x   = xx[0]
      N   = pp[0]
      mu  = pp[1]
      sig = pp[2]
      a1  = pp[3]
      p1  = pp[4]
      a2  = pp[5]
      p2  = pp[6]
      
      t = (x-mu)/sig
     
      if (t < -a1):
        a = TMath.Exp(-0.5*a1*a1)
        b = p1/a1 - a1
        return N*a/TMath.Power(a1/p1*(b - t), p1)
      elif (t > a2 ):
        a = TMath.Exp(-0.5*a2*a2)
        b = p2/a2 - a2
        return N*a/TMath.Power(a2/p2*(b + t), p2)
      return N*TMath.Exp(-0.5*t*t)

def generate_asimov_H(N,_min,_max):
  faH = TF1("faH","[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",_min,_max)
  faH.FixParameter(0, 3.24618e+02) # norm
  faH.FixParameter(1, 1.15838e+02) # mean1
  faH.FixParameter(2, 1.50755e+01) # sigma1
  faH.FixParameter(3, 1.30503e+00) # N2/norm
  faH.FixParameter(4, 1.25840e+02) # mean2
  faH.FixParameter(5, 8.33481e+00) # sigma2
  faH.FixParameter(6, 6.02804e+02) # N3/norm
  faH.FixParameter(7, -5.19023e+03) # mean3
  faH.FixParameter(8, 2.32496e+03) # sigma3
  faH.FixParameter(9, -1.95470e+00) # Const.

  SignalInjected = TH1F("SignalInjected","SignalInjected",N,_min,_max)

  for _bin in range(1,N+1):

    x1 = SignalInjected.GetBinLowEdge(_bin)
    x2 = x1 + SignalInjected.GetBinWidth(_bin)
    eventsGenerated = faH.Integral(x1,x2)
    SignalInjected.SetBinContent(_bin,eventsGenerated)
  SignalInjected.Scale(216/SignalInjected.Integral())

  return SignalInjected

def generate_asimov_V(N,_min,_max):
  faZ =  TF1("faZ","[0]* (TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Gaus(x,[7],[8],1)) + [9]",_min,_max)
  faZ.FixParameter(0, 1.98793e+04) #  norm
  faZ.FixParameter(1, 9.26091e+01) #  mean1
  faZ.FixParameter(2, 9.16318e+00) #  sigma1
  faZ.FixParameter(3, 1.64196e+00) #  N2/norm
  faZ.FixParameter(4, 1.65687e+01) #  mean2
  faZ.FixParameter(5, 6.21259e+01) #  sigma2
  faZ.FixParameter(6, 3.06609e+02) #  N3/norm
  faZ.FixParameter(7, -1.22241e+04) #  mean3
  faZ.FixParameter(8, 9.45498e+03) #  sigma3
  faZ.FixParameter(9, -8.33126e+01) # Const.

  ZInjected = TH1F("ZInjected","ZInjected",N,_min,_max)
  
  for _bin in range(1,N+1):

    x1 = ZInjected.GetBinLowEdge(_bin)
    x2 = x1 + ZInjected.GetBinWidth(_bin)
    eventsGenerated = faZ.Integral(x1,x2)
    ZInjected.SetBinContent(_bin,eventsGenerated)
  ZInjected.Scale(7700/ZInjected.Integral())

  return ZInjected

def generate_asimov_ttbar(N,_min,_max):
  fattbar = TF1("fattbar",fnc_dscb(),_min,_max,7);
  fattbar.FixParameter(0, 6.37124e+02) # norm
  fattbar.FixParameter(1, 1.71570e+02) # mean
  fattbar.FixParameter(2, 1.48571e+01) # sigma
  fattbar.FixParameter(3, 4.23357e-01) # a1
  fattbar.FixParameter(4, 3.51360e+00) # p1
  fattbar.FixParameter(5, 1.20526e+00) # a2
  fattbar.FixParameter(6, 3.57687e-01) # p2
  # bins for template  
  ttbarInjected = TH1F("ttbarInjected","ttbarInjected",N,_min,_max)
  
  for _bin in range(1,N+1):
    x1 = ttbarInjected.GetBinLowEdge(_bin)
    x2 = x1 + ttbarInjected.GetBinWidth(_bin)
    eventsGenerated = fattbar.Integral(x1,x2)
    ttbarInjected.SetBinContent(_bin,eventsGenerated)
  ttbarInjected.Scale(10550/ttbarInjected.Integral())
  return ttbarInjected

def generate_asimov_QCD(N,_min,_max,param):

  assert(len(param)>=2)

  # making function string
  f_qcd = "[0]*TMath::Exp("
  for i,p in enumerate(param):   
    if not i : continue
    if i == 1 : 
      f_qcd = f_qcd + "([1]*(x-140.)/70.)"
    else :
      f_qcd = f_qcd + "+(["+str(i)+"]*pow((x-140.)/70.,"+str(i)+"))"
  f_qcd = f_qcd + ")"
  
  # print to screen
  print "........"
  print "Asimov QCD:"
  print f_qcd
  print "........"

  # QCD parameteric
  faQCD = TF1("faQCD",f_qcd,_min,_max)
  
  # setting parameters
  for i,p in enumerate(param):   
    if i == 0: faQCD.FixParameter(i,1)
    else : faQCD.FixParameter(i,p)
  
  # adding yield
  norm_par = param[0]/faQCD.Integral(_min,_max)
  faQCD.FixParameter(0,norm_par)

  BestBackgroundPrediction = TH1F("BestBackgroundPrediction","BestBackgroundPrediction",N,_min,_max)

  for _bin in range(1,N+1):

    x1 = BestBackgroundPrediction.GetBinLowEdge(_bin)
    x2 = x1 + BestBackgroundPrediction.GetBinWidth(_bin)
    eventsGenerated = (faQCD.Integral(x1,x2))
    BestBackgroundPrediction.SetBinContent(_bin,eventsGenerated)

  return BestBackgroundPrediction
