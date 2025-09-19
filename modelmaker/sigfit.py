# ==========================================
# Signal Fitting Lib for XML Model Making
#
# Migle Stankaityte
# 1 April 2019

#modified for Z->bb calibration
#QT PhD Project 2019-2020
# ==========================================

#import
import sys, re

import sigfit_init as init

import ROOT
from ROOT import gROOT, gStyle, TH1F, TF1, TCanvas, TFile, TAxis

gROOT.SetBatch(True)
#gStyle.SetOptStat(0)
#gStyle.SetPalette(1)
#gROOT.LoadMacro("../style/AtlasStyle.C")
#gROOT.LoadMacro("../style/AtlasUtils.C")
#SetAtlasStyle()

# ==========================================
# Reading json file with model information
# returns: histograms produced with parametric fitting

def getFittedHistos(template,h,h_up,h_down,wbin_fit,nbins_fit,mmin_fit,mmax_fit,nbins_out,mmin_out,mmax_out,ccpath,conf,cand) :

  print 'Running template fiting to parametric functions.'

  #------------- initialising -------------
  # dh - dummy for fit
  # ph - final pfit histo

  h_name = re.sub(r'\_orig','',h.GetTitle())
  dh = h.Clone(h_name+'_dummy')
  ph = TH1F(h_name,h_name,nbins_out,mmin_out,mmax_out)
  
  dh_up = {}
  ph_up = {}
  dh_down = {}
  ph_down = {}
  
  for syst in h_up :

    h_up_name = re.sub(r'\_orig','',h_up[syst].GetTitle())
    h_down_name = re.sub(r'\_orig','',h_down[syst].GetTitle())

    dh_up[syst] = h_up[syst].Clone(h_up_name+'_dummy')
    ph_up[syst] = TH1F(h_up_name,h_up_name,nbins_out,mmin_out,mmax_out)
    
    dh_down[syst] = h_down[syst].Clone(h_down_name+'_dummy')
    ph_down[syst] = TH1F(h_down_name,h_down_name,nbins_out,mmin_out,mmax_out)
  

  # ------------ fitting procedure -------------
  
  # --- nominal ---
  # getting fit function
  fname = 'f_'+h_name
  f = init.getFitFn(template,fname,mmin_fit,mmax_fit,wbin_fit,conf,cand)
  
  # fitting
  doFit(f,fname,dh,template,nbins_fit,mmin_fit,mmax_fit,ccpath)

  
  # --- systematics ---
  f_up = {}
  f_down = {}

  for syst in ph_up :
    
    fname_up= 'f_'+ph_up[syst].GetTitle()
    fname_down = 'f_'+ph_down[syst].GetTitle()
    
    # getting fit function
    f_up[syst] = init.getFitFn(template,fname_up,mmin_fit,mmax_fit,wbin_fit,conf,cand) 
    f_down[syst] = init.getFitFn(template,fname_down,mmin_fit,mmax_fit,wbin_fit,conf,cand)
    
    # fitting
    doFit(f_up[syst],fname_up,dh_up[syst],template,nbins_fit,mmin_fit,mmax_fit,ccpath)
    doFit(f_down[syst],fname_down,dh_down[syst],template,nbins_fit,mmin_fit,mmax_fit,ccpath)
    

  #------------- making parametric histograms -------------
  
  # --- nominal ---
  # making/writing histo
  makeFitHist(ph,f,h.Integral(),nbins_out)

  # --- systematics ---
  for syst in ph_up :  
    # up
    makeFitHist(ph_up[syst],f_up[syst],h_up[syst].Integral(),nbins_out)
    # down
    makeFitHist(ph_down[syst],f_down[syst],h_down[syst].Integral(),nbins_out)
  
  #------------- root output -------------
  # making output file
  fout = TFile(ccpath+h_name+"_fit.root","RECREATE")
  
  # --- nominal ---
  ph.Write()
  f.Write()
  
  # --- systematics ---
  for syst in ph_up :
    
    # up
    ph_up[syst].Write()
    f_up[syst].Write()
    
    # down
    ph_down[syst].Write()
    f_down[syst].Write()

  # closing
  fout.Close()

  #----------------------------------------
  # out: hist/dict/dict
  return ph, ph_up, ph_down

# End of: getFittedHistos
# ==========================================

# Performing fit
def doFit(f,fname,dh,template,nbins,mmin,mmax,ccpath) : 
 
  print '********** FITTING ',(mmax-mmin)/nbins,' GeV **********'
 
  # gStyle Options 
  gStyle.SetOptFit(1111)
  
  if 'ttbar' in template:
    gStyle.SetStatX(0.4)
    gStyle.SetStatY(0.9)
    gStyle.SetStatW(0.15)
    gStyle.SetStatH(0.12)
  else :
    gStyle.SetStatX(0.9)
    gStyle.SetStatY(0.9)
    gStyle.SetStatW(0.15)
    gStyle.SetStatH(0.12)

  # printing fit stats
  dh.Fit(fname)
  fres = dh.GetFunction(fname)
  print 'Chi2:\t',fres.GetChisquare()
  print 'NDF:\t',f.GetNDF()
  print 'Chi2/NDF:\t',fres.GetChisquare()/f.GetNDF()
  print '1-Chi2/NDF:\t',abs(1-fres.GetChisquare()/f.GetNDF()),'\n'
  
  # making fit plots 
  #dh.GetYaxis().SetRangeUser(0.1,dh.GetMaximum()*2)
  dh.GetXaxis().SetTitle('Large-R Jet Mass [GeV]')
  dh.GetYaxis().SetTitle('Events / 1 GeV')
  dh.SetTitle("")
  dh.GetYaxis().SetTitleOffset(1.0)
  dh.GetXaxis().SetTitleOffset(1.0)
  dh.GetXaxis().SetRangeUser(mmin,mmax)
  #TString text_Chi2=Form("#chi2/ndf = %0.3f", fres.GetChisquare()/f.GetNDF())
  c = TCanvas("c_"+fname,"c_"+fname,800,600)
  latex = ROOT.TLatex ()
  latex.SetNDC ()
  latex.SetTextSize (0.1)
  latex.DrawText (0.15 ,0.83 ,"#bf{ATLAS Internal}")
  latex.SetTextSize (0.1)
  latex.DrawText (0.15 ,0.80 ,"Xbb 60% W.P. pT=[450,1000] GeV")
  #latex.SetTextSize (0.06)
  #latex.DrawText (0.7 ,0.77 , "chi2/ndf" = fres.GetChisquare()/f.GetNDF())
  dh.Draw()

  c.Print(ccpath+fname+".pdf")

  print '********** Fit Done! **********'

# End of: doFit()
# ==========================================

# Make histogram from fit

def makeFitHist(ph,f,h_int,nbins) :

  # looping over bins
  ibin = 1  
  while ibin<nbins+1 :
    
    w = ph.GetBinWidth(ibin)
    x0 = ph.GetBinLowEdge(ibin)
    x = x0 + w
    
    # filling histo with fit integral over bin range
    y  = f.Integral(x0,x)
   
   # making sure there are no negative bins
    if y<0 :
      y = 0.00000000001
      print 'WARNING!: in ',ph.GetTitle(),' fit value < 0 in range ',x0,' - ',x,', setting y = 0.00000000001!'

    ph.SetBinContent(ibin,y)
    
    ibin += 1

  # rescaling to match the integral of original histogram
  ph.Scale(h_int/ph.Integral())

# End of: makeFitHist()
# ==========================================


