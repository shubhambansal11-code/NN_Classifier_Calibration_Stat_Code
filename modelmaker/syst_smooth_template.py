# ==========================================
# Systematic Smoothing for XML Model Making
#
# Original author: Migle Stankaityte
# 4 February 2020

#Refurbished for QT PhD 2019-2020
# ==========================================

#import
import sys, re

import sigfit_init as init

import ROOT
from ROOT import gROOT, gStyle, TH1F, TF1, TCanvas, TFile, TAxis
from ROOT import TSpline3

sys.path.insert(0, 'plotting/')
import plot_helpers

gROOT.SetBatch(True)


# ==========================================
# TODO add wider plots!!

def getSmoothedSysts(template,rebin,outpath,nbins,mmin,mmax): 

  #------------ coarse plots -------------
  # nominal
  hname = template['hist'].GetTitle()
  template['h_coarse'] = template['hist'].Clone(hname+'_coarse')
  template['h_coarse'].SetName(hname+'_coarse') 
  template['h_coarse'].Rebin(rebin)

  # systematics 
  template['hsyst_up_coarse'] = {}
  template['hsyst_down_coarse'] = {}

  # syst/nom ratios
  hup_ratio = {}
  hdown_ratio = {}

  # splines
  #sup_ratio = {}
  #sdown_ratio = {}
  sup = {}
  sdown = {}
  
  # smoothed / fine bin histos
  #hup_rsmooth = {}
  #hdown_rsmooth = {}
  
  hup_smooth = {}
  hdown_smooth = {}

  for syst in template['hsyst_up'] :
    
    hup_name = template['hsyst_up'][syst].GetTitle()
    hdown_name = template['hsyst_down'][syst].GetTitle()

    template['hsyst_up_coarse'][syst] = template['hsyst_up'][syst].Clone(hup_name+'_coarse')
    template['hsyst_down_coarse'][syst] = template['hsyst_down'][syst].Clone(hdown_name+'_coarse')
    
    template['hsyst_up_coarse'][syst].SetName(hup_name+'_coarse') 
    template['hsyst_down_coarse'][syst].SetName(hdown_name+'_coarse') 
    
    template['hsyst_up_coarse'][syst].Rebin(rebin)
    template['hsyst_down_coarse'][syst].Rebin(rebin)

    # ratios
    #hup_ratio[syst] = template['hsyst_up_coarse'][syst].Clone(hup_name+'_ratio')
    #hdown_ratio[syst] = template['hsyst_down_coarse'][syst].Clone(hdown_name+'_ratio')
    
    #hup_ratio[syst].SetName(hup_name+'_ratio') 
    #hdown_ratio[syst].SetName(hdown_name+'_ratio') 

    #hup_ratio[syst].Divide(template['h_coarse'])
    #hdown_ratio[syst].Divide(template['h_coarse'])

    # splines
    #sup_ratio[syst] = TSpline3(hup_ratio[syst])
    #sdown_ratio[syst] = TSpline3(hdown_ratio[syst])
    
    sup[syst] = TSpline3(template['hsyst_up_coarse'][syst])
    sdown[syst] = TSpline3(template['hsyst_down_coarse'][syst])
    
    # initialising histos for smoothed ratios and temps
    #hup_rsmooth[syst] = TH1F(hup_name+'_rsmooth',hup_name+'_rsmooth',nbins,mmin,mmax)
    #hdown_rsmooth[syst] = TH1F(hdown_name+'_rsmooth',hdown_name+'_rsmooth',nbins,mmin,mmax)
    
    hup_smooth[syst] = TH1F(hup_name+'_smooth',hup_name+'_smooth',nbins,mmin,mmax)
    hdown_smooth[syst] = TH1F(hdown_name+'_smooth',hdown_name+'_smooth',nbins,mmin,mmax)

    # smooth systematics (initially clones of nominal)
    #hup_smooth[syst] = template['hist'].Clone(hup_name)
    #hdown_smooth[syst] = template['hist'].Clone(hdown_name)
    
    #hup_smooth[syst].SetName(hup_name) 
    #hdown_smooth[syst].SetName(hdown_name) 

  #---------- draw slplines -----------
  
  for syst in template['hsyst_up'] : 
    #drawSpline(hup_ratio[syst],sup_ratio[syst],outpath)
    #drawSpline(hdown_ratio[syst],sdown_ratio[syst],outpath)
    drawSpline(template['hsyst_up_coarse'][syst],sup[syst],outpath)
    drawSpline(template['hsyst_down_coarse'][syst],sdown[syst],outpath)

  #-------- smooth systematics --------
  
  for syst in template['hsyst_up'] :
    #makeSmoothRatio(hup_rsmooth[syst],sup_ratio[syst],nbins)
    #makeSmoothRatio(hdown_rsmooth[syst],sdown_ratio[syst],nbins)
    makeSmoothRatio(hup_smooth[syst],sup[syst],nbins,template['hsyst_up'][syst].Integral())
    makeSmoothRatio(hdown_smooth[syst],sdown[syst],nbins,template['hsyst_down'][syst].Integral())

    #hup_smooth[syst].Multiply(hup_rsmooth[syst])   
    #hdown_smooth[syst].Multiply(hdown_rsmooth[syst])   

  # TODO: plot
  #------------- root output -------------
  # making output file
  fout = TFile(outpath+hname+"_systsmooth.root","RECREATE")
  
  # --- nominal ---
  template['hist'].Write()
  template['h_coarse'].Write()
  
  # --- systematics ---
  for syst in template['hsyst_up'] :
   
    # up
    template['hsyst_up'][syst].Write()
    template['hsyst_up_coarse'][syst].Write()
    #hup_ratio[syst].Write()
    #sup_ratio[syst].Write()
    #hup_rsmooth[syst].Write()
    hup_smooth[syst].Write()

    # down
    template['hsyst_down'][syst].Write()
    template['hsyst_down_coarse'][syst].Write()
    #hdown_ratio[syst].Write()
    #sdown_ratio[syst].Write()
    #hdown_rsmooth[syst].Write()
    hdown_smooth[syst].Write()

  # closing
  fout.Close()

  sys.exit()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def drawSpline(h,spline,outpath) :

  # gStyle Options 
  gStyle.SetOptStat(0)
 
  # title
  hname = h.GetTitle()

  # legend
  #pos =[]
  pos=[0.65,0.75,0.85,0.87]

  l = plot_helpers.makeLegend([h,spline],
                              ['template','spline'],
                              pos=pos
                              #ncols=2
                              )

  # cosmetics 
  #h.SetMaximum(h.GetMaximum()+((h.GetMaximum()-h.GetMinimum())*0.75))
  h.SetMaximum(h.GetMaximum()*1.2)
  h.GetXaxis().SetTitle('Large-R Jet Mass [GeV]')
  h.GetYaxis().SetTitle('Events / 5 GeV')
  
  h.GetXaxis().SetLabelSize(0.045)
  h.GetYaxis().SetLabelSize(0.045)
  h.GetXaxis().SetTitleSize(0.045)
  h.GetYaxis().SetTitleSize(0.045)

  # making plot
  c = TCanvas("c"+hname,"c"+hname,1600,1000) 
  h.Draw()
  spline.Draw('samec')
  l.Draw()

  # atlas label
  plot_helpers.drawAtlasLabel(0.13,0.85,
                              p_type="Preliminary Internal",
                              s="13TeV",
                              lumi="136fb^{-1}")

  # printing plot
  c.Print(outpath+'spline_'+hname+".pdf")


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def makeSmoothRatio(h,spline,nbins,integral):
  # looping over bins
  ibin = 1  
  while ibin<nbins+1 :

    w = h.GetBinWidth(ibin)
    x = h.GetBinLowEdge(ibin)+(w/2)

    # filling histogram with mid values of spline
    y  = spline.Eval(x)
    #print('coords:\t',x,' ',y)

    h.SetBinContent(ibin,y)
    ibin += 1

  # rescaling to match the integral of original histogram
  h.Scale(integral/h.Integral())





