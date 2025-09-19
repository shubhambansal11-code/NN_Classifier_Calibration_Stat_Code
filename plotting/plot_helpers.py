# ==========================================
# Subtracted data cross-check plot
#
# Original Author: Migle Stankaityte
# 15 May 2019
# Refurbished for QT PhD project 2019-2020
# ==========================================

#import
import sys, os, re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

import ROOT
from ROOT import gROOT, gPad 
from ROOT import TLatex, TText, TFile, TTree, TH1F, TCanvas, THStack, TLegend
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kPink
from ROOT import kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet
from ROOT import kTRUE
gROOT.SetBatch(True)

# ==========================================
def fColors():
    return [kAzure+1,kSpring-5,kGreen+2,kYellow-4,kRed-7,kOrange+1,kGray+1,kViolet-4,kRed-6,kYellow-6]
# ==========================================
def lColors():
    return [kAzure+2,kSpring-6,kGreen+3,kOrange,kOrange+2,kRed+1,kGray+2,kViolet-5,kRed-2,kYellow+2]

# ==========================================
# TH1
def getColorOptions(colors,nhists,fill=False) :
  
  # colors is an index list of colors you want

  if fill :
    cl_o = fColors() 
  else :
    cl_o = lColors()

  cl = []
  
  if colors :
    assert(len(colors) < cl)
    assert(len(colors) == nhists) 
    for ic in colors :
      cl.append(cl_o[ic])
  else :
    cl = cl_o

  return cl

# ==========================================
# TH1
def setHistoColors(hists,fill=True,fillt=False,colors=[]) :

  fillcolors = getColorOptions(colors,len(hists),fill=True) 
  linecolors = getColorOptions(colors,len(hists))
  
  for i,h in enumerate(hists) : 
    h.SetLineColor(linecolors[i])
    h.SetLineWidth(2)
    if fillt : h.SetFillColorAlpha(fillcolors[i], 0.35)
    elif fill : h.SetFillColor(fillcolors[i])
  
  return
# ==========================================
# TH1
def setMarkerColors(hists,colors=[]) :

  markercolors = getColorOptions(colors,len(hists))
  for i,h in enumerate(hists) : 
    h.SetMarkerColor(markercolors[i])
    h.SetMarkerSize(2)
  
  return

# ==========================================
# TH1
def atlasLabel(x,y,color=1,text="",tsize=0.04,mt=0.08):
  l = TLatex()
  l.SetNDC(kTRUE)
  l.SetTextColor(color)
  l.SetTextFont(72)
  l.SetTextSize(tsize)
  l.DrawLatex(x, y, "ATLAS")
  if text:
    atlasText(text,x+mt,y,color=1,tsize=tsize)

# ==========================================
# TH1
def atlasText(text,x,y,color=1,tsize=0.04) :
  l = TLatex()
  l.SetNDC(kTRUE)
  l.SetTextFont(42)
  l.SetTextColor(color)
  l.SetTextSize(tsize)
  l.DrawLatex(x, y, text)

# ==========================================
# TH1
def drawAtlasLabel(x,y,p_type="",s="",lumi="",slabel="",tlabel="",color=1, space=0.05, tsize=0.04, mt=0.08) :
  
  # ATLAS
  atlasLabel(x,y,color=color,
             text=p_type,tsize=tsize,mt=mt)
  
  # info
  inf = ""
  
  if s: 
    inf = inf+"#sqrt{s}="+s
    if lumi: inf = inf+", "
  if lumi: inf = inf+lumi

  y2 = y-space
  y3 = y2-space
  y4 = y3-space

  if inf : 
    atlasText(inf, x, y2,color=color,tsize=tsize)
 
  # sublabel 
  if slabel :
    if inf : atlasText(slabel,x,y3,color=color,tsize=tsize)
    else : atlasText(slabel,x,y2,color=color,tsize=tsize)

    # thirdlabel (only if slabel exists)
    if tlabel :
      if inf : atlasText(tlabel,x,y4,color=color,tsize=tsize)
      else : atlasText(tlabel,x,y3,color=color,tsize=tsize)


# ==========================================
# TH1
def setDataAtt(h, rebin=0) :

  h.SetLineColor(kBlack)
  h.SetMarkerStyle(20)
  if rebin : h.Rebin(rebin)
  
  return

# ==========================================
# TH1
def makeLegend(hists,tags,pos=[],loc='topright',ammend=[], ncols=1, edge=False) :
  
  # getting coordinates
  if pos :
    assert(len(pos)==4),'Legend coordinate number is not 4!'
    loc = None
  
  # converting location name to coordinates
  if loc :
    if loc == 'topright' : 
      pos = [0.55,0.55,0.85,0.85]
    elif loc == 'topleft' : 
      pos = [0.15,0.55,0.45,0.85]
    else : raise ValueError('Unknown legend location name: ',loc,
                              '\nKnown values: topleft, topright')

  # ammending coordinates
  if ammend :
    assert(len(ammend)==4),'Legend coordinate ammendment number is not 4!'
    pos = list(np.array(pos)+np.array(ammend))
    
  x0, y0, x, y = pos
  
  l = TLegend(x0,y0,x,y)

  for i,h in enumerate(hists):
    l.AddEntry(h,tags[i])  
  
  # columns
  l.SetNColumns(ncols)

  # edge
  if not edge :
    l.SetBorderSize(0)

  return l

# ==========================================
# TH2
# making text labels

def addLabels (h,npars_x,pars_x,npars_y,pars_y,ts=0.018) :
  h.GetXaxis().SetLabelOffset(99)
  h.GetYaxis().SetLabelOffset(99)

  #draw labels along X
  y = gPad.GetUymin() - 0.2*h.GetYaxis().GetBinWidth(1)
  t = TLatex() 

  t.SetTextAngle(60)
  t.SetTextSize(float(ts))
  t.SetTextFont(42)
  t.SetTextAlign(33)

  for i in range(npars_x):
    x = h.GetXaxis().GetBinCenter(i+1)
    t.DrawLatex(x,y,pars_x[i])

  # draw labels along y
  x = gPad.GetUxmin() - 0.1*h.GetXaxis().GetBinWidth(1)
  t.SetTextAlign(32)
  t.SetTextAngle(0)
  
  for i in range(npars_y):
    y = h.GetYaxis().GetBinCenter(i+1)
    t.DrawLatex(x,y,pars_y[i])

# ==========================================
# pyplot
def setHistoColorsByHeight(colmap,n,bins,patches) :

  cm = plt.cm.get_cmap(colmap)
  # normalising values
  col = (n-n.min())/(n.max()-n.min())
  
  # setting colors
  for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))
  
  return

# ==========================================
# pyplot
def setTex() :
  rc('text', usetex=True)
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
