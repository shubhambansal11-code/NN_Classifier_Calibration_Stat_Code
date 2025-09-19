# ==========================================
# Subtracted data cross-check plot
#
# Original Author: Migle Stankaityte
# June 2019
# Refurbished for QT PhD project 2019-2020
# ==========================================

#import
import sys, os, re
#import json
import argparse
import ROOT
from ROOT import gROOT, TFile, TTree, TH2D, TCanvas, THStack, TLegend, TText
from ROOT import kBlack, gStyle, gPad
from ROOT import RooFitResult, RooArgSet, RooArgList, RooDataHist
gROOT.SetBatch(True)

import helpers

sys.path.insert(0, 'plotting/')
import plot_helpers

gStyle.SetOptStat(0)
# ==========================================

# salutation
print 'Plotting fit to data.'

# getting argumets

parser = argparse.ArgumentParser(description='Printing fit-result information')
parser.add_argument('infile', help='Path of input file')
parser.add_argument('--debug', help='debug', action='store_true')
parser.add_argument('-g','--gammas', help='Show gamma pars', action='store_true')
parser.add_argument('-d','--details', help='details', action='store_true')
parser.add_argument('-t','--text', help='text in corr plot', action='store_true')
parser.add_argument('-r','--reverse', help='reverse', action='store_true')
parser.add_argument('--tag', default=None, help='tag')
args = parser.parse_args()

# debug
debug = args.debug

# gammas
gammas = args.gammas

# detailed
details = args.details

# text in plot
text = args.text

# tag
tag = args.tag

# getting output name
outname = args.infile[args.infile.rfind('/')+1:args.infile.rfind('.root')]

# -------------------------------

# reading file
f = TFile(args.infile,'READ')

# -------------------------------
# get workspace details mode

if details :
  fit = f.Get('combWS')
  fit.Print()
  sys.exit()

# -------------------------------

# getting fit result
fres = f.Get('fitResult')

# print fit result
print '*****************************************'
print 'Fit Result: '
fres.Print()
print '*****************************************'

# -------------------------------
# getting parameter list

# all parameters
init_pars = fres.floatParsFinal()
if debug : 
  print 'All parameters:'
  init_pars.Print()
  print

# initialising plotted paramater set
pars_set = RooArgSet("param_list")

# looping over all parameters
pars_exist = True
i = 0
npars = 0

while pars_exist :
  par = init_pars.at(i)
  
  if par : 
    i += 1
    par_t = par.GetTitle()

    # skipping gammas ( -g to keep ) 
    if not gammas and "gamma_stat" in par_t :
      if debug : print "-- Skipping ",par_t
      continue

    if debug: print "++ Saving ",par_t
    pars_set.add(par)
    npars += 1
  
  # exiting
  else : pars_exist = False

# debug : printing par list
if debug : 
  print "\nParameter set (",npars,"/",i,") for Corr Matrix:"
  pars_set.Print()
  print

# convert to list
pars = RooArgList(pars_set)
if debug :
  print "Converting set to list: ",pars

# reverse
if args.reverse: pars.sort(True)

# -------------------------------
# making plot

h_corr = TH2D("corrhist",outname,npars,0,npars,npars,0,npars)
p_titles = []

# filling
for i in range(npars):
  
  p1 = pars.at(i)
  
  # saving title
  p1_t = helpers.getPrettyParName(p1.GetTitle(),latex=True)
  p_titles.append(p1_t)
  
  for j in range(npars):
    
    p2 = pars.at(j)

    # getting correlation
    corr = fres.correlation(p1,p2)
    
    # print correlations
    if debug:
      p1_tp = helpers.getPrettyParName(p1.GetTitle())
      p2_tp = helpers.getPrettyParName(p2.GetTitle())
      print "Correlation : ",p1_tp,' & ',p2_tp,": ",corr
    
    # fill histo
    h_corr.SetBinContent(i+1,j+1,corr)

# settings
h_corr.SetMinimum(-1)
h_corr.SetMaximum(1)

# space
print
# ------ drawing ------

# canvas
if '_inc_' in outname or '_fid_' in outname : c_ch = TCanvas("c_ch",outname, 2000, 1200)
elif 'confNote' in outname : c_ch = TCanvas("c_ch",outname, 800, 600)
elif 'ttbar' in outname : c_ch = TCanvas("c_ch",outname, 1200, 800)
else: c_ch = TCanvas("c_ch",outname, 2000, 900)

# plot cosmetics
gStyle.SetPaintTextFormat('4.2f')
gStyle.SetPalette(57)

if 'ttbar' in outname : gPad.SetLeftMargin(0.2)
else: gPad.SetLeftMargin(0.15)

if '_inc_' in outname or '_fid_' in outname : gPad.SetBottomMargin(0.25)
elif 'confNote' in outname : gPad.SetBottomMargin(0.2)
else : gPad.SetBottomMargin(0.3)

gPad.SetTopMargin(0.08)

# draw 
if text : h_corr.Draw('colztext')
else: h_corr.Draw('colz')

# add labels

if '_inc_' in outname or '_fid_' in outname : plot_helpers.addLabels(h_corr,npars,p_titles,npars,p_titles,ts=0.015)
elif 'confNote' in outname : plot_helpers.addLabels(h_corr,npars,p_titles,npars,p_titles,ts=0.02)
else: plot_helpers.addLabels(h_corr,npars,p_titles,npars,p_titles)

# print
helpers.mkdir('plots/fitresults/')
if gammas : outname = outname+"_gammas"
if tag: outname = outname+"_"+tag
helpers.mkdir('plots/fitresults/'+outname)
c_ch.Print('plots/fitresults/correlations/'+outname+'_correlation.pdf')

sys.exit()



