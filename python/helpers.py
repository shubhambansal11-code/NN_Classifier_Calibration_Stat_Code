# ==========================================
# XML Model Creation Helpers
#
# Original Author: Migle Stankaityte
# 30 June 2019
# Refurbished for QT PhD project 2019-2020
# ==========================================

#import
import sys, os, re
import json

import numpy as np

import ROOT
from ROOT import gROOT, TFile, TTree, TH1F, TLorentzVector
gROOT.SetBatch(True)

# ==========================================
# Making Directories
# void 
def mkdir(path):

  print 'Creating Directory: ',path
  if not os.path.exists(path):
    os.mkdir(path)
    print 'Directory ',path,' created.\n'
  else:    
    print 'Directory ',path,' already exists.\n'
  return

# ==========================================
# Getting Printable Parameter Name

def getPrettyParName(name,latex=False) :

  if 'alpha_' in name : 
    name = re.sub('^alpha_','',name)
    name = re.sub('_',' ',name)
  elif 'gamma_stat_' in name : 
    name = re.sub('^gamma_stat_','stat error: ',name)
    name = re.sub('_',' ',name)
    name = re.sub('bin','(bin',name)
    name = re.sub('$',')',name)
  
  elif 'dnll' in name :
    if latex : name = re.sub('dnll','2#times#Delta NNL',name)
    else : name = re.sub('dnll','2 Delta NLL ',name)
  elif 'nll' in name :
    name = re.sub('nll','NLL',name)
  
  elif latex and 'mu' in name :  
    name = re.sub('mu','#mu',name)
    if '_' in name :
      name = re.sub('_','(',name)
      name = re.sub('$',')',name)
  else:
    name = re.sub('_',' ',name)
    if latex : 
      if 'ttbar' in name :  
        name = re.sub('ttbar','t#bar{t}',name)

  return name

# ==========================================
# Getting Printable Systematic Name

def getPrettySystName(name,tname,latex=False) :

  name = re.sub('_',' ',name)
  if tname in name : 
    name = re.sub(tname,tname+' ',name)
  if '1down' in name : 
    name = re.sub('1down','(Down Var.)',name)
  if '1up' in name : 
    name = re.sub('1up','(Up Var.)',name)

  return name

# ==========================================

def round_down(n, decimals=0):
  multiplier = 10 ** decimals
  return math.floor(n * multiplier) / multiplier

