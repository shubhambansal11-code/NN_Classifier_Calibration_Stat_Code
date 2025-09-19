# ==========================================
# Helper functions for xml card generation
#
# Original Author: Migle Stankaityte
# 24 September 2019
#
# Modifications made for Z->bb calibration
# QT PhD Project 2019-2020
# ==========================================

# import

import sys, os, re
import json
from lxml import etree
from decimal import Decimal
# ==========================================
# reading json file
def readInfo(jpath) :

  print 'Reading json file: '+jpath

  assert (os.path.exists(jpath)), 'Cannot find '+jpath+'!'
  
  # init
  inf = {}

  # opening json file
  with open(jpath, 'r') as f:
    # reading json file
    inf = json.load(f)

  return inf


# ==========================================
# creating symbolic links
def mkLink(inpath,outpath,filename):

  if os.path.exists(outpath+'/'+filename) :
    print outpath+'/'+filename+' already exists.'
  else:
    os.symlink(inpath+'/'+filename,outpath+'/'+filename)
    print inpath+'/'+filename+' linked to '+outpath+'/'+filename+'.'

# ==========================================
# printing xml to file
def printToFile(elem,name,outpath) :

  # getting xml string from element
  out_string = etree.tostring(elem,
                              pretty_print=True,
                              doctype='<!DOCTYPE '+name+' SYSTEM \'AnaWSBuilder.dtd\'>')
  
  # opening file
  out_file = open(outpath, 'w')
  
  # printing to file
  out_file.write(out_string)

  print 'Top level card created at: ',outpath

# ==========================================
# adding normalisation factors

def addDefaultNormFactors(tname,reg,samp_elem,data_type,fixt,spurious,fixed,gaus):
# predetermined (default) norm factors

  norm_name = ""
  att = {}

  # ------------ Fixed Templates -----------
  # setting k-factors and normalisations
  if fixed or (fixt and not spurious) :
    # ttbar in SR_confNote
    if (tname == "ttbar") and (data_type == "data") and ("SR_confNote" in reg):
      norm_name = "NormFactor"
      att['Name']='mu_ttbar[0.83]'
    else :
      return

  # ------------- SR/VRqcd -----------------
  if ("SR" in reg) or ("VRqcd" in reg):
    # --- Higgs ---
    if (tname == "Higgs"):
      norm_name = "NormFactor"
      if 'data' in data_type : att['Name']='mu[1,-100,100]'
      else : att['Name']='mu[1,-30,30]'

    # --- Vboson ---
    elif (tname == "Vboson"):
      norm_name = "NormFactor"
      att['Name']='mu_V[1,-20,20]'
    # --- Zboson ---
    elif (tname == "Zboson"):
      norm_name = "NormFactor"
      att['Name']='mu_Z[1,-30,30]'
    # --- Wboson ---
    elif (tname == "Wboson"):
      norm_name = "NormFactor"
      att['Name']='mu_W[1,-20,20]'

  # ------- ttbar in SR_confNote -----------
  if (tname == "ttbar") and (("SR_confNote" in reg) or gaus) :
      norm_name = "Systematic"
      att['Name']='ttbarSF_unc'
      att['Constr']='gaus'
      if (data_type == "data") :
        if "SR_confNote" in reg : att['CentralValue']='0.83'
        else : att['CentralValue']='0.85'
      else :
        att['CentralValue']='1'
      if "SR_confNote" in reg : att['Mag']='0.11'
      else : att['Mag']='0.02'
      att['WhereTo']='yield'
  
  # ---------- ttbar generally ------------
  elif (tname == "ttbar") :
    norm_name = "NormFactor"
    att['Name']='mu_ttbar[1,-50,50]'

  elif (tname == "Wboson") :
    norm_name = "NormFactor"
    att['Name']='mu_Wboson[1,-50,50]'

  elif (tname == "Zboson") :
    norm_name = "NormFactor"
    att['Name']='mu_Zboson[1,-50,50]'  

  elif (tname == "dijets") :
    norm_name = "NormFactor"
    att['Name']='mu_dijets[1,-50,50]'  

  # --------------------
  if norm_name :
    norm = etree.SubElement(samp_elem,norm_name,att)

# ==========================================
# adding yield systematics
# depending on template

def addYieldSyst(tname,reg,samp_elem,att):

  names = ["ATLAS_bEff",
           "ATLAS_cMistag",
           "ATLAS_lMistag"]

  # general attributes
  att['Constr']='gaus'
  att['CentralValue']='1'
  att['WhereTo']='yield'

  # magnitudes
  mags = []

  # Higgs in SR
  if (tname == "Higgs") and ("SR" in reg) :
    mags = ['0.12','0.01','0.01'] # same order as names
  
  # Vboson in SR
  if (tname == "Vboson") and ("SR" in reg) :
    mags = ['0.11','0.03','0.04'] # same order as names
  
  # ttbar in CRttbar 1tag
  if (tname == "ttbar") and ("CRttbar_1tag" in reg) :
    mags = ['0.12','0.01','0.01'] # same order as names

  if mags:
    for i,n in enumerate(names):
      att['Name']=n
      att['Mag']=mags[i]
      syst = etree.SubElement(samp_elem,"Systematic",att)

# ==========================================
# checking if rounded number == 0
# in that case - do not round
# needed for getMagStr() below
# returns val string

def getRound(val,acc):
  s_val = ""
  r_val = round(val,acc)
  if r_val :
    s_val = str(r_val)
  else :
    s_val ='%.2e' % Decimal(str(val))

  return s_val 

# ==========================================
# getting template syst magnitude string 

def getMagStr(syst,t,mag_up,mag_down):

  # convert to string
  s_up = getRound(mag_up,4)
  s_down = getRound(mag_down,4)
        
  
  # var == 0
  if mag_up == 0 : s_up = '1e-10'
  if mag_down == 0 : s_down = '1e-10'
  
  # set string 
  mag_str=s_down+','+s_up
  
  # WARNINGS
  
  # no yield unc.
  if (mag_down == 0) and (mag_up == 0) :
      print 'WARNING: ',t+' '+syst+' has 0 yield magnitude!' 
  
  # opposite effect
  if (mag_down > 0) and (mag_up < 0):
    print 'WARNING: '+t+' ',syst,' up yield eff <0; down yield eff > 0' 

  # both of the same sign
  if (mag_up*mag_down > 0):
    print 'WARNING: '+t+' '+syst+': up and down yield components of the same sign -- up: '+s_up+'; down: '+s_down
  
  return mag_str

# ==========================================
