# ==========================================
# XML Model Creation Helpers:
#   Template Manipulation
#
# Migle Stankaityte
# 16 May 2019
#modified for Z->bb calibration
# ==========================================

#import
import sys, os, re

import numpy as np

import model_helpers as helpers
import syst_smooth 

import ROOT
from ROOT import gROOT, TFile, TTree, TH1F 
gROOT.SetBatch(True)

# ==========================================
# Return maximum difference between normalised 
# histos to check if shape pruning applies
def maxShapeVariation(nom_hist,up_hist,down_hist):
  # do not modify original histograms
  nom_hist_clone  = nom_hist.Clone("nom_hist")
  up_hist_clone   = up_hist.Clone("up_hist")
  down_hist_clone = down_hist.Clone("down_hist")
  # make sure all histograms are normalised to 1
  nom_hist_clone.Scale(1/nom_hist_clone.Integral())
  up_hist_clone.Scale(1/up_hist_clone.Integral())
  down_hist_clone.Scale(1/down_hist_clone.Integral())
  # sys/nom ratios to find maximum bin variation
  up_hist_clone.Divide(nom_hist_clone)
  down_hist_clone.Divide(nom_hist_clone)
  
  maxBinVariation = max(abs(up_hist_clone.GetMaximum()-1),abs(1-up_hist_clone.GetMinimum()),abs(down_hist_clone.GetMaximum()-1),abs(1-down_hist_clone.GetMinimum()))
  #print 'maxBinVariation = '+str(maxBinVariation)

  return maxBinVariation

# end of: maxShapeVariation()
# ==========================================

# ==========================================
# Getting list of bins affected by a gamma parameter for a given template
# returns list of bin numbers
def gammaBins(h,threshold):

  binlist = []
  
  i_int = h.Integral()
  nbins = h.GetNbinsX()

  # looping over bins
  ibin = 1

  while ibin<nbins+1 :

    # getting bin content
    y = h.GetBinContent(ibin)
    yerr = h.GetBinError(ibin)

    # if negative set to infinitecimal ammount
    if y<=0 : 
      x0 = h.GetBinLowEdge(ibin)
      x = x0 + h.GetBinWidth(ibin)
      print 'WARNING!: in ',h.GetTitle(),' bin ',ibin,' value <= 0 in range ',x0,' - ',x,', setting y = 0.00000000001!'
      y = 0.00000000001
      h.SetBinContent(ibin,y)

    elif (yerr/y)>threshold :
        # gamma parameter for 1st bin is gamma_stat_channel_bin_0
        binlist.append(ibin-1)

    ibin += 1

  # rescaling to match the integral of original histogram
  h.Scale(i_int/h.Integral())

  return binlist

# end of: gammaBins()
# ==========================================

# Combining tagging uncertainties flavour-by-flavour summing in quadrature
# returns up and down systematics for each flavour
def getTaggingUncertainties(template,weightlist,tree,nbins,blow,bmax,inpath,truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,region,conf,sepuncerts=False,tag=None) :

  b_uncertainties={}
  if not sepuncerts :
    buncertainties_up={'B':[],'C':[],'Light':[]}
    buncertainties_down={'B':[],'C':[],'Light':[]}
  else :
    buncertainties_up={}
    buncertainties_down={}

  infile = TFile(inpath,"READ")
  histos=[]
  histos2=[]

  # Define all histos for final uncertainties if sepuncerts=False 
  if not sepuncerts :
    histo_Bup=TH1F('Bup','Bup',nbins,blow,bmax)
    histo_Cup=TH1F('Cup','Cup',nbins,blow,bmax)
    histo_Lightup=TH1F('Lightup','Lightup',nbins,blow,bmax)
    histo_Bdown=TH1F('Bdown','Bdown',nbins,blow,bmax)
    histo_Cdown=TH1F('Cdown','Cdown',nbins,blow,bmax)
    histo_Lightdown=TH1F('Lightdown','Lightdown',nbins,blow,bmax)

  for w in weightlist :
    if tag : name = template+w+'_'+tag
    else : name = template+w
    b_uncertainties[name]=TH1F(name,name,nbins,blow,bmax)
    helpers.fillHisto(infile,b_uncertainties[name],truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,tree,region,conf,w)
    histos.append(b_uncertainties[name])
  
  # Storing individual histograms for cross check
  ccpath = 'crosscheck/'
  outfile = ccpath+template+'_btagUncerts'
  outfile2 = ccpath+template+'_btagUncertsCombined'
  helpers.makeHistoRoot(outfile,histos)

  if tag : nominal = b_uncertainties[template+'w_'+tag].Clone('nominal')
  else : nominal = b_uncertainties[template+'w'].Clone('nominal')

  # Do calculation only if sepuncerts=False 
  if not sepuncerts : 
    # Loop over bins
    ibin = 1 
    while ibin < nbins+1 : 
      # Initiate sum of squares variable to start at 0 for each loop
      sum_of_squares_Bup = 0
      sum_of_squares_Bdown = 0
      sum_of_squares_Cup = 0
      sum_of_squares_Cdown = 0
      sum_of_squares_Lightup = 0
      sum_of_squares_Lightdown = 0

      # Loop over histograms in b_uncertainties and calculate sum of squares of each bin minus the nominal
      for name in b_uncertainties : 
        if 'B' in name and '1up' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Bup += bin_content**2

        elif 'B' in name and '1down' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Bdown += bin_content**2

        elif 'C' in name and '1up' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Cup += bin_content**2
          
        elif 'C' in name and '1down' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Cdown += bin_content**2

        elif 'Light' in name and '1up' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Lightup += bin_content**2

        elif 'Light' in name and '1down' in name :
          h = b_uncertainties[name]
          bin_content = h.GetBinContent(ibin)-nominal.GetBinContent(ibin)
          sum_of_squares_Lightdown += bin_content**2


      # Calculate root sum of squares  
      root_sum_of_squares_Bup = sum_of_squares_Bup**0.5
      root_sum_of_squares_Bdown = sum_of_squares_Bdown**0.5
      root_sum_of_squares_Cup = sum_of_squares_Cup**0.5
      root_sum_of_squares_Cdown = sum_of_squares_Cdown**0.5
      root_sum_of_squares_Lightup = sum_of_squares_Lightup**0.5
      root_sum_of_squares_Lightdown = sum_of_squares_Lightdown**0.5

      # Fill new histogram with the root sum of squares of bins of all histograms
      histo_Bup.SetBinContent(ibin,root_sum_of_squares_Bup)
      histo_Bdown.SetBinContent(ibin,root_sum_of_squares_Bdown)
      histo_Cup.SetBinContent(ibin,root_sum_of_squares_Cup)
      histo_Cdown.SetBinContent(ibin,root_sum_of_squares_Cdown)
      histo_Lightup.SetBinContent(ibin,root_sum_of_squares_Lightup)
      histo_Lightdown.SetBinContent(ibin,root_sum_of_squares_Lightdown)

      ibin += 1

    # Add or subtract syst. histograms from nominal as direction requires
    histo_Bup.Scale(-1)
    histo_Bup.Add(nominal)
    histo_Bdown.Add(nominal)
    histo_Cup.Scale(-1)
    histo_Cup.Add(nominal)
    histo_Cdown.Add(nominal)
    histo_Lightup.Scale(-1)
    histo_Lightup.Add(nominal)
    histo_Lightdown.Add(nominal)
  

    # Collect up histograms and down histograms together
    buncertainties_up['B'] = histo_Bup.Clone('B')
    buncertainties_up['C'] = histo_Cup.Clone('C')
    buncertainties_up['Light'] = histo_Lightup.Clone('Light')
    buncertainties_down['B'] = histo_Bdown.Clone('B')
    buncertainties_down['C'] = histo_Cdown.Clone('C')
    buncertainties_down['Light'] = histo_Lightdown.Clone('Light')

    if buncertainties_up['B'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_up['B'])
    if buncertainties_up['C'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_up['C'])
    if buncertainties_up['Light'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_up['Light'])
    if buncertainties_down['B'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_down['B'])
    if buncertainties_down['C'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_down['C'])
    if buncertainties_down['Light'].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_down['Light'])

    # Store combined btagging histos for crosscheck
    histos2.append(histo_Bup)
    histos2.append(histo_Cup)
    histos2.append(histo_Lightup)
    histos2.append(histo_Bdown)
    histos2.append(histo_Cdown)
    histos2.append(histo_Lightdown)
    helpers.makeHistoRoot(outfile2,histos2)

  # If keeping uncertainties separate 
  else : 
    for name in b_uncertainties : 
      if '1up' in name :
        buncertainties_up[name] = b_uncertainties[name].Clone(name)
      if '1down' in name : 
        buncertainties_down[name] = b_uncertainties[name].Clone(name)

    for h in buncertainties_up : 
      if buncertainties_up[h].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_up[h])
    for h in buncertainties_down : 
      if buncertainties_down[h].GetMinimum() <= 0 : helpers.zeroNegBins(buncertainties_down[h])

  return buncertainties_up,buncertainties_down
# end of: getTaggingUncertainties
# ==========================================

# Getting Template Systematic Histograms
# returns nominal histogram and array with sys histograms
#TODO: add choice of variation
#TODO: we can simplify this, as using .SetDirectory(0) will prevent the code from breaking 
     # if histogram is initialised after a file is opened and then closed
     # i.e. we can move histogram creation into fillHisto but it is a bit of work which is not vital

def getTemplateHistos(template, systlist, weightlist, Vjets_wsyst_up, Vjets_wsyst_down, inpath, tree, nbins, blow, bmax, pt_min, pt_max, truth_pt_min, truth_pt_max, w_min, w_max, dsids, region, truth_RW, conf, sepuncerts=False, fit=False, smooth=False, tag=None, basicsymm=False, updown=False):
  
  # if the histogram wil be smoothed set tag
  if fit : tag = 'orig'
  if smooth : tag = 'wide'

  # getting all systematic variations and checking if they exist in file
  allsyst = { 'up' : [], 'down' : []}

  if systlist : 
    for syst in systlist : 
      #allsyst['up'].append(tree+syst+'__1up')
      #allsyst['down'].append(tree+syst+'__1down')
      allsyst['up'].append(syst+'__1up')
      allsyst['down'].append(syst+'__1down')

    exists = helpers.existsInFile(inpath,allsyst['up']+allsyst['down'],region) 
    variations = ['up','down']

  elif weightlist and not smooth :
    allsyst = { 'up' : [], 'down' : []}
    variations = ['up','down']
    print 'No systematics, but alternative b weights are available.'

  else : 
    allsyst = {}
    variations = []
    print 'Nominal only!'

  # --------------------------------
  # initialising histograms
  
  if tag : t_name = template+'_'+tag
  else : t_name = template
  h_nominal = TH1F(t_name,t_name,nbins,blow,bmax)
  
  histos_up = {}
  histos_down = {}
  
  for var in variations : 
    for syst in allsyst[var] :
      
      if tag : name = template+re.sub(tree,'',syst)+'_'+tag
      else : name = template+re.sub(tree,'',syst)

      base = syst.split('__')[0] 
      sname = re.sub(tree,'',base)
      
      assert (exists[base+'__1up']), base+' is not in input file!'
      
      if var == 'up' : histos_up[sname] = TH1F(name,name,nbins,blow,bmax)
      else : histos_down[sname] = TH1F(name,name,nbins,blow,bmax)

  # --------------------------------
  # opening file
  infile = TFile(inpath,"READ")
   # --------------------------------
  # filling existing histograms
  print '\tHistogram input:'
  # getting histograms in input file
  histlist = infile.GetListOfKeys()
  # single histogram input
  print '\t\t Single histogram.'
  if(template == "Zboson") :
    assert(histlist.Contains("Zqq")),'No Zqq histogram in file.'
    h_nominal = infile.Get('Zqq')

  elif(template == "Wboson") :
    assert(histlist.Contains("Wqq")),'No Wqq histogram in file.'
    h_nominal = infile.Get('Wqq')

  elif(template == "dijets") :
    assert(histlist.Contains("dijets")),'No dijets histogram in file.'
    h_nominal = infile.Get('dijets')
  elif(template == "ttbar") :
    assert(histlist.Contains("ttbar")),'No ttbar histogram in file.'
    h_nominal = infile.Get('ttbar')
  elif(template == "ttbar1") :
    assert(histlist.Contains("ttbar")),'No ttbar histogram in file.'
    h_nominal = infile.Get('ttbar')  
  elif(template == "ttbar2") :
    assert(histlist.Contains("ttbar")),'No ttbar histogram in file.'
    h_nominal = infile.Get('ttbar')
  elif(template == "ttbar2W") :
    assert(histlist.Contains("ttbar")),'No ttbar histogram in file.'
    h_nominal = infile.Get('ttbar')  
  elif(template == "ttbar2O") :
    assert(histlist.Contains("ttbar")),'No ttbar histogram in file.'
    h_nominal = infile.Get('ttbar')     
  h_nominal.GetXaxis().SetRangeUser(blow,bmax)
  h_nominal.SetTitle(t_name)
  # --------------------------------
  # filling existing histograms 
  #helpers.fillHisto(infile,h_nominal,truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,tree,region,conf,'w',w_min=w_min,w_max=w_max,truth_RW=truth_RW)

  if not fit and h_nominal.GetMinimum() <= 0 : helpers.zeroNegBins(h_nominal)

  for var in variations : 
    for syst in allsyst[var] :
    
      sname = re.sub(tree,'',(syst.split('__')[0]))
    
      if exists[syst] : 
        
        if var == 'up' : 
          helpers.fillHisto(infile,histos_up[sname],truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,str(syst),region,conf,'w',w_min=w_min,w_max=w_max)
          if not fit and histos_up[sname].GetMinimum() <= 0 : helpers.zeroNegBins(histos_up[sname])
        elif not "JET_MassRes" in sname :
          helpers.fillHisto(infile,histos_down[sname],truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,str(syst),region,conf,'w',w_min=w_min,w_max=w_max)
          if not fit and histos_down[sname].GetMinimum() <= 0 : helpers.zeroNegBins(histos_down[sname])
  
  # --------------------------------
  # return if nominal 
  h_nominal.SetDirectory(0)
  
  if not allsyst and not weightlist : return h_nominal, histos_up, histos_down

  # --------------------------------

  # adding b-tagging systematics histos to histos_up and histos_down
  if weightlist and not smooth :
    buncertainties_up,buncertainties_down = getTaggingUncertainties(template,weightlist,tree,nbins,blow,bmax,inpath,truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,region,conf,sepuncerts,tag)

    # if combining uncertainties, add combined uncertainties to systematics
    if not sepuncerts :
      histos_up['B'] = buncertainties_up['B'].Clone('B')
      histos_up['C'] = buncertainties_up['C'].Clone('C')
      histos_up['Light'] = buncertainties_up['Light'].Clone('Light')
      histos_down['B'] = buncertainties_down['B'].Clone('B')
      histos_down['C'] = buncertainties_down['C'].Clone('C')
      histos_down['Light'] = buncertainties_down['Light'].Clone('Light')

    # if not, add all separate uncertainties with name from weightlist
    else :
      for w in weightlist : 
        if w == 'w' : continue
        else : 
          if tag : name = template+w+'_'+tag
          else : name = template+w
          sname = re.sub(template,'',(name.split('__')[0]))
          if '1up' in name : histos_up[sname] = buncertainties_up[name].Clone(sname)
          elif '1down' in name : histos_down[sname] = buncertainties_down[name].Clone(sname)
        
  # --------------------------------

  # adding weight modelling systematics histos to histos_up and histos_down
  if Vjets_wsyst_up and Vjets_wsyst_down and not smooth :
    print '\nYou defined weight modelling systematics, going to evaluate and include them\n'
    for wmodel_var in Vjets_wsyst_up :
      print wmodel_var
      if tag : name = template+wmodel_var+tag
      else : name = template+wmodel_var
      sname = re.sub(template,'',(name.split('__')[0]))
      histos_up[sname] = TH1F(name,name,nbins,blow,bmax)
      helpers.fillHisto(infile,histos_up[sname],truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,tree,region,conf,wmodel_var)
      if not fit and histos_up[sname].GetMinimum() <= 0 : helpers.zeroNegBins(histos_up[sname])
    for idx,wmodel_var in enumerate(Vjets_wsyst_down) :
      print wmodel_var
      if tag : name = template+Vjets_wsyst_up[idx]+tag
      else : name = template+Vjets_wsyst_up[idx]
      sname = re.sub(template,'',(name.split('__')[0]))
      histos_down[sname] = TH1F(name,name,nbins,blow,bmax)
      helpers.fillHisto(infile,histos_down[sname],truth_pt_min,truth_pt_max,pt_min,pt_max,dsids,tree,region,conf,wmodel_var)
      if not fit and histos_down[sname].GetMinimum() <= 0 : helpers.zeroNegBins(histos_down[sname])

  # --------------------------------
  # closing file
  h_nominal.SetDirectory(0)
  infile.Close()
  print h_nominal.GetTitle()
  # --------------------------------

  # symmetrising onesided histograms
  # currently dummy procedure
  for syst in allsyst['down']: 
    if not exists[syst] or 'JET_MassRes' in syst :
      
      if tag : name = template+re.sub(tree,'',syst)+'_'+tag
      else : name = template+re.sub(tree,'',syst)
      
      base = syst.split('__')[0]
      sname = re.sub(tree,'',base)

      # adding down histos for one-sided if basicsymm=true
      # one-sided symmetrisation
      if basicsymm :
        print '\t This is one-sided, I am going to symmetrise it!\n'
        histos_down[sname] = h_nominal.Clone(name)
        histos_down[sname].Add(histos_up[sname],-1)
        histos_down[sname].Scale(2)
        histos_down[sname].Add(histos_up[sname])

      elif updown : histos_down[sname] = histos_up[sname].Clone(name) 
      else :  histos_down[sname] = h_nominal.Clone(name) 
      histos_down[sname].SetTitle(name)
      histos_down[sname].SetDirectory(0) 
      
      print '\t'+base+' is onesided, '+name+' created by symmetrising.\n'  

  # --------------------------------
  # output
  return h_nominal, histos_up, histos_down

# end of: getTemplateHistos()
# ==========================================
# Adding shape systematics that are templates to main template
# void

def addTemplateSysts(templates,nomtemps,basicsymm=False,updown=False):

  # looping over all templates that have systs. from templates 
  for nt in nomtemps :
    
    print '--------------------------------------'
    print 'Adding shape systematics from other templates to: ',nt
   
    # itterate over templates
    for t in templates:
      if nt == str(t['sample']) :
        
        # iterate over syst-templates
        for st in templates :
          if nt in st['syst_for'] :
            
            systname = str(st['sample'])
            print'\tAdding ',systname,' to ',t['sample'] 
           
            # adding new systematic
            t['syst'].append(systname)
            
            # adding up histos
            t['hsyst_up'][systname] = st['hist'].Clone(systname+'__1up') 
            t['hsyst_up'][systname].SetTitle(systname+'__1up') 
            
            # adding down histos
            # one-sided symmetrisation
            if basicsymm : 
              t['hsyst_down'][systname] = t['hist'].Clone(systname+'__1down')
              t['hsyst_down'][systname].Add(st['hist'],-1)
              t['hsyst_down'][systname].Scale(2)
              t['hsyst_down'][systname].Add(st['hist'])
            elif updown: 
              t['hsyst_down'][systname] = st['hist'].Clone(systname+'__1down')
            else :
              t['hsyst_down'][systname] = t['hist'].Clone(systname+'__1down')
            t['hsyst_down'][systname].SetTitle(systname+'__1down') 

  return

# end of: addTemplateSysts()
# ==========================================
# Making combined templates
# void

def makeCombTemp(templates,combtemps,nominal,nbins,m_min,m_max,nbins_fit,m_fitmin,m_fitmax,ccpath,reg,rebin,sepuncerts=False):

  for ct in combtemps :
    
    print '--------------------------------------'
    print 'Making merged template: ',ct

    # decorrelate JET systematics?
    decorrelate_JET = False

    # intialising containers
    ct_h = TH1F(ct,ct,nbins,m_min,m_max)
    ct_wh = TH1F(ct+'_wide',ct+'_wide',nbins_fit,m_fitmin,m_fitmax)
    ct_syst = []
    ct_h_up = {}
    ct_h_down = {}
    ct_h_up_wide = {}
    ct_h_down_wide = {}
    ct_fix_l = []
    ct_spurious_l = []
    ct_smoothsyst = False
    ct_smoothth = 0
    ct_mu = None

    # b-tagging array to add to syst list
    b_names = ['B','C','Light']
   
    # --------------------------------------
    # looping  
    for t in templates:
      # Save name of v+jets systs
      # syst names saved as 'up', despite weights having different names
      
      if not nominal :
        if 'Vjets_wsyst_up' in t :
          for w in t['Vjets_wsyst_up'] :
            t['syst'].append(w)

        if t['weightlist'] and not sepuncerts : 
          for bsyst in b_names : 
            if bsyst not in t['syst'] : t['syst'].append(bsyst)
        elif t['weightlist'] and sepuncerts : 
          for w in t['weightlist'] : 
            if w == 'w' : continue
            elif w not in t['syst'] : 
              weightname = re.split('__',w)[0]
              t['syst'].append(weightname)

      if ct in t['combine'] :
        # decorrelate JET systematics?
        # yes, if one of the combined templates has decorrelate_JET=true
        if 'decorrelate_JET' in t:
          if t['decorrelate_JET']: decorrelate_JET=True
       
        # adding up nominal
        print '\tAdding sub-template: ',t['sample']
        ct_h.Add(t['hist'])
        if 'whist' in t : ct_wh.Add(t['whist'])
        
        # getting mu
        if 'mu' in t and not ct_mu : ct_mu = t['mu']

        # systematics smoothed?
        if t['smoothsyst'] : 
          ct_smoothsyst = True
          # if variable binning - set threshold to highest given
          if t['smooth_th'] :  
            if ct_smoothth < t['smooth_th'] : ct_smoothth = t['smooth_th']
        
        # making systematic list    
        for syst in t['syst'] :
          
          # for shape syst from other generators, replace with comb-temp name
          if t['sample'] in syst : 
            nsyst = ct+re.sub(t['sample'],'',syst)
            print '\t\tHas shape systematics from other templates! ',syst,' added as: ',nsyst
            syst = nsyst
          # add systematics to the list
          if syst not in ct_syst : ct_syst.append(syst)
        
        # combined template will be fixed / spurious if all combined t's have flag
        ct_fix_l.append(t['fixed'])
        ct_spurious_l.append(t['fixed'])
  

    # --------------------------------------
    # spurious & fixed
    ct_fix = all(ct_fix_l)
    ct_spurious = all(ct_spurious_l)
    
    # --------------------------------------
    # filling systematics
    print '\tSystematics in combined sample:'
    print '\t(in sub-samples where systematic does not exist nominal is added)'
    for s in ct_syst:
      print '\t\t',s
  
      # initialising combined histo
  
      sh_name = ct+s
      if ct in s : sh_name = s
  
      ct_h_up[s] =  TH1F(sh_name+'__1up',sh_name+'__1up',nbins,m_min,m_max)
      ct_h_down[s] =  TH1F(sh_name+'__1down',sh_name+'__1down',nbins,m_min,m_max)
  
      # initialise smoothed histograms
      for t in templates:
        if ct in t['combine']:
          # if syst smooth, initialise
          if s in t['syst'] :
            if 'whsyst_up' in t : 
              if s in t['whsyst_up']: 
                ct_h_up_wide[s] =  TH1F(sh_name+'__1up_wide',sh_name+'__1up_wide',nbins_fit,m_fitmin,m_fitmax)
            if 'whsyst_down' in t : 
              if s in t['whsyst_down']:
                ct_h_down_wide[s] =  TH1F(sh_name+'__1down_wide',sh_name+'__1down_wide',nbins_fit,m_fitmin,m_fitmax)
  
      # looping over templates to be combined
      for t in templates:
        if ct in t['combine']:
          # if syst exists - add
          if s in t['syst'] :
            ct_h_up[s].Add(t['hsyst_up'][s])
            ct_h_down[s].Add(t['hsyst_down'][s])
            if 'whsyst_up' in t : 
              if s in t['whsyst_up']: ct_h_up_wide[s].Add(t['whsyst_up'][s])
            if 'whsyst_down' in t : 
              if s in t['whsyst_down']: ct_h_down_wide[s].Add(t['whsyst_down'][s])
          
          # if shape syst exists - add
          elif (ct in s) and ((t['sample']+re.sub(ct,'',s)) in t['syst']) : 
            ct_h_up[s].Add(t['hsyst_up'][t['sample']+re.sub(ct,'',s)])
            ct_h_down[s].Add(t['hsyst_down'][t['sample']+re.sub(ct,'',s)])
            print '\t\t\tAdded: ',t['sample']+re.sub(ct,'',s),' to ',s
          
          # if syst does not exist - add nominal
          else :
            print '\t\t\t(Not had by ',t['sample'],'! Adding nominal.)'
            ct_h_up[s].Add(t['hist'])
            ct_h_down[s].Add(t['hist'])
            if 'ohist' in t :
              ct_h_up_wide[s].Add(t['whist'])
              ct_h_down_wide[s].Add(t['whist'])
        
    
    # --------------------------------------
    # make template dictionary
    new_temp = {'sample':ct,
               'decorrelate_JET':decorrelate_JET,
               'mu':ct_mu,
               'syst':ct_syst,
               'hist':ct_h,
               'whist':ct_wh,
               'hsyst_up':ct_h_up,
               'hsyst_down':ct_h_down,
               'whsyst_up':ct_h_up_wide,
               'whsyst_down':ct_h_down_wide,
               'fit':False,
               'sub':False,
               'fixed':ct_fix,
               'spurious':ct_spurious,
               'combine':[],
               'injection':None,
               'ignoreWS':False,
               'smoothsyst':ct_smoothsyst,
               'smooth_th':ct_smoothth
               }
    
    # --------------------------------------
    # if systematics smoothed in one of the sub-templates - smooth the whole combination
    if ct_smoothsyst :
      
      # output
      smoothpath = ccpath+'/smoothsysts_'+ct+'/'
      helpers.mkdir(smoothpath)
      
      print "\t\tSystematics smoothed!"
      syst_smooth.getSmoothedSysts(new_temp,rebin,smoothpath,
                                   nbins,m_min,m_max,
                                   nbins_fit,m_fitmin,m_fitmax,
                                   reg)
   
    
    # --------------------------------------
    # add combination to template list
    templates.append(new_temp)


  return

# end of: makeCombTemp()
# ==========================================

# end of file
# ==========================================
# ==========================================
